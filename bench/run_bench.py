#!/usr/bin/env python3
"""Benchmark and correctness harness for fast-mafft (standard library only).

Runs `<prefix>/bin/mafft --thread N data.fa` for every combination of install prefix,
dataset and thread count, times each FFT-NS-2 stage from MAFFT's stderr progress
messages, records wall time and peak RSS (/usr/bin/time -l), and checks that each
candidate's stdout md5 equals the baseline's for the same dataset and thread count.

Every timed run holds the shared timing lock (bench/with_timing_lock.sh). The baseline
(baseline-install) is always re-timed in the same session, interleaved with the
candidates, as PLAN.md requires.

Noise guard: other processes may run unlocked correctness jobs. Inside the lock, each run
first waits (up to --max-bg-wait s) until other processes use at most --max-bg cores,
then samples background CPU every 2 s while MAFFT runs. A run whose mean background
load exceeds --max-bg is marked noisy and retried (up to --retries times); medians use
clean attempts when there are any.

Examples
  # candidate vs baseline on the default suite (sim1k sim2k sim5k prot2k x 1,4,10 threads)
  python3 bench/run_bench.py work/dist/install

  # two candidates, 3 repetitions (median reported), chosen datasets and threads
  python3 bench/run_bench.py dist=work/dist/install work/prog/install -d sim5k sim10k -t 4 10 -r 3

  # correctness only: no lock, baseline not re-run (md5s from bench/expected_md5.json)
  python3 bench/run_bench.py --check work/dist/install -t 10

Outputs go to bench/results/<timestamp>/ (or --out): runs.csv (one row per attempt),
summary.csv (medians), summary.md (tables), logs/ (timestamped stderr per run).
Exit status is 1 if any candidate output differs from the baseline.
"""
import argparse, csv, datetime, hashlib, json, os, platform, re, shlex, statistics
import subprocess, sys, threading, time

BENCH = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(BENCH)
LOCK_SCRIPT = os.path.join(BENCH, "with_timing_lock.sh")
LOCK_DIR = os.path.join(BENCH, ".timing.lock")
EXPECTED = os.path.join(BENCH, "expected_md5.json")
DEFAULT_BASELINE = os.path.join(ROOT, "baseline-install")
DEFAULT_DATASETS = ["sim1k", "sim2k", "sim5k", "prot2k"]
DEFAULT_THREADS = [1, 4, 10]
DEFAULT_MAX_BG = 2.5       # cores; the idle machine sits around 1-2 (system daemons)
DEFAULT_MAX_BG_WAIT = 300  # seconds to wait for quiet before running anyway

# Stage columns, in pipeline order. Each stage runs from its start marker to the next marker.
STAGES = ["startup", "dist6mer", "upgma1", "prog1", "distmsa", "upgma2", "prog2", "finish"]
# Stages shown in the markdown tables (upgma = upgma1 + upgma2, other = startup + finish).
MD_STAGES = [("dist6mer", "6-mer dist"), ("upgma", "UPGMA"), ("prog1", "prog 1/2"),
             ("distmsa", "dist from msa"), ("prog2", "prog 2/2"), ("other", "other")]
PROGRESS_RE = re.compile(r"^\s*(STEP\s+)?\d+\s*/\s*\d+")


# ---------------------------------------------------------------------------------------
# Single measured run (executed inside with_timing_lock.sh via the hidden --_one mode)
# ---------------------------------------------------------------------------------------
def background_cpu(exclude_root=None):
    """Cores busy in processes outside our own tree (sum of ps %cpu / 100). ps %cpu is a
    decaying average, so a job that just started is under-counted for a few seconds."""
    try:
        out = subprocess.run(["ps", "-A", "-o", "pid=,ppid=,%cpu="], capture_output=True,
                             text=True).stdout
        procs = [l.split() for l in out.splitlines() if l.strip()]
        root = exclude_root or os.getpid()
        mine, grew = {root}, True
        while grew:
            grew = False
            for pid, ppid, _ in procs:
                if int(ppid) in mine and int(pid) not in mine:
                    mine.add(int(pid)); grew = True
        return sum(float(c) for pid, _, c in procs if int(pid) not in mine) / 100.0
    except (OSError, ValueError):
        return None


def measure_one(prefix, data, threads, outfile, logfile):
    """Run MAFFT once; return a dict of timings. Stage times come from the arrival time
    of the first byte of each stderr line (MAFFT writes stderr unbuffered)."""
    env = dict(os.environ)
    env.pop("MAFFT_BINARIES", None)  # would silently redirect every prefix to one binary set
    max_bg = float(env.get("FASTMAFFT_MAX_BG", DEFAULT_MAX_BG))
    max_wait = float(env.get("FASTMAFFT_MAX_BG_WAIT", DEFAULT_MAX_BG_WAIT))
    cmd = ["/usr/bin/time", "-l", os.path.join(prefix, "bin", "mafft"), "--thread", str(threads), data]

    bg_pre, waited = background_cpu(), 0.0
    while bg_pre is not None and bg_pre > max_bg and waited < max_wait:
        time.sleep(5); waited += 5
        bg_pre = background_cpu()
    load1 = os.getloadavg()[0]

    samples, stop = [], threading.Event()
    def sampler():
        while not stop.wait(2.0):
            v = background_cpu()
            if v is not None:
                samples.append(v)

    events = []  # (seconds since start, line)
    with open(outfile, "wb") as fo:
        t0 = time.perf_counter()
        p = subprocess.Popen(cmd, stdout=fo, stderr=subprocess.PIPE, env=env)
        th = threading.Thread(target=sampler, daemon=True); th.start()
        fd = p.stderr.fileno()
        buf, line_t = b"", None
        while True:
            chunk = os.read(fd, 65536)
            now = time.perf_counter() - t0
            if not chunk:
                break
            for piece in re.split(rb"([\r\n])", chunk):
                if piece in (b"\r", b"\n"):
                    if buf.strip():
                        events.append((line_t, buf.decode(errors="replace").strip()))
                    buf, line_t = b"", None
                elif piece:
                    if line_t is None:
                        line_t = now
                    buf += piece
        if buf.strip():
            events.append((line_t, buf.decode(errors="replace").strip()))
        rc = p.wait()
        total = time.perf_counter() - t0
        stop.set(); th.join()

    bg_mean = statistics.mean(samples) if samples else bg_pre
    bg_max = max(samples) if samples else bg_pre
    kept = [(t, l) for t, l in events if not PROGRESS_RE.match(l)]
    with open(logfile, "w") as fl:
        fl.write(f"# {shlex.join(cmd)}\n# exit={rc} total={total:.3f}s load1={load1:.2f} "
                 f"bg_pre={bg_pre} bg_mean={bg_mean} bg_max={bg_max} bg_wait={waited:.0f}s\n")
        for t, l in kept:
            fl.write(f"{t:9.3f}  {l}\n")

    res = {"exit_code": rc, "total_s": total, "load1": load1, "bg_cpu": bg_pre,
           "bg_cpu_mean": bg_mean, "bg_cpu_max": bg_max, "bg_wait_s": waited}
    res.update(parse_stages(kept, total))
    res.update(parse_time_l([l for _, l in kept]))
    return res


def parse_stages(events, total):
    """Map stderr marker lines to stage durations (seconds). Missing markers give None."""
    marks = {}
    seen_distmsa = False
    for t, l in events:
        if l.startswith("Making a distance matrix from msa"):
            marks.setdefault("distmsa", t); seen_distmsa = True
        elif l.startswith("Making a distance matrix"):
            marks.setdefault("dist6mer", t)
        elif l.startswith("Constructing a UPGMA tree"):
            marks.setdefault("upgma2" if seen_distmsa else "upgma1", t)
        elif l.startswith("Progressive alignment 1/2"):
            marks.setdefault("prog1", t)
        elif l.startswith("Progressive alignment 2/2"):
            marks.setdefault("prog2", t)
        elif l == "done." and "prog2" in marks:
            marks.setdefault("finish", t)
    order = ["dist6mer", "upgma1", "prog1", "distmsa", "upgma2", "prog2", "finish"]
    out = {f"{s}_s": None for s in STAGES}
    if "dist6mer" in marks:
        out["startup_s"] = marks["dist6mer"]
    for a, b in zip(order, order[1:]):
        if a in marks and b in marks:
            out[f"{a}_s"] = marks[b] - marks[a]
    if "finish" in marks:
        out["finish_s"] = total - marks["finish"]
    return out


def parse_time_l(lines):
    """Parse macOS /usr/bin/time -l output (real/user/sys, maximum resident set size in bytes).
    maxrss covers the whole process tree (the mafft shell script and disttbfast)."""
    out = {"time_real_s": None, "user_s": None, "sys_s": None, "maxrss_mb": None}
    for l in lines:
        m = re.match(r"^([\d.]+)\s+real\s+([\d.]+)\s+user\s+([\d.]+)\s+sys$", l)
        if m:
            out["time_real_s"], out["user_s"], out["sys_s"] = map(float, m.groups())
        m = re.match(r"^(\d+)\s+maximum resident set size$", l)
        if m:
            out["maxrss_mb"] = int(m.group(1)) / 2**20
    return out


# ---------------------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------------------
def md5_file(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def resolve_prefixes(args_prefixes, baseline_path):
    prefixes, labels = [], set()
    base_real = os.path.realpath(baseline_path)
    items = [("baseline", baseline_path)]
    for spec in args_prefixes:
        label, _, path = spec.partition("=") if "=" in spec else ("", "", spec)
        path = os.path.abspath(path)
        if not label:
            if os.path.realpath(path) == base_real:
                label = "baseline"
            elif os.path.basename(path.rstrip("/")) == "install":
                label = os.path.basename(os.path.dirname(path.rstrip("/")))
            else:
                label = os.path.basename(path.rstrip("/"))
        items.append((label, path))
    for label, path in items:
        path = os.path.abspath(path)
        if any(os.path.realpath(p["path"]) == os.path.realpath(path) for p in prefixes):
            continue  # same install listed twice (e.g. baseline passed explicitly)
        mafft = os.path.join(path, "bin", "mafft")
        if not os.access(mafft, os.X_OK):
            sys.exit(f"error: {mafft} is not an executable")
        check_script_prefix(mafft, path)
        base, k = label, 2
        while label in labels:
            label, k = f"{base}{k}", k + 1
        labels.add(label)
        prefixes.append({"label": label, "path": path,
                         "is_baseline": os.path.realpath(path) == base_real})
    return prefixes


def check_script_prefix(mafft, path):
    """The mafft driver script hard-codes its libexec directory. A script copied from another
    install would silently run that install's binaries, so warn when it points elsewhere."""
    try:
        with open(mafft, errors="replace") as f:
            for line in f:
                m = re.match(r"\s*prefix=(\S+libexec/mafft)\s*$", line)
                if m:
                    lib = os.path.realpath(m.group(1).strip("\"'"))
                    if not lib.startswith(os.path.realpath(path) + os.sep):
                        print(f"WARNING: {mafft} runs binaries from {lib}, outside its own prefix",
                              file=sys.stderr)
                    return
    except OSError:
        pass


def resolve_datasets(names):
    out = []
    for n in names:
        path = n if os.path.isfile(n) else os.path.join(BENCH, "data", n + ".fa")
        if not os.path.isfile(path):
            sys.exit(f"error: dataset {n} not found (looked for {path})")
        out.append({"name": re.sub(r"\.(fa|fasta|fas)$", "", os.path.basename(path)),
                    "path": os.path.abspath(path), "md5": md5_file(path)})
    return out


def load_expected():
    try:
        with open(EXPECTED) as f:
            return json.load(f)
    except (OSError, ValueError):
        return {}


def save_expected(new_entries):
    """Add baseline md5s to bench/expected_md5.json (never overwrite; warn on conflict)."""
    exp = load_expected()
    changed = False
    for (ds, thr), (in_md5, out_md5) in new_entries.items():
        e = exp.setdefault(ds, {"input_md5": in_md5, "output_md5": {}})
        if e.get("input_md5") != in_md5:
            print(f"WARNING: {ds}: input md5 {in_md5} differs from recorded {e.get('input_md5')};"
                  " not updating expected_md5.json", file=sys.stderr)
            continue
        old = e["output_md5"].get(str(thr))
        if old is None:
            e["output_md5"][str(thr)] = out_md5; changed = True
        elif old != out_md5:
            print(f"WARNING: baseline {ds} t={thr} md5 {out_md5} != recorded {old}", file=sys.stderr)
    if changed:
        tmp = EXPECTED + f".tmp{os.getpid()}"
        with open(tmp, "w") as f:
            json.dump(exp, f, indent=1, sort_keys=True)
            f.write("\n")
        os.replace(tmp, EXPECTED)


def run_one(prefix, ds, threads, tag, outdir, locked, env):
    alnfile = os.path.join(outdir, "aln", tag + ".fa")
    logfile = os.path.join(outdir, "logs", tag + ".log")
    inner = [sys.executable, os.path.abspath(__file__), "--_one", prefix["path"], ds["path"],
             str(threads), alnfile, logfile]
    cmd = [LOCK_SCRIPT] + inner if locked else inner
    if locked and os.path.isdir(LOCK_DIR):
        print("    (waiting for the timing lock held by another run...)", file=sys.stderr, flush=True)
    t0 = time.perf_counter()
    p = subprocess.run(cmd, stdout=subprocess.PIPE, text=True, env=env)
    elapsed = time.perf_counter() - t0
    if p.returncode != 0 or not p.stdout.strip():
        sys.exit(f"error: measurement failed for {tag} (exit {p.returncode}); see {logfile}")
    res = json.loads(p.stdout.strip().splitlines()[-1])
    res["lock_wait_s"] = (max(0.0, elapsed - res["total_s"] - res.get("bg_wait_s", 0) - 0.1)
                          if locked else 0.0)
    res["md5"] = md5_file(alnfile)
    res["aln"] = alnfile
    return res


def fmt(x, nd=2):
    return "" if x is None else f"{x:.{nd}f}"


def median(vals):
    vals = [v for v in vals if v is not None]
    return statistics.median(vals) if vals else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefixes", nargs="*", help="install prefixes to test, optionally LABEL=PATH")
    ap.add_argument("-d", "--datasets", nargs="+", default=DEFAULT_DATASETS,
                    help="names in bench/data (without .fa) or FASTA paths (default: %(default)s)")
    ap.add_argument("-t", "--threads", nargs="+", type=int, default=DEFAULT_THREADS,
                    help="thread counts (default: %(default)s)")
    ap.add_argument("-r", "--reps", type=int, default=1, help="repetitions per combination; medians reported")
    ap.add_argument("-o", "--out", help="output directory (default: bench/results/<timestamp>)")
    ap.add_argument("--baseline", default=DEFAULT_BASELINE, help="baseline prefix (default: %(default)s)")
    ap.add_argument("--check", action="store_true",
                    help="correctness only: no timing lock, no noise guard, baseline re-run only "
                         "where bench/expected_md5.json has no entry; times are not comparable")
    ap.add_argument("--max-bg", type=float, default=DEFAULT_MAX_BG,
                    help="background cores above which a run counts as noisy (default: %(default)s)")
    ap.add_argument("--max-bg-wait", type=float, default=DEFAULT_MAX_BG_WAIT,
                    help="seconds to wait for a quiet machine before each run (default: %(default)s)")
    ap.add_argument("--retries", type=int, default=2,
                    help="extra attempts for a noisy run (default: %(default)s)")
    ap.add_argument("--keep-outputs", action="store_true",
                    help="keep every alignment (default: keep only outputs that differ from baseline)")
    ap.add_argument("--_one", nargs=5, help=argparse.SUPPRESS)
    a = ap.parse_args()

    if a._one:
        prefix, data, threads, outfile, logfile = a._one
        print(json.dumps(measure_one(prefix, data, int(threads), outfile, logfile)))
        return

    locked = not a.check
    env = dict(os.environ, FASTMAFFT_MAX_BG=str(a.max_bg),
               FASTMAFFT_MAX_BG_WAIT=str(a.max_bg_wait if locked else 0))
    retries = a.retries if locked else 0
    datasets = resolve_datasets(a.datasets)
    expected = load_expected()
    prefixes = resolve_prefixes(a.prefixes, a.baseline)
    if a.check:
        # Only re-run the baseline where no recorded md5 exists for some combination.
        need_base = any(expected.get(ds["name"], {}).get("input_md5") != ds["md5"] or
                        str(t) not in expected.get(ds["name"], {}).get("output_md5", {})
                        for ds in datasets for t in a.threads)
        if not need_base and len(prefixes) > 1:
            prefixes = [p for p in prefixes if not p["is_baseline"]]
    base_labels = {p["label"] for p in prefixes if p["is_baseline"]}

    stamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    outdir = os.path.abspath(a.out or os.path.join(BENCH, "results", stamp + ("-check" if a.check else "")))
    os.makedirs(os.path.join(outdir, "aln"), exist_ok=True)
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)

    print(f"results -> {outdir}", file=sys.stderr)
    for p in prefixes:
        print(f"  {p['label']:<12} {p['path']}", file=sys.stderr)

    rows = []
    n_total = len(datasets) * len(a.threads) * a.reps * len(prefixes)
    k = 0
    for ds in datasets:
        for thr in a.threads:
            for rep in range(1, a.reps + 1):
                # Alternate the order each repetition so drift does not favour one prefix.
                order = prefixes if rep % 2 == 1 else list(reversed(prefixes))
                for p in order:
                    k += 1
                    attempts = []
                    for att in range(1, retries + 2):
                        tag = f"{p['label']}.{ds['name']}.t{thr}.r{rep}" + (f".a{att}" if att > 1 else "")
                        r = run_one(p, ds, thr, tag, outdir, locked, env)
                        noisy = locked and r.get("bg_cpu_mean") is not None and r["bg_cpu_mean"] > a.max_bg
                        row = {"label": p["label"], "prefix": p["path"], "dataset": ds["name"],
                               "threads": thr, "rep": rep, "attempt": att, "locked": int(locked),
                               "noisy": int(noisy), "used": 0, **r}
                        attempts.append(row)
                        print(f"[{k}/{n_total}] {p['label']:<10} {ds['name']:<10} t={thr:<3} rep{rep}"
                              + (f" try{att}" if att > 1 else "") +
                              f"  total {r['total_s']:7.2f}s  rss {fmt(r['maxrss_mb'], 0):>5} MB  "
                              f"md5 {r['md5'][:12]}  bg {fmt(r.get('bg_cpu_mean'), 1)} cores"
                              + ("  NOISY" if noisy else "")
                              + (f"  lock wait {r['lock_wait_s']:.0f}s" if r['lock_wait_s'] > 1 else "")
                              + (f"  quiet wait {r['bg_wait_s']:.0f}s" if r.get('bg_wait_s') else ""),
                              file=sys.stderr, flush=True)
                        if r["exit_code"] != 0:
                            print(f"WARNING: exit code {r['exit_code']}", file=sys.stderr)
                        if not noisy:
                            break
                    clean = [x for x in attempts if not x["noisy"]]
                    chosen = clean[0] if clean else min(attempts, key=lambda x: x["bg_cpu_mean"] or 0)
                    chosen["used"] = 1
                    rows.extend(attempts)

    # Reference md5 per (dataset, threads): this session's baseline, else expected_md5.json.
    ref, new_expected = {}, {}
    for ds in datasets:
        for thr in a.threads:
            b = {r["md5"] for r in rows if r["dataset"] == ds["name"] and r["threads"] == thr
                 and r["label"] in base_labels and r["exit_code"] == 0}
            if len(b) > 1:
                print(f"WARNING: baseline output is not deterministic for {ds['name']} t={thr}: {b}",
                      file=sys.stderr)
            if len(b) == 1:
                ref[(ds["name"], thr)] = b.pop()
                if os.path.dirname(ds["path"]) == os.path.join(BENCH, "data"):  # shared sets only
                    new_expected[(ds["name"], thr)] = (ds["md5"], ref[(ds["name"], thr)])
            else:
                e = expected.get(ds["name"], {})
                if e.get("input_md5") == ds["md5"] and str(thr) in e.get("output_md5", {}):
                    ref[(ds["name"], thr)] = e["output_md5"][str(thr)]
    if new_expected:
        save_expected(new_expected)
    exp_now = load_expected()
    for (dsn, thr), m in ref.items():
        rec = exp_now.get(dsn, {}).get("output_md5", {}).get(str(thr))
        if rec and rec != m:
            print(f"WARNING: {dsn} t={thr}: session baseline md5 {m} != recorded {rec}", file=sys.stderr)

    for r in rows:
        m = ref.get((r["dataset"], r["threads"]))
        r["identical"] = "NO-REF" if m is None else ("IDENTICAL" if r["md5"] == m else "DIFFERENT")
        if not a.keep_outputs and (r["identical"] != "DIFFERENT" or r["label"] in base_labels):
            try:
                os.remove(r["aln"])
            except OSError:
                pass

    write_outputs(rows, prefixes, datasets, a, outdir, locked)
    n_diff = sum(r["identical"] == "DIFFERENT" for r in rows)
    n_noisy = sum(r["noisy"] and r["used"] for r in rows)
    print(f"\n{'ALL IDENTICAL' if n_diff == 0 else f'{n_diff} RUN(S) DIFFERENT'}"
          + (f"; {n_noisy} reported run(s) still noisy after retries" if n_noisy else "")
          + f"; wrote {outdir}/summary.md", file=sys.stderr)
    sys.exit(1 if n_diff else 0)


def write_outputs(rows, prefixes, datasets, a, outdir, locked):
    run_cols = ["label", "dataset", "threads", "rep", "attempt", "used", "noisy", "locked",
                "lock_wait_s", "bg_wait_s", "bg_cpu", "bg_cpu_mean", "bg_cpu_max", "load1",
                "exit_code", "total_s"] + [f"{s}_s" for s in STAGES] + \
               ["time_real_s", "user_s", "sys_s", "maxrss_mb", "md5", "identical", "prefix"]
    with open(os.path.join(outdir, "runs.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=run_cols, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({c: (f"{r[c]:.3f}" if isinstance(r.get(c), float) else r.get(c)) for c in run_cols})

    summ = {}  # medians over the used attempt of each rep; ordered by dataset, threads, prefix
    for ds in datasets:
        for thr in a.threads:
            for p in prefixes:
                every = [r for r in rows if r["label"] == p["label"] and r["dataset"] == ds["name"]
                         and r["threads"] == thr]
                rs = [r for r in every if r["used"]]
                if not rs:
                    continue
                s = {c: median([r.get(c) for r in rs]) for c in
                     ["total_s", "user_s", "maxrss_mb", "bg_cpu_mean"] + [f"{x}_s" for x in STAGES]}
                s["upgma_s"] = (None if s["upgma1_s"] is None or s["upgma2_s"] is None
                                else s["upgma1_s"] + s["upgma2_s"])
                s["other_s"] = (None if s["startup_s"] is None or s["finish_s"] is None
                                else s["startup_s"] + s["finish_s"])
                ids = {r["identical"] for r in every}  # any differing attempt is a failure
                s["identical"] = "DIFFERENT" if "DIFFERENT" in ids else ids.pop() if len(ids) == 1 else "MIXED"
                s["md5"] = ",".join(sorted({r["md5"] for r in every}))
                s["n"] = len(rs)
                s["noisy"] = sum(r["noisy"] for r in rs)
                summ[(p["label"], ds["name"], thr)] = s
    base = next((p["label"] for p in prefixes if p["is_baseline"]), None)

    def speedup(bv, cv):
        return None if bv is None or cv is None or cv <= 0 else bv / cv

    sum_cols = ["label", "dataset", "threads", "n_reps", "noisy_reps", "total_s", "speedup_total"] + \
               [f"{s}_s" for s, _ in MD_STAGES] + [f"speedup_{s}" for s, _ in MD_STAGES] + \
               ["user_s", "maxrss_mb", "bg_cpu_mean", "identical", "md5"]
    with open(os.path.join(outdir, "summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(sum_cols)
        for (lab, dsn, thr), s in summ.items():
            b = summ.get((base, dsn, thr)) if base else None
            row = [lab, dsn, thr, s["n"], s["noisy"], fmt(s["total_s"], 3),
                   fmt(speedup(b["total_s"], s["total_s"]), 3) if b else ""]
            row += [fmt(s[f"{x}_s"], 3) for x, _ in MD_STAGES]
            row += [fmt(speedup(b[f"{x}_s"], s[f"{x}_s"]), 3) if b else "" for x, _ in MD_STAGES]
            row += [fmt(s["user_s"], 2), fmt(s["maxrss_mb"], 1), fmt(s["bg_cpu_mean"], 2),
                    s["identical"], s["md5"]]
            w.writerow(row)

    L = []
    L.append(f"# MAFFT benchmark {datetime.datetime.now():%Y-%m-%d %H:%M}\n")
    L.append(f"- Command: `{shlex.join([os.path.relpath(os.path.abspath(sys.argv[0]), ROOT)] + sys.argv[1:])}`")
    L.append(f"- Host: {platform.node()} ({platform.machine()}, {os.cpu_count()} cores), "
             f"Python {platform.python_version()}")
    L.append(f"- Repetitions: {a.reps} (medians shown); timing lock: "
             + ("held for every run" if locked else "NOT held (--check): times not comparable"))
    for p in prefixes:
        L.append(f"- `{p['label']}`: `{p['path']}`")
    L.append("- Datasets: " + ", ".join(f"{d['name']} (input md5 {d['md5'][:12]})" for d in datasets))
    L.append("\nStage times in seconds, from MAFFT stderr markers. UPGMA = both guide trees; "
             "other = startup before the 6-mer stage plus output after progressive 2/2. "
             "user = user CPU seconds. bg = mean cores used by other processes during the run "
             f"(system daemons alone give about 1-2); runs above {a.max_bg} are counted as noisy.\n")

    L.append("## Absolute times (median)\n")
    L.append("| dataset | threads | build | " + " | ".join(n for _, n in MD_STAGES) +
             " | total | user | peak RSS MB | bg | noisy | output md5 | vs baseline |")
    L.append("|---|---:|---|" + "---:|" * (len(MD_STAGES) + 5) + "---|---|")
    for (lab, dsn, thr), s in summ.items():
        L.append(f"| {dsn} | {thr} | {lab} | " + " | ".join(fmt(s[f'{x}_s']) for x, _ in MD_STAGES) +
                 f" | **{fmt(s['total_s'])}** | {fmt(s['user_s'], 1)} | {fmt(s['maxrss_mb'], 0)} | "
                 f"{fmt(s['bg_cpu_mean'], 1)} | {s['noisy']}/{s['n']} | `{s['md5'][:12]}` | {s['identical']} |")

    cands = [p["label"] for p in prefixes if not p["is_baseline"]]
    if base and cands:
        L.append("\n## Speedup vs baseline (baseline median / candidate median; >1 is faster)\n")
        L.append("| dataset | threads | build | " + " | ".join(n for _, n in MD_STAGES) +
                 " | base total s | cand total s | **total speedup** | RSS ratio | output |")
        L.append("|---|---:|---|" + "---:|" * (len(MD_STAGES) + 4) + "---|")
        for (lab, dsn, thr), s in summ.items():
            if lab not in cands or (base, dsn, thr) not in summ:
                continue
            b = summ[(base, dsn, thr)]
            sp = [speedup(b[f"{x}_s"], s[f"{x}_s"]) for x, _ in MD_STAGES]
            rss = None if not b["maxrss_mb"] or s["maxrss_mb"] is None else s["maxrss_mb"] / b["maxrss_mb"]
            L.append(f"| {dsn} | {thr} | {lab} | " + " | ".join(fmt(v) + ("x" if v else "") for v in sp) +
                     f" | {fmt(b['total_s'])} | {fmt(s['total_s'])} | "
                     f"**{fmt(speedup(b['total_s'], s['total_s']))}x** | {fmt(rss)} | {s['identical']} |")
    elif cands:
        L.append("\n## Correctness vs recorded baseline md5 (bench/expected_md5.json)\n")
        L.append("| dataset | threads | build | output |")
        L.append("|---|---:|---|---|")
        for (lab, dsn, thr), s in summ.items():
            L.append(f"| {dsn} | {thr} | {lab} | {s['identical']} |")
    with open(os.path.join(outdir, "summary.md"), "w") as f:
        f.write("\n".join(L) + "\n")


if __name__ == "__main__":
    main()
