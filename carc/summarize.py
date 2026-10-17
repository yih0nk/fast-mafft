#!/usr/bin/env python3
# Python 3.6-compatible (CARC). Usage: summarize.py JOB1.log JOB2.log ...
# Reads the result lines that carc_run.py prints inside carc_v9_final.sh / carc_v10_final.sh jobs:
#   LABEL dataset t=N wall=S maxrss_gb=G md5=XXXXXXXXXXXX exit=E stages=...
# with LABEL = baseAggr-rK (baseline) or v9aggr-rK / v10aggr-rK (fast-mafft), K = repetition number.
# Prints the two tables of RESULTS.md (summary and raw runs) as Markdown, plus any reruns after an md5
# difference (baseAgain / v10Again lines).
# Per-run speedup = baseline wall / fast-mafft wall of the same repetition (same job, same node).
# Mean speedup = arithmetic mean of the per-run speedups; min = the smallest of them.
import re, sys

ORDER = ["sim1k", "rnasim1k", "prot2k", "sim2k", "rnasim2k",
         "rnasim5k", "rnasim10k", "16S.T", "16S.3", "homfam_p450", "rnasim20k", "16S.B.ALL",
         "homfam_sdr", "homfam_Acetyltransf", "rnasim25k", "homfam_rvp", "rnasim50k",
         "sarscov2_500", "sarscov2_1000", "sarscov2_2000", "edge", "long100x30k", "long1000x30k",
         "long300x50k", "long500x30k"]
LINE = re.compile(r"^(baseAggr|v\d+aggr|baseAgain|v\d+Again)-r(\d+) (\S+) t=(\d+) wall=([0-9.]+) maxrss_gb=([0-9.]+) "
                  r"md5=([0-9a-f]+) exit=(-?\d+)")

runs = {}  # (dataset, rep) -> {"base": (wall, rss, md5, exit), "v9": ...}
for path in sys.argv[1:]:
    for line in open(path, errors="replace"):
        m = LINE.match(line.strip())
        if not m:
            continue
        lab = m.group(1)
        who = "base" if lab == "baseAggr" else "again-base" if lab == "baseAgain" else "again-v" if lab.endswith("Again") else "v9"
        key = (m.group(3), int(m.group(2)))
        runs.setdefault(key, {})[who] = (float(m.group(5)), float(m.group(6)), m.group(7), int(m.group(8)))

reps = sorted({r for (_, r) in runs})
sets = [d for d in ORDER if any((d, r) in runs for r in reps)]
sets += sorted({d for (d, _) in runs} - set(sets))

def fmt(x, nd=1):
    return "-" if x is None else ("%.*f" % (nd, x))

print("| dataset | baseline mean (s) | fast-mafft mean (s) | " + " | ".join("speedup r%d" % r for r in reps) +
      " | mean speedup | min speedup | md5 match | peak RSS baseline / fast-mafft (GB) |")
print("|---|---:|---:|" + "---:|" * len(reps) + "---:|---:|---|---|")
for d in sets:
    bw, vw, sp, ok, brss, vrss = [], [], [], True, [], []
    for r in reps:
        x = runs.get((d, r), {})
        b, v = x.get("base"), x.get("v9")
        if b: bw.append(b[0]); brss.append(b[1])
        if v: vw.append(v[0]); vrss.append(v[1])
        if b and v and v[0] > 0:
            sp.append(b[0] / v[0])
            if b[2] != v[2] or b[3] != 0 or v[3] != 0:
                ok = False
        else:
            sp.append(None); ok = False
    good = [s for s in sp if s is not None]
    md5 = ("yes (%d/%d)" % (len(good), len(reps))) if ok else "**NO / incomplete**"
    print("| %s | %s | %s | %s | %s | %s | %s | %s / %s |" % (
        d, fmt(sum(bw) / len(bw) if bw else None), fmt(sum(vw) / len(vw) if vw else None),
        " | ".join(fmt(s, 2) + ("x" if s is not None else "") for s in sp),
        fmt(sum(good) / len(good) if good else None, 2) + ("x" if good else ""),
        fmt(min(good) if good else None, 2) + ("x" if good else ""),
        md5, fmt(max(brss) if brss else None, 2), fmt(max(vrss) if vrss else None, 2)))

print()
print("| dataset | rep | baseline wall (s) | fast-mafft wall (s) | speedup | baseline md5 | fast-mafft md5 | baseline RSS (GB) | fast-mafft RSS (GB) |")
print("|---|---:|---:|---:|---:|---|---|---:|---:|")
for d in sets:
    for r in reps:
        x = runs.get((d, r), {})
        b, v = x.get("base"), x.get("v9")
        s = (b[0] / v[0]) if (b and v and v[0] > 0) else None
        print("| %s | %d | %s | %s | %s | %s | %s | %s | %s |" % (
            d, r, fmt(b[0] if b else None), fmt(v[0] if v else None),
            fmt(s, 2) + ("x" if s is not None else ""), b[2] if b else "-", v[2] if v else "-",
            fmt(b[1] if b else None, 2), fmt(v[1] if v else None, 2)))

again = [(d, r, x) for (d, r), x in sorted(runs.items()) if "again-base" in x or "again-v" in x]
if again:
    print()
    print("Reruns after an md5 difference (does MAFFT give the same output twice?):")
    print()
    print("| dataset | rep | baseline md5, run 1 / run 2 | fast-mafft md5, run 1 / run 2 |")
    print("|---|---:|---|---|")
    for d, r, x in again:
        b1, b2, v1, v2 = x.get("base"), x.get("again-base"), x.get("v9"), x.get("again-v")
        print("| %s | %d | %s / %s | %s / %s |" % (d, r, b1[2] if b1 else "-", b2[2] if b2 else "-", v1[2] if v1 else "-", v2[2] if v2 else "-"))
