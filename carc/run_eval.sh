#!/bin/bash
#SBATCH --job-name=fast-mafft-eval
#SBATCH --partition=main
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=fast-mafft-eval-%j.log
#
# Evaluate fast-mafft against MAFFT v7.526, both built from this repository on the node that runs
# the job. For every input file it runs MAFFT and fast-mafft back to back at --thread N, checks that
# the two outputs are byte-identical, and reports wall time, speedup and peak memory.
#
# Usage (submit from the repository root):
#   sbatch --account=<account> carc/run_eval.sh                          # bench/data (bench/get_data.sh)
#   sbatch --account=<account> carc/run_eval.sh --data /path/to/fastas --reps 3
#   bash carc/run_eval.sh --data DIR --threads 8                         # without Slurm
#
# Options:
#   --data DIR       input directory, searched recursively for *.fa *.fasta *.fas *.fna *.faa *.tfa
#                    (default: bench/data). Duplicate files are run once.
#   --threads N      --thread value for both programs (default: $SLURM_CPUS_PER_TASK, else 32)
#   --reps R         repetitions per input, interleaved (default: 1)
#   --flags "..."    CFLAGS for both builds (default: "-O3 -march=native -flto")
#   --check-default  also build both with the Makefile's default flags, and confirm that the
#                    default-flag outputs are byte-identical too (one extra run of each per input)
#   --selfcheck      also run fast-mafft once per input with its internal self-checks on
#                    (every 6-mer distance, MSA distance and guide-tree step recomputed the original way)
#   --work DIR       build and output directory (default: <repo>/eval-work/<job id or date>)
set -uo pipefail

REPO=${FASTMAFFT_REPO:-${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}}
REPO=$(cd "$REPO" && pwd)
[ -d "$REPO/fast-mafft-src/core" ] || { echo "error: run from the fast-mafft repository root (no $REPO/fast-mafft-src)"; exit 2; }
DATA=$REPO/bench/data; T=${SLURM_CPUS_PER_TASK:-32}; REPS=1; FLAGS="-O3 -march=native -flto"
CHECKDEF=0; SELFCHECK=0; WORK=""
while [ $# -gt 0 ]; do
  case $1 in
    --data) DATA=$2; shift 2;; --threads) T=$2; shift 2;; --reps) REPS=$2; shift 2;;
    --flags) FLAGS=$2; shift 2;; --check-default) CHECKDEF=1; shift;; --selfcheck) SELFCHECK=1; shift;;
    --work) WORK=$2; shift 2;; -h|--help) sed -n '/^# Evaluate fast-mafft/,/^set -uo/p' "$0" | sed -e '$d' -e 's/^#//' -e 's/^ //'; exit 0;;
    *) echo "unknown option $1"; exit 2;;
  esac
done
DATA=$(cd "$DATA" 2>/dev/null && pwd) || { echo "error: no data directory (run bench/get_data.sh, or pass --data DIR)"; exit 2; }
WORK=${WORK:-$REPO/eval-work/${SLURM_JOB_ID:-$(date +%Y%m%d-%H%M%S)}}
mkdir -p "$WORK"
if [ -f /etc/profile.d/modules.sh ]; then source /etc/profile.d/modules.sh; module load gcc/13.3.0 >/dev/null 2>&1 || true; fi
CC=${CC:-gcc}
RUN="python3 $REPO/bench/carc_run.py"

echo "== fast-mafft evaluation  $(date)"
echo "host $(hostname)  cpu: $( (lscpu 2>/dev/null | grep 'Model name' | sed 's/.*: *//') || sysctl -n machdep.cpu.brand_string 2>/dev/null)"
echo "compiler: $($CC --version | head -1)  CFLAGS: $FLAGS  threads: $T  reps: $REPS"
echo "repo: $REPO ($(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo 'no git'))  data: $DATA  work: $WORK"

# ---- sources: MAFFT v7.526 = this repository's first commit; fast-mafft = fast-mafft-src/
SRC_BASE=$WORK/mafft-7.526; SRC_FAST=$WORK/fast-mafft-src
rm -rf "$SRC_BASE" "$SRC_FAST" "$WORK/archive"; cp -r "$REPO/fast-mafft-src" "$SRC_FAST"
if ROOT=$(git -C "$REPO" rev-list --max-parents=0 HEAD 2>/dev/null); then
  mkdir -p "$WORK/archive" && git -C "$REPO" archive "$ROOT" fast-mafft-src | tar -x -C "$WORK/archive" \
    && mv "$WORK/archive/fast-mafft-src" "$SRC_BASE" && rm -rf "$WORK/archive"
else
  echo "no git history: recreating MAFFT v7.526 by reversing fast-mafft-v9.patch"
  cp -r "$REPO/fast-mafft-src" "$SRC_BASE" && (cd "$SRC_BASE" && patch -R -p1 -s < "$REPO/fast-mafft-v9.patch")
fi
md5f() { python3 -c "import hashlib,sys;print(hashlib.md5(open(sys.argv[1],'rb').read()).hexdigest())" "$1"; }
[ "$(md5f "$SRC_BASE/core/disttbfast.c")" = 564c5519e44b43d888f14b8badb23032 ] && [ "$(md5f "$SRC_BASE/core/mltaln9.c")" = c4c2055eb66c16905ef3c6090d74e580 ] \
  || { echo "error: baseline source is not pristine MAFFT v7.526"; exit 3; }
echo "baseline source verified as pristine MAFFT v7.526 (sysimm/mafft ee97999)"

build() { # name src [cflags]: stock Makefile, same PREFIX on both make calls
  local name=$1 src=$2 flags=${3:-}
  make -C "$src/core" clean >/dev/null 2>&1
  if [ -n "$flags" ]; then
    make -j"$T" -C "$src/core" CC="$CC" CFLAGS="$flags" PREFIX="$WORK/$name" >"$WORK/build-$name.log" 2>&1 \
      && make -C "$src/core" CC="$CC" CFLAGS="$flags" PREFIX="$WORK/$name" install >>"$WORK/build-$name.log" 2>&1
  else
    make -j"$T" -C "$src/core" CC="$CC" PREFIX="$WORK/$name" >"$WORK/build-$name.log" 2>&1 \
      && make -C "$src/core" CC="$CC" PREFIX="$WORK/$name" install >>"$WORK/build-$name.log" 2>&1
  fi || { echo "error: build of $name failed, see $WORK/build-$name.log"; exit 4; }
  echo "built $name"
}
cp -r "$SRC_BASE" "$WORK/src-base-def"; cp -r "$SRC_FAST" "$WORK/src-fast-def"
build mafft "$SRC_BASE" "$FLAGS"; build fast "$SRC_FAST" "$FLAGS"
if [ $CHECKDEF = 1 ]; then build mafft-default "$WORK/src-base-def"; build fast-default "$WORK/src-fast-def"; fi

# ---- inputs (sorted by size, duplicates by content skipped)
INPUTS=(); while IFS= read -r line; do INPUTS+=("$line"); done < <(python3 - "$DATA" <<'EOF'
import hashlib, os, sys
root = sys.argv[1]; seen = set(); files = []
for d, _, fs in os.walk(root):
    for f in fs:
        if f.lower().endswith((".fa", ".fasta", ".fas", ".fna", ".faa", ".tfa")):
            files.append(os.path.join(d, f))
for p in sorted(files, key=lambda p: (os.path.getsize(p), p)):
    h = hashlib.md5(open(p, "rb").read()).hexdigest()
    if h not in seen:
        seen.add(h); print(p)
EOF
)
echo "inputs: ${#INPUTS[@]}"
RES=$WORK/runs.txt; : > "$RES"
for f in "${INPUTS[@]}"; do
  rel=${f#$DATA/}
  for r in $(seq 1 "$REPS"); do
    $RUN "$WORK/mafft/bin/mafft" "$T" "$f" "mafft|$r|$rel" | tee -a "$RES"
    $RUN "$WORK/fast/bin/mafft"  "$T" "$f" "fast|$r|$rel"  | tee -a "$RES"
  done
  if [ $CHECKDEF = 1 ]; then
    $RUN "$WORK/mafft-default/bin/mafft" "$T" "$f" "mafftdef|1|$rel" | tee -a "$RES"
    $RUN "$WORK/fast-default/bin/mafft"  "$T" "$f" "fastdef|1|$rel"  | tee -a "$RES"
  fi
  if [ $SELFCHECK = 1 ]; then
    out=$(MAFFT_SEXTET=check MAFFT_MSADIST=check MAFFT_UPG=check "$WORK/fast/bin/mafft" --thread "$T" "$f" 2>&1 >/dev/null); rc=$?
    if [ $rc -ne 0 ] || echo "$out" | grep -qiE "FAILED|mismatch"; then echo "SELFCHECK FAILED $rel (exit $rc)"; else echo "selfcheck passed $rel"; fi | tee -a "$RES"
  fi
done

echo; echo "== SUMMARY  (speedup = MAFFT wall / fast-mafft wall; identical = byte-identical output)"
python3 - "$RES" "$CHECKDEF" <<'EOF'
import re, sys, collections
rows = collections.OrderedDict(); fails = []
for line in open(sys.argv[1]):
    if line.startswith(("selfcheck", "SELFCHECK")):
        if line.startswith("SELFCHECK"): fails.append(line.strip())
        continue
    parts = line.split()
    who, rep, rel = parts[0].split("|", 2)
    kv = dict(re.findall(r"(\w+)=(\S+)", line))
    rows.setdefault(rel, {}).setdefault(rep, {})[who] = (float(kv["wall"]), float(kv["maxrss_gb"]), kv["md5"], kv["exit"])
print("| input | MAFFT mean (s) | fast-mafft mean (s) | speedup per rep | mean | min | identical | peak RSS MAFFT / fast (GB) |")
print("|---|---:|---:|---|---:|---:|---|---|")
allok = True; mins = []
for rel, reps in rows.items():
    sp, b, v, ok, rb, rv = [], [], [], True, [], []
    for rep, d in sorted(reps.items()):
        if "mafft" in d and "fast" in d:
            b.append(d["mafft"][0]); v.append(d["fast"][0]); sp.append(d["mafft"][0] / max(d["fast"][0], 1e-9))
            ok &= d["mafft"][2] == d["fast"][2] and d["mafft"][3] == "0" and d["fast"][3] == "0"
            rb.append(d["mafft"][1]); rv.append(d["fast"][1])
        if "mafftdef" in d and "fastdef" in d:
            ok &= d["mafftdef"][2] == d["fastdef"][2] == reps["1"]["mafft"][2]
    if not sp: continue
    allok &= ok; mins.append(min(sp))
    print("| %s | %.1f | %.1f | %s | %.1fx | %.1fx | %s | %.2f / %.2f |" % (rel, sum(b)/len(b), sum(v)/len(v),
          " / ".join("%.1fx" % s for s in sp), sum(sp)/len(sp), min(sp), "yes" if ok else "**NO**", max(rb), max(rv)))
print()
print("all outputs identical: %s;  lowest single-run speedup: %s" % ("YES" if allok else "NO", ("%.1fx" % min(mins)) if mins else "n/a"))
if sys.argv[2] == "1": print("default-flag builds checked: both outputs equal the aggressive builds' (included in 'identical')")
for f in fails: print(f)
EOF
echo "== done $(date)  (raw runs: $RES)"
