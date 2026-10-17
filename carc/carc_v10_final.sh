#!/bin/bash
#SBATCH --job-name=fastmafft-v10-final
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=v10final-%j.log
#
# Final measurement of fast-mafft v10: one repetition of one part per job. Submit from the repository
# root (a git checkout: MAFFT v7.526 is built from its first commit), six jobs in parallel:
#   for rep in 1 2 3; do for part in std long; do
#     sbatch --account=<account> carc/carc_v10_final.sh $rep $part; done; done
#
#   PART=std   the 12 datasets of the v9 evaluation, plus 5 small ones (1,000-2,000 sequences)
#   PART=long  real SARS-CoV-2 collections and long-genome sets (alignments near or over 30,000 columns)
#
# Builds MAFFT v7.526 and fast-mafft on this node with gcc 13.3 and -O3 -march=native -flto, and times them
# interleaved at --thread 32 (carc_run.py: wall time, md5 of stdout, peak RSS). On an md5 difference both
# programs run again, to tell whether MAFFT itself gives the same output twice on that input. REP=1 also
# runs the identity suite of its part against MAFFT v7.526 built with the Makefile's default flags.
# Data: bench/data (bash bench/get_data.sh --balibase) and bench/data/long (bench/DATA.md).
set -u
REP=$1; PART=$2
REPO=${SLURM_SUBMIT_DIR:-$PWD}
cd $REPO
D=$REPO/bench/data; L=$D/large; G=$D/long; BB=$D/balibase/bb3_release
OUT=$REPO/carc-results/v10final-${SLURM_JOB_ID:-local}; mkdir -p $OUT
source /etc/profile.d/modules.sh; module load gcc/13.3.0
echo "== rep $REP part $PART node $(hostname) $(lscpu | grep 'Model name' | sed 's/  */ /g') gcc=$(gcc -dumpfullversion) start $(date)"
echo "== repo $(git rev-parse --short HEAD)  Cpus_allowed_list: $(grep Cpus_allowed_list /proc/self/status | cut -f2)"

# ---- sources: MAFFT v7.526 = the repository's first commit (checked); fast-mafft = fast-mafft-src/
rm -rf $OUT/src-base; mkdir -p $OUT/archive
git archive $(git rev-list --max-parents=0 HEAD) fast-mafft-src | tar -x -C $OUT/archive && mv $OUT/archive/fast-mafft-src $OUT/src-base && rm -rf $OUT/archive
[ "$(md5sum < $OUT/src-base/core/disttbfast.c | cut -c1-32)" = 564c5519e44b43d888f14b8badb23032 ] || { echo "baseline source is not MAFFT v7.526"; exit 3; }
build() { # NAME SRC [CFLAGS]: a fresh copy, make clean, the same PREFIX on both make calls
  rm -rf $OUT/bsrc-$1; cp -r $2 $OUT/bsrc-$1
  make -C $OUT/bsrc-$1/core clean >/dev/null 2>&1
  make -j32 -C $OUT/bsrc-$1/core CC=gcc ${3:+CFLAGS="$3"} PREFIX=$OUT/$1 >$OUT/build-$1.log 2>&1 && \
  make -C $OUT/bsrc-$1/core CC=gcc ${3:+CFLAGS="$3"} PREFIX=$OUT/$1 install >>$OUT/build-$1.log 2>&1
  echo "built $1 rc=$? $(grep -m1 '^	prefix=/' $OUT/$1/bin/mafft)"
}
AGGR="-O3 -march=native -flto"
build base-aggr $OUT/src-base "$AGGR"
build v10-aggr $REPO/fast-mafft-src "$AGGR"
BA=$OUT/base-aggr/bin/mafft; VA=$OUT/v10-aggr/bin/mafft
if [ "$REP" = 1 ]; then
  build base-default $OUT/src-base
  build v10-default $REPO/fast-mafft-src
fi
BD=$OUT/base-default/bin/mafft; VD=$OUT/v10-default/bin/mafft
R="python3 $REPO/carc/carc_run.py"
m() { "$1" --thread "$2" "$3" 2>/dev/null | md5sum | cut -c1-12; }

if [ "$REP" = 1 ]; then
  echo "== IDENTITY, part $PART (reference = MAFFT v7.526 built with the default flags)"
  n=0; bad=0
  three() { # FILE THREADS: reference, fast-mafft default build, fast-mafft aggressive build
    local r y z s; r=$(m $BD $2 $1); y=$(m $VD $2 $1); z=$(m $VA $2 $1); n=$((n+2))
    s=ok; [ "$y" != "$r" ] && { s="$s V10DEF-DIFF"; bad=$((bad+1)); }; [ "$z" != "$r" ] && { s="$s V10AGG-DIFF"; bad=$((bad+1)); }
    echo "$(basename $1) t=$2 ref=$r v10def=$y v10aggr=$z $s"
  }
  checked() { # FILE THREADS ENV...: the aggressive build with self-checks on (exit on any MISMATCH), vs the reference
    local d=$1 t=$2 r c s; shift 2
    r=$(m $BD $t $d); c=$(env "$@" $VA --thread $t $d 2>$OUT/chk-stderr.txt | md5sum | cut -c1-12); n=$((n+1))
    s=ok; [ "$c" != "$r" ] && { s=V10AGG-DIFF; bad=$((bad+1)); }
    echo "$(basename $d) t=$t $* ref=$r v10aggr=$c $s | $(tr '\r' '\n' < $OUT/chk-stderr.txt | grep -aE 'CHECK|MISMATCH|check passed|FAILED' | tr '\n' ' ')"
  }
  if [ "$PART" = std ]; then
    for d in $D/sim1k.fa $D/rnasim579.fa $D/rnasim1k.fa $D/prot2k.fa $D/rnasim2k.fa; do for t in 1 8 32; do three $d $t; done; done
    for d in $D/edge.fa $D/edge_noshort.fa $L/16S.3.fa; do for t in 1 8; do three $d $t; done; done
    echo "-- self-checks (6-mer distances, MSA distances, guide tree; blockAlign2, FFT lags, anchors)"
    for d in $D/rnasim1k.fa $L/16S.T.fa $L/homfam_aldosered.fa; do
      checked $d 32 MAFFT_SEXTET=check MAFFT_MSADIST=check MAFFT_UPG=check MAFFT_BA2CHECK=1 MAFFT_ARLCHECK=1 MAFFT_FFTCHECK=1
    done
    echo "-- BAliBASE 3.0 (fast-mafft default build vs reference, t=1 and t=8)"
    bn=0; bb=0
    for f in $(find $BB -name '*.tfa' | sort); do
      for t in 1 8; do a=$(m $BD $t $f); b=$(m $VD $t $f); bn=$((bn+1)); [ "$a" != "$b" ] && { bb=$((bb+1)); echo "DIFF $f t=$t"; }; done
    done
    echo "== IDENTITY SUMMARY: $n checks, $bad differ; BAliBASE $bn runs, $bb differ"
  else
    echo "-- strict one-thread identity, real SARS-CoV-2 collections"
    three $G/sarscov2_500.fa 1; three $G/sarscov2_1000.fa 1
    echo "-- long-genome sets (switch to MAFFT's memsave mode in the first steps of pass 1, or mid-pass for viral100)"
    for d in $G/long10x30k.fa $G/long30x30k.fa $G/viral100.fa; do for t in 1 8; do three $d $t; done; done
    echo "-- self-checks: blockAlign2, the long FFTs, anchors and lags; Myers-Miller calls and recursion nodes"
    checked $G/sarscov2_500.fa 8 MAFFT_BA2CHECK=1 MAFFT_FFTCHECK=1 MAFFT_ARLCHECK=1
    checked $G/sarscov2_1000.fa 8 MAFFT_BA2CHECK=1 MAFFT_FFTCHECK=1 MAFFT_ARLCHECK=1
    checked $G/sarscov2_500.fa 32 MAFFT_BA2CHECK=1 MAFFT_FFTCHECK=1 MAFFT_ARLCHECK=1
    checked $G/long30x30k.fa 32 MAFFT_MSKCHECK=1 MAFFT_MSCHECK=1 MAFFT_ARLCHECK=1 MAFFT_FFTCHECK=1 MAFFT_BA2CHECK=1
    echo "-- every SIMD level (runtime dispatch), long10x30k at one thread"
    for s in scalar sse2 avx2; do checked $G/long10x30k.fa 1 MAFFT_SIMD=$s MAFFT_MSCHECK=1; done
    echo "== IDENTITY SUMMARY: $n checks, $bad differ"
  fi
fi

echo "== TIMING rep $REP part $PART, --thread 32, baseline then fast-mafft on each input (same node)"
if [ "$PART" = std ]; then
  SETS="$D/sim1k.fa $D/rnasim1k.fa $D/prot2k.fa $D/sim2k.fa $D/rnasim2k.fa
        $D/rnasim5k.fa $D/rnasim10k.fa $L/16S.T.fa $L/16S.3.fa $L/homfam_p450.fa $L/rnasim20k.fa $L/16S.B.ALL.fa
        $L/homfam_sdr.fa $L/homfam_Acetyltransf.fa $L/rnasim25k.fa $L/homfam_rvp.fa $L/rnasim50k.fa"
else
  SETS="$G/sarscov2_500.fa $G/sarscov2_1000.fa $D/edge.fa $G/long100x30k.fa $G/long1000x30k.fa
        $G/long300x50k.fa $G/long500x30k.fa $G/sarscov2_2000.fa"
fi
for d in $SETS; do
  [ -s $d ] || { echo "MISSING $d"; continue; }
  lb=$($R $BA 32 $d baseAggr-r$REP); echo "$lb"
  lc=$($R $VA 32 $d v10aggr-r$REP); echo "$lc"
  mb=$(echo "$lb" | sed 's/.*md5=\([0-9a-f]*\).*/\1/'); mc=$(echo "$lc" | sed 's/.*md5=\([0-9a-f]*\).*/\1/')
  if [ "$mb" != "$mc" ]; then
    $R $BA 32 $d baseAgain-r$REP; $R $VA 32 $d v10Again-r$REP
  fi
done
echo "== end $(date)"
