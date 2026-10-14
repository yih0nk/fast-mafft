#!/bin/bash
#SBATCH --job-name=fastmafft-v9-final
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch2/yihanhon/fast-mafft/carc-results/v9final-%j.log
# Final measurement, one repetition per job (run 3 jobs in parallel, REP=1..3 as $1).
# Builds v9 (aggressive + default) on this node; times baseline (aggressive, whose output equals the default
# build's) and v9 (aggressive) interleaved on the same node at --thread 32. REP=1 also runs the identity suite.
# Account passed on the sbatch line.
set -u
REP=$1
P=/scratch2/yihanhon/fast-mafft
cd $P
source /etc/profile.d/modules.sh; module load gcc/13.3.0
OUT=$P/carc-results/v9final-$SLURM_JOB_ID; mkdir -p $OUT
echo "== rep $REP node $(hostname) $(lscpu | grep 'Model name' | sed 's/  */ /g') gcc=$(gcc -dumpfullversion) start $(date)"
build() { # prefix cflags
  rm -rf $OUT/bsrc-$(basename $1); cp -r $P/work/v9/src $OUT/bsrc-$(basename $1)
  make -C $OUT/bsrc-$(basename $1)/core clean >/dev/null 2>&1
  if [ -n "$2" ]; then F="CFLAGS=$2"; else F=""; fi
  make -j32 -C $OUT/bsrc-$(basename $1)/core CC=gcc ${2:+CFLAGS="$2"} PREFIX=$1 >$OUT/build-$(basename $1).log 2>&1 && \
  make -C $OUT/bsrc-$(basename $1)/core CC=gcc ${2:+CFLAGS="$2"} PREFIX=$1 install >>$OUT/build-$(basename $1).log 2>&1
  echo "built $(basename $1) rc=$? $(grep -m1 '^	prefix=/' $1/bin/mafft)"
}
build $OUT/v9-aggr "-O3 -march=native -flto"
build $OUT/v9-default ""
BD=$P/baseline-default/bin/mafft
BA=$P/carc-results/v7A-12999374/base-aggr/bin/mafft   # same flags, built on an EPYC 7513 node
VA=$OUT/v9-aggr/bin/mafft; VD=$OUT/v9-default/bin/mafft
R="python3 $P/bench/carc_run.py"; D=$P/bench/data; L=$P/bench/data/large
m() { "$1" --thread "$2" "$3" 2>/dev/null | md5sum | cut -c1-12; }

if [ "$REP" = 1 ]; then
  echo "== IDENTITY (reference = baseline built with default flags)"
  n=0; bad=0
  for d in $D/sim1k.fa $D/rnasim579.fa $D/rnasim1k.fa $D/prot2k.fa $D/rnasim2k.fa; do
    for t in 1 8 32; do
      r=$(m $BD $t $d); y=$(m $VD $t $d); z=$(m $VA $t $d); n=$((n+2))
      s=ok; [ "$y" != "$r" ] && { s="$s V9DEF-DIFF"; bad=$((bad+1)); }; [ "$z" != "$r" ] && { s="$s V9AGG-DIFF"; bad=$((bad+1)); }
      echo "$(basename $d) t=$t ref=$r v9def=$y v9aggr=$z $s"
    done
  done
  for d in $D/edge.fa $D/edge_noshort.fa $L/16S.3.fa; do
    for t in 1 8; do
      r=$(m $BD $t $d); y=$(m $VD $t $d); z=$(m $VA $t $d); n=$((n+2))
      s=ok; [ "$y" != "$r" ] && { s="$s V9DEF-DIFF"; bad=$((bad+1)); }; [ "$z" != "$r" ] && { s="$s V9AGG-DIFF"; bad=$((bad+1)); }
      echo "$(basename $d) t=$t ref=$r v9def=$y v9aggr=$z $s"
    done
  done
  echo "-- self-check modes"
  for d in $D/rnasim1k.fa $L/16S.T.fa $L/homfam_aldosered.fa; do
    MAFFT_SEXTET=check MAFFT_MSADIST=check $VA --thread 32 $d 2>&1 >/dev/null | grep -iE "check passed|FAILED|mismatch" | head -2
  done
  echo "-- BAliBASE (v9 default vs reference, t=1 and t=8)"
  bn=0; bb=0
  for f in $(find $P/datasets/balibase/bb3_release -name '*.tfa' | sort); do
    for t in 1 8; do a=$(m $BD $t $f); b=$(m $VD $t $f); bn=$((bn+1)); [ "$a" != "$b" ] && { bb=$((bb+1)); echo "DIFF $f t=$t"; }; done
  done
  echo "== IDENTITY SUMMARY: $n checks, $bad differ; BAliBASE $bn runs, $bb differ"
fi

echo "== TIMING rep $REP, --thread 32 (md5 of every timed run must match between baseline and v9)"
for d in $D/rnasim5k.fa $D/rnasim10k.fa $L/16S.T.fa $L/16S.3.fa $L/homfam_p450.fa $L/rnasim20k.fa $L/16S.B.ALL.fa \
         $L/homfam_sdr.fa $L/homfam_Acetyltransf.fa $L/rnasim25k.fa $L/homfam_rvp.fa $L/rnasim50k.fa; do
  $R $BA 32 $d baseAggr-r$REP
  $R $VA 32 $d v9aggr-r$REP
done
echo "== end $(date)"
