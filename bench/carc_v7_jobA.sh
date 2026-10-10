#!/bin/bash
#SBATCH --job-name=fastmafft-v7-A
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513|epyc-7542
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=10:00:00
#SBATCH --output=/scratch2/yihanhon/fast-mafft/carc-results/v7A-%j.log
# Account is passed on the sbatch command line (--account=...), never hard-coded here.
#
# Job A: (1) build aggressive variants on this node, (2) identity checks, (3) timing at --thread 32
# on large datasets up to ~28k sequences. Reference output = MAFFT 7.526 built with default flags.
set -u
P=/scratch2/yihanhon/fast-mafft
cd $P
source /etc/profile.d/modules.sh; module load gcc/13.3.0
OUT=$P/carc-results/v7A-$SLURM_JOB_ID; mkdir -p $OUT
echo "== node $(hostname) cpus=$SLURM_CPUS_PER_TASK $(lscpu | grep 'Model name' | sed 's/  */ /g') gcc=$(gcc -dumpfullversion)"
echo "== start $(date)"

AGG="-O3 -march=native -flto"
build() { # src prefix cflags
  rm -rf $OUT/bsrc-$(basename $2); cp -r $1 $OUT/bsrc-$(basename $2)
  make -C $OUT/bsrc-$(basename $2)/core clean >/dev/null 2>&1
  if [ -n "$3" ]; then
    make -j32 -C $OUT/bsrc-$(basename $2)/core CC=gcc CFLAGS="$3" PREFIX=$2 >$OUT/build-$(basename $2).log 2>&1 && \
    make -C $OUT/bsrc-$(basename $2)/core CC=gcc CFLAGS="$3" PREFIX=$2 install >>$OUT/build-$(basename $2).log 2>&1
  else
    make -j32 -C $OUT/bsrc-$(basename $2)/core CC=gcc PREFIX=$2 >$OUT/build-$(basename $2).log 2>&1 && \
    make -C $OUT/bsrc-$(basename $2)/core CC=gcc PREFIX=$2 install >>$OUT/build-$(basename $2).log 2>&1
  fi
  echo "built $(basename $2) rc=$? $(grep -m1 '^	prefix=/' $2/bin/mafft)"
}
build $P/baseline-src $OUT/base-aggr "$AGG"
build $P/work/v7/src  $OUT/v7-aggr   "$AGG"
BD=$P/baseline-default/bin/mafft          # reference (default flags, built on the login node)
VD=$P/work/v7/install-default/bin/mafft
BA=$OUT/base-aggr/bin/mafft
VA=$OUT/v7-aggr/bin/mafft
R="python3 $P/bench/carc_run.py"
D=$P/bench/data; L=$P/bench/data/large

m() { "$1" --thread "$2" "$3" 2>/dev/null | md5sum | cut -c1-12; }
echo "== IDENTITY (reference = baseline built with default flags)"
n=0; bad=0
for d in $D/sim1k.fa $D/rnasim579.fa $D/rnasim1k.fa $D/prot2k.fa; do
  for t in 1 8 32; do
    r=$(m $BD $t $d); x=$(m $BA $t $d); y=$(m $VD $t $d); z=$(m $VA $t $d); n=$((n+3))
    s="ok"; [ "$x" != "$r" ] && s="BASEAGGR-DIFFERS"; [ "$y" != "$r" ] && { s="$s V7DEF-DIFF"; bad=$((bad+1)); }; [ "$z" != "$r" ] && { s="$s V7AGG-DIFF"; bad=$((bad+1)); }
    echo "$(basename $d) t=$t ref=$r baseAggr=$x v7def=$y v7aggr=$z $s"
  done
done
for d in $D/edge.fa $D/edge_noshort.fa; do
  for t in 1 8; do
    r=$(m $BD $t $d); y=$(m $VD $t $d); z=$(m $VA $t $d); n=$((n+2))
    s="ok"; [ "$y" != "$r" ] && { s="$s V7DEF-DIFF"; bad=$((bad+1)); }; [ "$z" != "$r" ] && { s="$s V7AGG-DIFF"; bad=$((bad+1)); }
    echo "$(basename $d) t=$t ref=$r v7def=$y v7aggr=$z $s"
  done
done
echo "-- self-check modes (every 6-mer and msa distance compared bitwise with the original code)"
for d in $D/rnasim1k.fa $L/16S.T.fa $L/homfam_aldosered.fa; do
  MAFFT_SEXTET=check MAFFT_MSADIST=check $VA --thread 32 $d 2>&1 >/dev/null | grep -iE "check passed|check FAILED|mismatch|msadist fast" | head -3
done
echo "-- BAliBASE (v7 default flags vs reference, t=1 and t=8)"
bn=0; bb=0
for f in $(find $P/datasets/balibase/bb3_release -name '*.tfa' | sort); do
  for t in 1 8; do a=$(m $BD $t $f); b=$(m $VD $t $f); bn=$((bn+1)); [ "$a" != "$b" ] && { bb=$((bb+1)); echo "DIFF $f t=$t"; }; done
done
echo "BAliBASE runs=$bn differ=$bb"
echo "== IDENTITY SUMMARY: $n checks on standard sets, $bad differ; BAliBASE $bn runs, $bb differ"

echo "== TIMING --thread 32, 1 rep (baseline = aggressive build if its output matched default; v7 aggressive)"
for d in $D/rnasim5k.fa $D/rnasim10k.fa $L/16S.T.fa $L/homfam_p450.fa $L/rnasim20k.fa $L/rnasim25k.fa $L/16S.B.ALL.fa; do
  $R $BA 32 $d baseAggr
  $R $VA 32 $d v7aggr
  $R $VD 32 $d v7def
done
echo "== end $(date)"
