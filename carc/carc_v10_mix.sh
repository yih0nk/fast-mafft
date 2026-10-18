#!/bin/bash
#SBATCH --job-name=fm-v10-mix
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=v10mix-%j.log
#
# Inputs whose alignment crosses 30,000 columns in the middle of a pass (VERIFICATION.md §3): runs MAFFT
# v7.526 and fast-mafft repeatedly on the same node, with builds that log every step's aligner and every
# Falign call split between the two aligners (diag_mix.py; stdout unchanged), and keeps each run's stderr.
#   sbatch --account=<account> carc/carc_v10_mix.sh          (from the repository root; bench/data/long)
set -u
REPO=${SLURM_SUBMIT_DIR:-$PWD}; cd $REPO
OUT=$REPO/carc-results/v10mix-${SLURM_JOB_ID:-local}; mkdir -p $OUT/err
source /etc/profile.d/modules.sh; module load gcc/13.3.0
echo "== node $(hostname) start $(date)  repo $(git rev-parse --short HEAD)"
mkdir -p $OUT/archive && git archive $(git rev-list --max-parents=0 HEAD) fast-mafft-src | tar -x -C $OUT/archive
for b in base v10; do
  rm -rf $OUT/src-$b
  if [ $b = base ]; then mv $OUT/archive/fast-mafft-src $OUT/src-$b; else cp -r $REPO/fast-mafft-src $OUT/src-$b; fi
  make -C $OUT/src-$b/core clean >/dev/null 2>&1
  python3 $REPO/carc/diag_mix.py $OUT/src-$b
  make -j32 -C $OUT/src-$b/core CC=gcc CFLAGS="-O3 -march=native -flto" PREFIX=$OUT/$b >$OUT/build-$b.log 2>&1 && \
  make -C $OUT/src-$b/core CC=gcc CFLAGS="-O3 -march=native -flto" PREFIX=$OUT/$b install >>$OUT/build-$b.log 2>&1
  echo "built diag $b rc=$? $(grep -m1 '^	prefix=/' $OUT/$b/bin/mafft)"
done
G=$REPO/bench/data/long
run() { # BUILD FILE THREADS TAG
  local t0=$(date +%s.%N) m
  m=$($OUT/$1/bin/mafft --thread $3 $2 2>$OUT/err/$4.err | md5sum | cut -c1-12)
  local t1=$(date +%s.%N)
  gzip -f $OUT/err/$4.err
  echo "$4 $1 $(basename $2 .fa) t=$3 md5=$m wall=$(python3 -c "print('%.1f' % ($t1-$t0))")"
}
for i in 1 2 3 4 5 6 7 8; do run base $G/viral100.fa 8 v100t8-base-$i; run v10 $G/viral100.fa 8 v100t8-v10-$i; done
for i in 1 2 3 4; do run base $G/viral100.fa 32 v100t32-base-$i; run v10 $G/viral100.fa 32 v100t32-v10-$i; done
for i in $(seq 1 12); do run v10 $G/sarscov2_2000.fa 32 sc2k-v10-$i; done
run base $G/sarscov2_2000.fa 32 sc2k-base-1
echo "== end $(date)"
echo "== the switch steps and split calls of each run: zcat \$OUT/err/<tag>.err.gz | tr '\\r' '\\n' | grep -E 'SWITCH|MIXED|TOP'"
