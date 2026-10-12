#!/bin/bash
#SBATCH --job-name=fastmafft-v7-C
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch2/yihanhon/fast-mafft/carc-results/v7C-%j.log
# Job C: repeat timing (2 more reps -> 3 total with jobs A/B), --thread 32, baseline (aggressive) vs v7 (aggressive),
# interleaved per dataset. Uses the builds from job A ($1 = A's output dir). Account passed on the sbatch line.
set -u
P=/scratch2/yihanhon/fast-mafft; A=$1
cd $P
source /etc/profile.d/modules.sh; module load gcc/13.3.0
echo "== node $(hostname) $(lscpu | grep 'Model name' | sed 's/  */ /g')  start $(date)"
R="python3 $P/bench/carc_run.py"; D=$P/bench/data; L=$P/bench/data/large
for rep in 2 3; do
  for d in $D/rnasim5k.fa $D/rnasim10k.fa $L/16S.T.fa $L/homfam_p450.fa $L/rnasim20k.fa $L/16S.B.ALL.fa \
           $L/homfam_sdr.fa $L/rnasim25k.fa $L/homfam_rvp.fa $L/rnasim50k.fa; do
    $R $A/base-aggr/bin/mafft 32 $d baseAggr-r$rep
    $R $A/v7-aggr/bin/mafft 32 $d v7aggr-r$rep
  done
done
echo "== end $(date)"
