#!/bin/bash
#SBATCH --job-name=fastmafft-v7-B
#SBATCH --partition=main
#SBATCH --constraint=epyc-7513|epyc-7542
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=16:00:00
#SBATCH --output=/scratch2/yihanhon/fast-mafft/carc-results/v7B-%j.log
# Account is passed on the sbatch command line (--account=...).
# Job B: the very large sets (50k-94k sequences), --thread 32, 1 rep. Uses the aggressive builds made by job A
# (pass A's output dir as $1). Memory is the risk here: 64 GB limit; ru_maxrss is reported per run.
set -u
P=/scratch2/yihanhon/fast-mafft; A=$1
cd $P
source /etc/profile.d/modules.sh; module load gcc/13.3.0
echo "== node $(hostname) $(lscpu | grep 'Model name' | sed 's/  */ /g')  start $(date)"
R="python3 $P/bench/carc_run.py"; L=$P/bench/data/large
for d in $L/homfam_sdr.fa $L/rnasim50k.fa $L/homfam_rvp.fa; do
  $R $A/base-aggr/bin/mafft 32 $d baseAggr
  $R $A/v7-aggr/bin/mafft 32 $d v7aggr
done
echo "== end $(date)"
