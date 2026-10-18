# Evaluating on CARC

**To evaluate fast-mafft yourself, use `run_eval.sh`.** It is self-contained: it builds MAFFT v7.526
and fast-mafft from this repository on the job's node, runs both on any directory of FASTA files, and
checks byte-identity and speedup. See the top-level README, "How to evaluate", or run
`bash carc/run_eval.sh --help`.

The rest of this file documents how the numbers in `../RESULTS.md` were produced, and the logs in
`logs/`.

## The v10 evaluation (RESULTS.md §3)

`carc_v10_final.sh REP PART` is self-contained like `run_eval.sh`: submitted from the repository root,
it builds MAFFT v7.526 (from the first commit) and fast-mafft on its node and reads the inputs from
`bench/data` (`bash bench/get_data.sh --long --balibase`). Six jobs ran in parallel on 2026-10-18 at
commit `b35964d`:

```
for rep in 1 2 3; do for part in std long; do
  sbatch --account=<account> carc/carc_v10_final.sh $rep $part; done; done
python3 carc/summarize.py carc/logs/v10final-*.log        # the tables of RESULTS.md §3
```

| job | REP | PART | node | what |
|---|---:|---|---|---|
| 13224806 | 1 | std | a03-19 | identity suite (45 checks, self-checks, BAliBASE 772 runs) + timing |
| 13224807 | 1 | long | a03-19 | identity suite (23 checks) + timing |
| 13224808 | 2 | std | a04-05 | timing |
| 13224809 | 2 | long | a04-05 | timing |
| 13224810 | 3 | std | b10-02 | timing |
| 13224811 | 3 | long | b10-03 | timing |

PART=std times sim1k, rnasim1k, prot2k, sim2k, rnasim2k and the 12 datasets of the v9 evaluation;
PART=long times sarscov2_500, sarscov2_1000, edge, long100x30k, long1000x30k, long300x50k,
long500x30k and sarscov2_2000. Each timed line has the format below (labels `baseAggr-rK` and
`v10aggr-rK`; `baseAgain` / `v10Again` for the reruns after an md5 difference).

## The v9 evaluation (RESULTS.md §4)

The `carc_v9_final.sh` jobs, whose paths point at the author's CARC scratch directory.

## Files

| file | what it is |
|---|---|
| `carc_v10_final.sh` | Slurm job script of the v10 evaluation (above). |
| `carc_v10_mix.sh`, `diag_mix.py` | The diagnostic job of RESULTS.md §3.4: repeated runs of both programs, with builds that log every step's aligner, on inputs that switch to the memsave mode in the middle of a pass. |
| `carc_v9_final.sh` | Slurm job script of the v9 evaluation, a verbatim copy of `bench/carc_v9_final.sh`. One repetition per job. |
| `carc_run.py` | Timing harness, a verbatim copy of `bench/carc_run.py`. It runs one `mafft --thread N` and prints one result line. |
| `summarize.py` | Turns the job logs (v9 or v10) into the tables of `../RESULTS.md`. |
| `run_eval.sh` | The portable evaluation job (top-level README). |
| `logs/` | `v10final-1322480{6..9}.log`, `v10final-1322481{0,1}.log`, `v10mix-13226265.log`, `v9final-1307213{6,7,8}.log`. |

### What one v9 job does (REP = 1, 2 or 3)

1. Builds fast-mafft v9 twice on its node with gcc 13.3: once with `-O3 -march=native -flto`
   (`v9-aggr`) and once with the default flags (`v9-default`).
2. **Only when REP = 1:** runs the identity suite against `baseline-default`, v7.526 built with the
   default flags:
   - 42 md5 checks;
   - the `MAFFT_SEXTET=check MAFFT_MSADIST=check` self-checks;
   - BAliBASE at 1 and 8 threads.
3. For each of the 12 datasets in order, times the aggressive baseline and then `v9-aggr` at
   `--thread 32`, interleaved on the same node.

One quirk: the script's header comment still says "Builds v8"; that is a leftover. The script builds
`$P/work/v9/src`.

### Result line format

Each timed run prints one line:

```
baseAggr-r1 rnasim10k t=32 wall=92.9 maxrss_gb=3.71 md5=xxxxxxxxxxxx exit=0 stages=6mer:..,prog1:..,distmsa:..,prog2:..
v9aggr-r1   rnasim10k t=32 wall=...
```

The stage windows run from one MAFFT stderr marker to the next:

- `6mer` covers the 6-mer distances and guide tree 1;
- `prog1` covers progressive pass 1;
- `distmsa` covers the msa distances and guide tree 2;
- `prog2` runs from pass 2 to the `Strategy:` line, so it includes writing the output.

`maxrss_gb` is Linux `ru_maxrss` (KiB) of the largest child process, converted to GB. On macOS the
same field is in bytes, so the number is meaningless there.

## Paths the script expects (`P=/scratch2/yihanhon/fast-mafft`)

| path | contents |
|---|---|
| `$P/work/v9/src/` | fast-mafft v9 source: this package's `fast-mafft-src/`. Sync it with `rsync -a --delete fast-mafft-src/ carc:$P/work/v9/src/`. Quick check: `md5sum core/disttbfast.c core/mltaln9.c` must give `5a2011f18d10…` and `6214cdc588d6…`. |
| `$P/baseline-default/` | MAFFT v7.526 installed with gcc 13.3 default flags (the identity reference) |
| `$P/carc-results/v7A-12999374/base-aggr/` | MAFFT v7.526 installed with gcc 13.3 `-O3 -march=native -flto`, built on an EPYC 7513 node (the timed baseline) |
| `$P/bench/carc_run.py` | the harness |
| `$P/bench/data/*.fa` | sim1k, rnasim579, rnasim1k, prot2k, rnasim2k, rnasim5k, rnasim10k, edge, edge_noshort |
| `$P/bench/data/large/*.fa` | 16S.T, 16S.3, 16S.B.ALL, homfam_p450, homfam_aldosered, homfam_sdr, homfam_Acetyltransf, homfam_rvp, rnasim20k, rnasim25k, rnasim50k |
| `$P/datasets/balibase/bb3_release/` | BAliBASE 3.0 core (`*.tfa`) |
| `$P/carc-results/` | the job writes `v9final-<jobid>.log` and a `v9final-<jobid>/` directory of builds here |

To rebuild the two baselines from a clean v7.526 tree, use the `build()` function in `../BUILD.md`
§2 with gcc 13.3, then point `BD` and `BA` in the script at them. For a different user or scratch
directory, edit `P` and the `#SBATCH --output` line.

## Submit

The script sets the partition, constraint (`epyc-7513`), 32 CPUs, 64 GB and an 8-hour limit. Pass
the account on the command line; the development runs used the course allocation. Submit the three
repetitions together so that they land on three nodes:

```
cd $P
for r in 1 2 3; do sbatch --account=<allocation> bench/carc_v9_final.sh $r; done
squeue -u $USER
```

Running time: the baseline dominates.

- Above 20,000 sequences it runs its progressive passes on one core. homfam_rvp alone took about
  950–1220 s in earlier jobs.
- REP = 1 also runs the identity suite and 2 × 772 BAliBASE runs before timing.

## Summarize

```
python3 carc/summarize.py $P/carc-results/v9final-<job1>.log $P/carc-results/v9final-<job2>.log $P/carc-results/v9final-<job3>.log
```

This prints the summary table (§3.1 of RESULTS.md) and the raw-run table (§3.2) as Markdown.

- **Speedups.** Per-run speedup = baseline wall / v9 wall within the same job. "Mean speedup" is the
  mean of the three; "min speedup" is the smallest.
- **md5 match.** A row shows "NO / incomplete" if any repetition is missing, exited non-zero, or
  printed a different md5 from its baseline.
- **Identity suite.** Copy its summary lines from the REP = 1 log into §3.3:

  ```
  grep -E '== IDENTITY SUMMARY|check passed|FAILED|mismatch|DIFF' v9final-<job1>.log
  ```
