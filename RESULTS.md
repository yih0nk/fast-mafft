# CARC evaluation results

**Result: on all 12 datasets, at `--thread 32` on a 32-core / 64 GB CARC node, fast-mafft v9 produced
output byte-identical to MAFFT v7.526 and ran 9.7x to 36.2x faster on average. The slowest single
run of any dataset was 8.3x** (16S.T, run 1). The baseline is MAFFT v7.526 built with
`-O3 -march=native -flto`; its output equals the default build's.

Source: three final CARC jobs of [`carc/carc_v9_final.sh`](carc/carc_v9_final.sh), REP = 1, 2 and 3.
The logs are in [`carc/logs/`](carc/logs/), and [`carc/summarize.py`](carc/summarize.py) regenerates
the tables below from them.

## 1. Methodology

| item | setting |
|---|---|
| Hardware | USC CARC `main` partition, `--constraint=epyc-7513` (AMD EPYC 7513, Zen 3, AVX2, no AVX-512). One node per job: `--nodes=1 --ntasks=1 --cpus-per-task=32 --mem=64G`, `--time=08:00:00`. The 32 cores of a job can be spread over both sockets of a shared node. |
| Compiler | gcc 13.3.0 (`module load gcc/13.3.0`). |
| Baseline (timed) | MAFFT v7.526 built with gcc 13.3 and `-O3 -march=native -flto` on an EPYC 7513 node (`carc-results/v7A-12999374/base-aggr`). Its output equals the default build's. |
| Identity reference | MAFFT v7.526 built with gcc 13.3 and the Makefile's default flags (`baseline-default`). The default build is the reference whenever the output could depend on how it was compiled. |
| fast-mafft v9 | Built inside each job, on that node, from `fast-mafft-src`. It is built with `-O3 -march=native -flto` (timed) and with default flags (identity suite). |
| Command | `bin/mafft --thread 32 in.fa`, run through `carc/carc_run.py`. That script records wall-clock time from launch to exit, the md5 of stdout, peak RSS (`ru_maxrss` of the largest child process, i.e. `disttbfast`) and stage windows from MAFFT's stderr markers. |
| Interleaving | For each of the 12 datasets in the order below, the baseline runs and then v9 runs, on the same node in the same job. |
| Repetitions | 3 jobs (REP = 1, 2, 3) submitted together, so 3 nodes. REP = 1 first runs the identity suite (§3.3). |
| Per-run speedup | Baseline wall time / v9 wall time within the same job (same node). |
| Mean speedup | Arithmetic mean of the three per-run speedups. |
| Min speedup | The smallest of the three per-run speedups. |
| md5 match | v9's md5 equals the baseline's in all three repetitions. |
| Peak RSS | The largest of the three runs, in GB. |

## 2. Datasets

Lengths are ungapped, computed from the input files. The input md5 is the first 12 hex digits of the
FASTA file.

| # | dataset | sequences | length min / avg / max | type | source | input md5 | notes |
|---:|---|---:|---|---|---|---|---|
| 1 | rnasim5k | 5,000 | 1,526 / 1,555 / 1,587 | RNA (ACGU) | 5,000 sequences of rnasim10k (random seed 5, original order) | 3d39ddd41994 | |
| 2 | rnasim10k | 10,000 | 1,516 / 1,555 / 1,590 | RNA (ACGU) | RNASim 10K (`RNASim/10000_R_1`) true alignment, gaps removed | 217a0f12193d | |
| 3 | 16S.T | 7,350 | 506 / 1,492 / 4,066 | rRNA, IUPAC codes | 16S.T from the Recursive MAGUS archive | 3b32c0689f16 | sequences up to 4,066 nt; alignment stays under 30,000 columns |
| 4 | 16S.3 | 6,323 | 732 / 1,557 / 4,066 | rRNA, IUPAC codes | 16S.3, same archive | eb6842523936 | |
| 5 | homfam_p450 | 21,013 | 21 / 332 / 557 | protein | HomFam family p450, same archive | f6414db14eac | > 20,000: baseline progressive passes are serial |
| 6 | rnasim20k | 20,000 | 1,522 / 1,555 / 1,592 | RNA (ACGU) | RNASim 20K (`RNASim/20000_R_1`), same archive | 247220301d41 | exactly 20,000: baseline stays threaded |
| 7 | 16S.B.ALL | 27,643 | 266 / 1,372 / 1,851 | rRNA, IUPAC codes | 16S.B.ALL, same archive | da3d28d3e1e8 | > 20,000 |
| 8 | homfam_sdr | 50,157 | 21 / 163 / 297 | protein | HomFam family sdr, same archive | d7f891074991 | > 20,000; about 10 GB distance matrix |
| 9 | homfam_Acetyltransf | 46,285 | 28 / 83 / 181 | protein | HomFam family Acetyltransf, same archive | 45db38a70017 | > 20,000 |
| 10 | rnasim25k | 25,000 | 1,518 / 1,555 / 1,592 | RNA (ACGU) | 25,000 sequences of rnasim50k (`random.Random(25).sample`, original order) | 45ffeb52a2ae | > 20,000 |
| 11 | homfam_rvp | 93,681 | 22 / 94 / 133 | protein | HomFam family rvp, same archive | 9ec8a71d2ea7 | > 20,000; about 35 GB distance matrix (baseline peak RSS 32.97 GB) |
| 12 | rnasim50k | 50,000 | 1,517 / 1,554 / 1,593 | RNA (ACGU) | RNASim 50K (`RNASim/50000_R_1`), same archive | 0620e62aaa34 | > 20,000; about 10 GB distance matrix |

**Sources.**

- **RNASim, 16S and HomFam.** Smirnov, "Recursive MAGUS" (PLoS Comput Biol 2021), datasets on the
  Illinois Data Bank, doi:10.13012/B2IDB-1048258_V1, license CC0 (`data.zip`, 3.69 GB). Unaligned
  FASTA was extracted from it.
  - RNASim is the Guo/Wang/Kim one-million-sequence simulation of 16S-like RNA (arXiv 0912.2326).
  - HomFam is the Clustal Omega benchmark of large Pfam families with Homstrad references.
- **BAliBASE 3.0** (Thompson, Koehl, Ripp and Poch, Proteins 61:127–136, 2005), core release R1–5:
  386 `.tfa` inputs (RV11–RV50). It is used for identity only (REP = 1), not timing, because each case
  runs in well under a second.
- **Smaller identity sets** (§3.3): synthetic sim1k and prot2k (`bench/simulate*.py`), rnasim579
  (TWILIGHT's GitHub sample, MIT), the rnasim1k and rnasim2k subsets of RNASim 10K, and edge.fa and
  edge_noshort.fa (alignments wider than 30,000 columns).

## 3. Final results

### 3.1 Summary, `--thread 32`

| dataset | baseline mean (s) | v9 mean (s) | speedup r1 | speedup r2 | speedup r3 | mean speedup | min speedup | md5 match | peak RSS baseline / v9 (GB) |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| rnasim5k | 45.8 | 4.7 | 8.72x | 10.47x | 10.02x | 9.74x | 8.72x | yes (3/3) | 3.03 / 3.23 |
| rnasim10k | 96.4 | 8.8 | 9.60x | 11.61x | 12.03x | 11.08x | 9.60x | yes (3/3) | 5.77 / 3.66 |
| 16S.T | 59.1 | 5.9 | 8.29x | 10.78x | 11.64x | 10.24x | 8.29x | yes (3/3) | 3.75 / 4.07 |
| 16S.3 | 47.6 | 4.7 | 8.84x | 10.35x | 11.65x | 10.28x | 8.84x | yes (3/3) | 3.10 / 3.02 |
| homfam_p450 | 144.2 | 7.7 | 14.95x | 20.56x | 22.58x | 19.36x | 14.95x | yes (3/3) | 2.46 / 2.59 |
| rnasim20k | 234.6 | 18.7 | 11.66x | 12.87x | 13.30x | 12.61x | 11.66x | yes (3/3) | 12.86 / 5.46 |
| 16S.B.ALL | 509.3 | 21.6 | 19.40x | 24.85x | 28.35x | 24.20x | 19.40x | yes (3/3) | 3.84 / 4.49 |
| homfam_sdr | 272.9 | 15.4 | 13.79x | 19.60x | 21.56x | 18.32x | 13.79x | yes (3/3) | 9.97 / 10.27 |
| homfam_Acetyltransf | 186.3 | 12.7 | 9.94x | 19.30x | 19.37x | 16.20x | 9.94x | yes (3/3) | 8.47 / 8.69 |
| rnasim25k | 1199.8 | 34.9 | 25.73x | 38.28x | 44.62x | 36.21x | 25.73x | yes (3/3) | 4.00 / 4.64 |
| homfam_rvp | 1002.4 | 32.9 | 24.96x | 32.39x | 38.92x | 32.09x | 24.96x | yes (3/3) | 33.00 / 33.20 |
| rnasim50k | 2523.6 | 88.2 | 23.46x | 28.26x | 37.15x | 29.62x | 23.46x | yes (3/3) | 12.59 / 12.77 |

Nodes, all AMD EPYC 7513 32-Core with gcc 13.3.0: r1 = a04-03 (job 13072136), r2 = a03-17 (job 13072137),
r3 = b09-04 (job 13072138). On each node the baseline and v9 ran back to back, one dataset at a time.
Run-to-run spread comes from sharing a 64-core node with other users' jobs. v9 keeps all 32 cores busy,
so it is more sensitive to a noisy neighbour than the baseline, whose long phases are single-threaded
for large inputs.

### 3.2 Raw runs

| dataset | rep | baseline wall (s) | v9 wall (s) | speedup | baseline md5 | v9 md5 | baseline RSS (GB) | v9 RSS (GB) |
|---|---:|---:|---:|---:|---|---|---:|---:|
| rnasim5k | 1 | 46.2 | 5.3 | 8.72x | 455da94efd24 | 455da94efd24 | 2.90 | 3.23 |
| rnasim5k | 2 | 47.1 | 4.5 | 10.47x | 455da94efd24 | 455da94efd24 | 2.98 | 3.12 |
| rnasim5k | 3 | 44.1 | 4.4 | 10.02x | 455da94efd24 | 455da94efd24 | 3.03 | 3.12 |
| rnasim10k | 1 | 97.9 | 10.2 | 9.60x | 5e7bb28c8e8e | 5e7bb28c8e8e | 5.77 | 3.63 |
| rnasim10k | 2 | 96.4 | 8.3 | 11.61x | 5e7bb28c8e8e | 5e7bb28c8e8e | 5.74 | 3.66 |
| rnasim10k | 3 | 95.0 | 7.9 | 12.03x | 5e7bb28c8e8e | 5e7bb28c8e8e | 5.68 | 3.44 |
| 16S.T | 1 | 59.7 | 7.2 | 8.29x | f487c1a56702 | f487c1a56702 | 3.75 | 4.00 |
| 16S.T | 2 | 59.3 | 5.5 | 10.78x | f487c1a56702 | f487c1a56702 | 3.53 | 3.99 |
| 16S.T | 3 | 58.2 | 5.0 | 11.64x | f487c1a56702 | f487c1a56702 | 3.73 | 4.07 |
| 16S.3 | 1 | 48.6 | 5.5 | 8.84x | 6b646c6b3c1f | 6b646c6b3c1f | 3.10 | 2.88 |
| 16S.3 | 2 | 47.6 | 4.6 | 10.35x | 6b646c6b3c1f | 6b646c6b3c1f | 3.09 | 3.02 |
| 16S.3 | 3 | 46.6 | 4.0 | 11.65x | 6b646c6b3c1f | 6b646c6b3c1f | 3.10 | 3.00 |
| homfam_p450 | 1 | 146.5 | 9.8 | 14.95x | a8da25a0a301 | a8da25a0a301 | 2.46 | 2.54 |
| homfam_p450 | 2 | 146.0 | 7.1 | 20.56x | a8da25a0a301 | a8da25a0a301 | 2.46 | 2.57 |
| homfam_p450 | 3 | 140.0 | 6.2 | 22.58x | a8da25a0a301 | a8da25a0a301 | 2.46 | 2.59 |
| rnasim20k | 1 | 234.3 | 20.1 | 11.66x | 89f21e3ddf15 | 89f21e3ddf15 | 12.07 | 5.03 |
| rnasim20k | 2 | 236.8 | 18.4 | 12.87x | 89f21e3ddf15 | 89f21e3ddf15 | 12.17 | 5.15 |
| rnasim20k | 3 | 232.8 | 17.5 | 13.30x | 89f21e3ddf15 | 89f21e3ddf15 | 12.86 | 5.46 |
| 16S.B.ALL | 1 | 506.3 | 26.1 | 19.40x | 8a066cfcbfbe | 8a066cfcbfbe | 3.83 | 4.44 |
| 16S.B.ALL | 2 | 516.8 | 20.8 | 24.85x | 8a066cfcbfbe | 8a066cfcbfbe | 3.83 | 4.49 |
| 16S.B.ALL | 3 | 504.7 | 17.8 | 28.35x | 8a066cfcbfbe | 8a066cfcbfbe | 3.84 | 4.41 |
| homfam_sdr | 1 | 269.0 | 19.5 | 13.79x | 4d0ab395b00f | 4d0ab395b00f | 9.96 | 10.18 |
| homfam_sdr | 2 | 282.3 | 14.4 | 19.60x | 4d0ab395b00f | 4d0ab395b00f | 9.97 | 10.24 |
| homfam_sdr | 3 | 267.4 | 12.4 | 21.56x | 4d0ab395b00f | 4d0ab395b00f | 9.96 | 10.27 |
| homfam_Acetyltransf | 1 | 185.8 | 18.7 | 9.94x | f9d7882de1bf | f9d7882de1bf | 8.47 | 8.64 |
| homfam_Acetyltransf | 2 | 200.7 | 10.4 | 19.30x | f9d7882de1bf | f9d7882de1bf | 8.47 | 8.69 |
| homfam_Acetyltransf | 3 | 172.4 | 8.9 | 19.37x | f9d7882de1bf | f9d7882de1bf | 8.47 | 8.65 |
| rnasim25k | 1 | 1196.6 | 46.5 | 25.73x | ce1f97d60307 | ce1f97d60307 | 4.00 | 4.51 |
| rnasim25k | 2 | 1198.1 | 31.3 | 38.28x | ce1f97d60307 | ce1f97d60307 | 4.00 | 4.64 |
| rnasim25k | 3 | 1204.8 | 27.0 | 44.62x | ce1f97d60307 | ce1f97d60307 | 4.00 | 4.60 |
| homfam_rvp | 1 | 1160.6 | 46.5 | 24.96x | a6774caadc06 | a6774caadc06 | 32.97 | 33.15 |
| homfam_rvp | 2 | 916.6 | 28.3 | 32.39x | a6774caadc06 | a6774caadc06 | 33.00 | 33.20 |
| homfam_rvp | 3 | 930.1 | 23.9 | 38.92x | a6774caadc06 | a6774caadc06 | 33.00 | 33.20 |
| rnasim50k | 1 | 2491.8 | 106.2 | 23.46x | 47bb47d38946 | 47bb47d38946 | 12.54 | 12.40 |
| rnasim50k | 2 | 2549.4 | 90.2 | 28.26x | 47bb47d38946 | 47bb47d38946 | 12.59 | 12.73 |
| rnasim50k | 3 | 2529.6 | 68.1 | 37.15x | 47bb47d38946 | 47bb47d38946 | 12.58 | 12.77 |

### 3.3 Identity suite (REP = 1, job 13072136)

The suite compares md5s against `baseline-default`, for v9 built with default flags and with
`-O3 -march=native -flto`:

- sim1k, rnasim579, rnasim1k, prot2k and rnasim2k at 1, 8 and 32 threads;
- edge, edge_noshort and 16S.3 at 1 and 8 threads.

That is 42 checks. The suite then runs the self-check modes on rnasim1k, 16S.T and homfam_aldosered at
32 threads, and BAliBASE (386 inputs × 1 and 8 threads).

| part | result |
|---|---|
| md5 checks | **42 checks, 0 differ** |
| self-check modes | **passed**: `sextetfast check passed` on all three sets (kinds 1, 1 and 4, SIMD level 2 = AVX2), every pair's distance bitwise identical; no `FAILED` |
| BAliBASE | **772 runs, 0 differ** |
| timed runs (§3.1/3.2) | **36 of 36** v9 outputs have the same md5 as the baseline's for the same input |

## 4. Earlier CARC measurements (orientation only, not the final evaluation)

These come from development jobs recorded in `notes/v8s-lever.md` and `notes/v9p-lever.md`. They are
all `--thread 32` with gcc 13.3 `-O3 -march=native -flto`. The baselines were not always measured in
the same job or on the same node, so the ratios are approximate.

| set | baseline-aggr (s) | v9 final source (s) | job | ratio |
|---|---|---|---|---|
| homfam_rvp | ~950–1220 (earlier jobs) | 24.3 / 24.3 | 13066113, 2 reps | 39–50x |
| homfam_sdr | ~294–322 (earlier jobs) | 13.4 / 12.7 | 13066113, 2 reps | 22–25x |
| homfam_Acetyltransf | 190.1 (job 13063938) | 9.6 / 10.1 | 13066113, 2 reps | 19–20x |
| homfam_p450 | ~147–157 (earlier jobs) | 6.9 / 9.5 (quiet node a03-12: 6.31–6.67 over 5 runs) | 13066113, 13071504 | 15–23x |
| rnasim10k | 92.9 (v7 job, node a04-08) | 8.3 (1 rep) | 13066113 | ~11x |
| 16S.T | 57.4 (v7 job, node a04-08) | 5.6 (1 rep) | 13066113 | ~10x |
| rnasim20k, 16S.B.ALL, rnasim50k | not measured in these jobs | 22.2, 23.3, 87.6 (1 rep each) | 13066113 | |

- **rnasim5k.** v8 took 4.3 / 4.3 s against the baseline's 43.8 s (job 13030789). v9 was not timed
  on this set in these jobs. On the DNA sets it did time, v9 was unchanged or faster than v8.
- **Peak RSS, v9 against the baseline (GB).**

  | set | v9 | baseline |
  |---|---|---|
  | homfam_rvp | 33.19–33.20 | 32.97 |
  | homfam_sdr | 10.24–10.28 | 9.96 |
  | homfam_Acetyltransf | 8.62–8.68 | 8.47 |
  | homfam_p450 | 2.55–2.59 | 2.46 |

  v9 is 0.7–5% above the baseline on these sets. With v8 (same code in this area), 16S.T was 4.2–4.4
  GB against 3.71 GB.
- **For scale.** The first CARC test of an early version (v3, gcc 13.3 `-O3 -g`) gave 2.9x on
  rnasim5k (45.7 → 15.5 s) and 3.5x on rnasim10k (99.3 → 28.1 s).
