# CARC evaluation results

**Result: on the 12 datasets of the v9 evaluation, at `--thread 32` on a 32-core / 64 GB CARC node,
fast-mafft v10 produced output byte-identical to MAFFT v7.526 in every run and was 10.1x to 43.9x
faster on average; its slowest single run was 9.98x** (rnasim5k). The baseline is MAFFT v7.526 built
with `-O3 -march=native -flto` (its output equals the default build's).

Beyond those 12:

- **Real SARS-CoV-2 collections:** 11.6x (500 genomes) and 11.8x (1,000), identical in every run;
  205x on 2,000 genomes, where MAFFT switches to its memsave mode in the middle of a pass and takes
  63 minutes (one of v10's three outputs there differed from MAFFT's; §3.4).
- **Alignments over 30,000 columns** (simulated 30–50 kb genomes, and the edge set): 14.9x to 72.8x,
  identical in every run.
- **Small inputs:** 7.0–7.7x at 2,000 sequences, 5.0–5.6x at 1,000, and 2.4x on prot2k, where MAFFT
  itself takes 1.5 s.
- **Identity suite** (REP = 1): 45 of 45 md5 checks, BAliBASE 772 of 772, and 21 of 23 on the
  long-genome sets; every self-check passed. These two differences and the one in the timed runs
  are all on inputs whose alignment crosses 30,000 columns in the middle of a pass, where MAFFT's own
  multi-threaded output depends on thread timing (§3.4). At `--thread 1` there is no such race, and
  fast-mafft was identical on every input tested.

Source: six final CARC jobs of [`carc/carc_v10_final.sh`](carc/carc_v10_final.sh) (REP = 1, 2, 3 ×
PART = std, long), 2026-10-18, repository commit `b35964d`. The logs are in [`carc/logs/`](carc/logs/)
(`v10final-*.log`), and [`carc/summarize.py`](carc/summarize.py) regenerates the tables of §3 from
them. The v9 evaluation is kept in §4.

## 1. Methodology

| item | setting |
|---|---|
| Hardware | USC CARC `main` partition, `--constraint=epyc-7513` (AMD EPYC 7513, Zen 3, AVX2, no AVX-512). One node per job: `--nodes=1 --ntasks=1 --cpus-per-task=32 --mem=64G`. A node has 64 cores, so two jobs can share one: REP 1's two jobs ran on a03-19 (cores 0–31 and 32–63), REP 2's on a04-05, REP 3's on b10-02 and b10-03. |
| Compiler | gcc 13.3.0 (`module load gcc/13.3.0`). |
| Baseline (timed) | MAFFT v7.526, built inside each job from the repository's first commit (`f064354`, checked against the md5 of the pristine `disttbfast.c`), with gcc 13.3 and `-O3 -march=native -flto`. Its output equals the default build's. |
| Identity reference | MAFFT v7.526 built with gcc 13.3 and the Makefile's default flags (REP = 1). |
| fast-mafft v10 | Built inside each job from `fast-mafft-src`, with `-O3 -march=native -flto` (timed) and, in REP = 1, with the default flags (identity suite). |
| Command | `bin/mafft --thread 32 in.fa`, run through `carc/carc_run.py`: wall-clock time from launch to exit, the md5 of stdout, peak RSS (`ru_maxrss` of the largest child process, i.e. `disttbfast`) and stage windows from MAFFT's stderr markers. |
| Interleaving | For each input in the order of §3, the baseline runs and then fast-mafft, on the same node in the same job. When the two md5s differ, both run again. |
| Repetitions | Three per input (REP = 1, 2, 3), each in its own job. |
| Per-run speedup | Baseline wall time / fast-mafft wall time within the same job. |
| Mean speedup | Arithmetic mean of the three per-run speedups. |
| Min speedup | The smallest of the three. |
| md5 match | fast-mafft's md5 equals the baseline's in all three repetitions. |
| Peak RSS | The largest of the three runs, in GB. |

## 2. Datasets

The exact input files are published in this repository's `data-v1` release. `bash bench/get_data.sh`
downloads them and verifies the checksums in `bench/DATA_MD5.txt`. Provenance and licenses are in
`bench/DATA.md`. `bash bench/get_data.sh --long` adds the long-genome sets (release `data-long-v1`).
To rerun the comparison on them or on any other inputs, use `carc/run_eval.sh`.

The 12 datasets of the v9 evaluation (lengths are ungapped; the input md5 is the first 12 hex digits of the FASTA file):

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
- **Smaller sets** (timed in §3.2, and in the identity suite): synthetic sim1k, sim2k and prot2k
  (`bench/simulate*.py`), rnasim579 (TWILIGHT's GitHub sample, MIT), and the rnasim1k and rnasim2k
  subsets of RNASim 10K.
- **Long-genome sets** (timed in §3.3; `bash bench/get_data.sh --long`, provenance in
  [bench/DATA.md](bench/DATA.md)):

  | dataset | sequences | length | what it is |
  |---|---:|---|---|
  | sarscov2_500 / 1000 / 2000 | 500 / 1,000 / 2,000 | 29,641 – 29,903 | real complete SARS-CoV-2 genomes from NCBI (881 of the 2,000 with runs of N); the alignments of 500 and 1,000 stay just under 30,000 columns, 2,000 crosses it in the middle of pass 1 |
  | edge | 300 | up to 2,952 | sim1k with long unaligned repeats: 59–66 k columns |
  | long100x30k, long500x30k, long1000x30k | 100 / 500 / 1,000 | ~30,000 | simulated genomes (`simulate.py`), over 30,000 columns from the first steps |
  | long300x50k | 300 | ~50,000 | simulated, 50 kb |

## 3. Final results (v10)

### 3.1 The 12 datasets of the v9 evaluation, `--thread 32`

| dataset | baseline mean (s) | fast-mafft mean (s) | speedup r1 | speedup r2 | speedup r3 | mean speedup | min speedup | md5 match | peak RSS baseline / fast-mafft (GB) |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| rnasim5k | 44.1 | 4.4 | 10.02x | 9.98x | 10.30x | 10.10x | 9.98x | yes (3/3) | 2.93 / 3.35 |
| rnasim10k | 94.1 | 8.1 | 12.05x | 11.19x | 11.80x | 11.68x | 11.19x | yes (3/3) | 5.81 / 3.68 |
| 16S.T | 58.6 | 5.1 | 11.60x | 11.41x | 11.71x | 11.57x | 11.41x | yes (3/3) | 3.69 / 4.02 |
| 16S.3 | 46.7 | 4.0 | 11.92x | 11.39x | 11.72x | 11.68x | 11.39x | yes (3/3) | 3.15 / 3.19 |
| homfam_p450 | 140.7 | 6.5 | 21.43x | 21.15x | 22.71x | 21.77x | 21.15x | yes (3/3) | 2.46 / 2.59 |
| rnasim20k | 231.9 | 17.3 | 13.38x | 13.39x | 13.37x | 13.38x | 13.37x | yes (3/3) | 12.93 / 5.19 |
| 16S.B.ALL | 505.1 | 17.6 | 28.86x | 28.26x | 28.99x | 28.70x | 28.26x | yes (3/3) | 3.84 / 4.46 |
| homfam_sdr | 269.1 | 12.1 | 22.51x | 22.63x | 21.44x | 22.19x | 21.44x | yes (3/3) | 9.96 / 10.26 |
| homfam_Acetyltransf | 172.9 | 9.0 | 19.45x | 19.20x | 19.00x | 19.22x | 19.00x | yes (3/3) | 8.47 / 8.66 |
| rnasim25k | 1189.3 | 27.1 | 42.57x | 43.96x | 45.02x | 43.85x | 42.57x | yes (3/3) | 4.00 / 4.51 |
| homfam_rvp | 892.6 | 25.0 | 33.66x | 40.52x | 33.95x | 36.04x | 33.66x | yes (3/3) | 33.00 / 33.20 |
| rnasim50k | 2517.5 | 72.1 | 35.58x | 33.93x | 35.34x | 34.95x | 33.93x | yes (3/3) | 12.59 / 12.94 |

All 36 timed runs gave the baseline's md5. Against v9 (§4.1), v10's mean speedup is higher on every
dataset, and its slowest run rose from 8.3x to 9.98x: the anchor and segment changes of CHANGES.md
§8.5 also apply to these sets.

### 3.2 Smaller inputs (1,000–2,000 sequences)

| dataset | baseline mean (s) | fast-mafft mean (s) | speedup r1 | speedup r2 | speedup r3 | mean speedup | min speedup | md5 match | peak RSS baseline / fast-mafft (GB) |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| sim1k | 5.7 | 1.1 | 4.83x | 5.09x | 5.09x | 5.01x | 4.83x | yes (3/3) | 0.69 / 1.94 |
| rnasim1k | 10.4 | 1.9 | 5.20x | 5.78x | 5.83x | 5.60x | 5.20x | yes (3/3) | 0.88 / 2.32 |
| prot2k | 1.5 | 0.6 | 2.67x | 2.14x | 2.50x | 2.44x | 2.14x | yes (3/3) | 0.32 / 0.84 |
| sim2k | 11.9 | 1.7 | 7.06x | 7.06x | 6.94x | 7.02x | 6.94x | yes (3/3) | 1.02 / 2.43 |
| rnasim2k | 19.8 | 2.6 | 7.30x | 7.92x | 8.00x | 7.74x | 7.30x | yes (3/3) | 1.80 / 2.54 |

At this size MAFFT itself finishes in 1.5–20 s, and the fixed costs (launcher, input, the serial
parts of each stage) are a larger share. Every run was identical.

### 3.3 Long genomes and SARS-CoV-2 collections

| dataset | baseline mean (s) | fast-mafft mean (s) | speedup r1 | speedup r2 | speedup r3 | mean speedup | min speedup | md5 match | peak RSS baseline / fast-mafft (GB) |
|---|---:|---:|---:|---:|---:|---:|---:|---|---|
| sarscov2_500 | 28.7 | 2.5 | 11.48x | 11.44x | 11.96x | 11.63x | 11.44x | yes (3/3) | 0.80 / 1.57 |
| sarscov2_1000 | 48.6 | 4.1 | 12.46x | 11.57x | 11.30x | 11.78x | 11.30x | yes (3/3) | 1.37 / 2.00 |
| sarscov2_2000 | 3810.4 | 18.6 | 201.07x | 204.89x | 209.90x | 205.29x | 201.07x | **no** (2/3), see §3.4 | 3.00 / 2.78 |
| edge | 79.8 | 5.4 | 14.69x | 15.37x | 14.56x | 14.87x | 14.56x | yes (3/3) | 6.11 / 5.64 |
| long100x30k | 326.3 | 9.0 | 37.28x | 37.00x | 34.96x | 36.41x | 34.96x | yes (3/3) | 0.44 / 0.78 |
| long1000x30k | 939.7 | 23.1 | 41.65x | 40.76x | 39.85x | 40.75x | 39.85x | yes (3/3) | 0.96 / 1.49 |
| long300x50k | 2060.5 | 28.5 | 78.49x | 76.54x | 63.50x | 72.84x | 63.50x | yes (3/3) | 0.68 / 1.50 |
| long500x30k | 1633.1 | 39.8 | 43.60x | 41.30x | 38.58x | 41.16x | 38.58x | yes (3/3) | 0.81 / 1.48 |

The baseline's 63 minutes on sarscov2_2000 are MAFFT's memsave mode: the alignment crosses 30,000
columns in the middle of pass 1, and every later merge of the pass, and of pass 2, runs MAFFT's
single-threaded Myers-Miller aligner. Peak memory: at most 1.5 GB on the simulated long sets (MAFFT:
0.4–1.0 GB), 1.6–2.8 GB on the SARS-CoV-2 sets (MAFFT: 0.8–3.0 GB).

### 3.4 Identity

**Identity suite** (REP = 1; reference = v7.526 built with the default flags; fast-mafft built with
the default flags and with `-O3 -march=native -flto`, both compared):

| part | what | result |
|---|---|---|
| std | sim1k, rnasim579, rnasim1k, prot2k, rnasim2k at 1, 8 and 32 threads; edge, edge_noshort and 16S.3 at 1 and 8 | **42 md5 checks, 0 differ** |
| std | self-checks at 32 threads on rnasim1k, 16S.T, homfam_aldosered (6-mer distances, MSA distances, guide trees, `blockAlign2`, anchors and lags, long FFTs) | **3 runs, all passed**, md5 identical |
| std | BAliBASE 3.0, 386 inputs at 1 and 8 threads | **772 runs, 0 differ** |
| long | strict one-thread identity on sarscov2_500 and sarscov2_1000 | **identical** (`e8be67d631bc`, `3b000a61e12a`) |
| long | long10x30k, long30x30k at 1 and 8 threads | **identical** |
| long | viral100 at 1 and 8 threads | identical at 1 thread; **different at 8** (see below) |
| long | self-checks: sarscov2_500 and sarscov2_1000 at 8 threads, sarscov2_500 at 32 (`blockAlign2`, long FFTs, anchors); long30x30k at 32 with every check (3,334 Myers-Miller calls and 22,012 recursion nodes verified) | **all passed**, md5 identical |
| long | long10x30k at one thread at each SIMD level (scalar, SSE2, AVX2) with `MAFFT_MSCHECK` | **identical**, checks passed |

**Timed runs:** 74 of 75 identical. The exception is sarscov2_2000, REP 1: the baseline gave
`e6734daa2634` twice, fast-mafft `300dbffd1556` and then, on the rerun, `e6734daa2634`.

**Inputs whose alignment crosses 30,000 columns in the middle of a pass** (viral100, sarscov2_2000).
There MAFFT switches to its memsave aligners while other steps are in flight, and which of those steps,
and which of their segments, still use the normal aligner depends on thread timing (CHANGES.md §8.4).
MAFFT's own output then depends on timing too. A diagnostic job (13226265, [`carc/carc_v10_mix.sh`](carc/carc_v10_mix.sh),
node a04-16; builds that log every step's aligner, stdout unchanged) measured both programs on the
same node:

| input | threads | MAFFT v7.526 | fast-mafft v10 |
|---|---:|---|---|
| viral100 | 8 | `1e80180224db` in 8 of 8 runs | `1e80180224db` in 6 of 8; `4f9775f500b2` (the one-thread output) in 2 |
| viral100 | 32 | `96dcc1cd7d39` in 3 of 4, `f47ba74470f8` in 1 | 4 different md5s in 4 runs, one of them `96dcc1cd7d39` |
| sarscov2_2000 | 32 | `e6734daa2634` in 5 of 5 (with the timed runs) | `e6734daa2634` in 11 of 16; `300dbffd1556` in 5 |

In MAFFT, the steps near the switch take about the same time from run to run, so one outcome of the
race usually wins on a given machine (but not always: viral100 at 32 threads). fast-mafft runs those
steps much faster (20x to 200x overall on these two inputs), by different factors for different steps, so the race has different odds
and sometimes a different winner. All 37 logged runs of both programs switched at the same steps (49 and 50; 1884 and 1981);
the runs differ only in which steps near the switch, and which of their segments, used which aligner.
At `--thread 1` there is no race, and fast-mafft was identical on every input tested. On every input whose
alignment stays under 30,000 columns, or crosses it in the first steps of a pass (the long simulated
sets, edge), the output was identical in every run.

### 3.5 Raw runs

| dataset | rep | baseline wall (s) | fast-mafft wall (s) | speedup | baseline md5 | fast-mafft md5 | baseline RSS (GB) | fast-mafft RSS (GB) |
|---|---:|---:|---:|---:|---|---|---:|---:|
| sim1k | 1 | 5.8 | 1.2 | 4.83x | 85a82397cf39 | 85a82397cf39 | 0.67 | 1.93 |
| sim1k | 2 | 5.6 | 1.1 | 5.09x | 85a82397cf39 | 85a82397cf39 | 0.69 | 1.92 |
| sim1k | 3 | 5.6 | 1.1 | 5.09x | 85a82397cf39 | 85a82397cf39 | 0.65 | 1.94 |
| rnasim1k | 1 | 10.4 | 2.0 | 5.20x | 4ec80fa72d83 | 4ec80fa72d83 | 0.87 | 2.32 |
| rnasim1k | 2 | 10.4 | 1.8 | 5.78x | 4ec80fa72d83 | 4ec80fa72d83 | 0.88 | 2.20 |
| rnasim1k | 3 | 10.5 | 1.8 | 5.83x | 4ec80fa72d83 | 4ec80fa72d83 | 0.88 | 2.27 |
| prot2k | 1 | 1.6 | 0.6 | 2.67x | 54420c021453 | 54420c021453 | 0.32 | 0.82 |
| prot2k | 2 | 1.5 | 0.7 | 2.14x | 54420c021453 | 54420c021453 | 0.30 | 0.84 |
| prot2k | 3 | 1.5 | 0.6 | 2.50x | 54420c021453 | 54420c021453 | 0.32 | 0.82 |
| sim2k | 1 | 12.0 | 1.7 | 7.06x | f9b5e4a733ca | f9b5e4a733ca | 0.97 | 2.39 |
| sim2k | 2 | 12.0 | 1.7 | 7.06x | f9b5e4a733ca | f9b5e4a733ca | 0.97 | 2.40 |
| sim2k | 3 | 11.8 | 1.7 | 6.94x | f9b5e4a733ca | f9b5e4a733ca | 1.02 | 2.43 |
| rnasim2k | 1 | 19.7 | 2.7 | 7.30x | af6ae72412b3 | af6ae72412b3 | 1.76 | 2.54 |
| rnasim2k | 2 | 19.8 | 2.5 | 7.92x | af6ae72412b3 | af6ae72412b3 | 1.78 | 2.40 |
| rnasim2k | 3 | 20.0 | 2.5 | 8.00x | af6ae72412b3 | af6ae72412b3 | 1.80 | 2.40 |
| rnasim5k | 1 | 44.1 | 4.4 | 10.02x | 455da94efd24 | 455da94efd24 | 2.85 | 3.13 |
| rnasim5k | 2 | 43.9 | 4.4 | 9.98x | 455da94efd24 | 455da94efd24 | 2.91 | 3.35 |
| rnasim5k | 3 | 44.3 | 4.3 | 10.30x | 455da94efd24 | 455da94efd24 | 2.93 | 3.10 |
| rnasim10k | 1 | 94.0 | 7.8 | 12.05x | 5e7bb28c8e8e | 5e7bb28c8e8e | 5.75 | 3.68 |
| rnasim10k | 2 | 94.0 | 8.4 | 11.19x | 5e7bb28c8e8e | 5e7bb28c8e8e | 5.57 | 3.61 |
| rnasim10k | 3 | 94.4 | 8.0 | 11.80x | 5e7bb28c8e8e | 5e7bb28c8e8e | 5.81 | 3.66 |
| 16S.T | 1 | 58.0 | 5.0 | 11.60x | f487c1a56702 | f487c1a56702 | 3.64 | 4.02 |
| 16S.T | 2 | 58.2 | 5.1 | 11.41x | f487c1a56702 | f487c1a56702 | 3.69 | 3.93 |
| 16S.T | 3 | 59.7 | 5.1 | 11.71x | f487c1a56702 | f487c1a56702 | 3.62 | 4.02 |
| 16S.3 | 1 | 46.5 | 3.9 | 11.92x | 6b646c6b3c1f | 6b646c6b3c1f | 3.15 | 2.95 |
| 16S.3 | 2 | 46.7 | 4.1 | 11.39x | 6b646c6b3c1f | 6b646c6b3c1f | 3.11 | 3.19 |
| 16S.3 | 3 | 46.9 | 4.0 | 11.72x | 6b646c6b3c1f | 6b646c6b3c1f | 3.07 | 2.96 |
| homfam_p450 | 1 | 139.3 | 6.5 | 21.43x | a8da25a0a301 | a8da25a0a301 | 2.46 | 2.59 |
| homfam_p450 | 2 | 139.6 | 6.6 | 21.15x | a8da25a0a301 | a8da25a0a301 | 2.46 | 2.59 |
| homfam_p450 | 3 | 143.1 | 6.3 | 22.71x | a8da25a0a301 | a8da25a0a301 | 2.46 | 2.59 |
| rnasim20k | 1 | 230.1 | 17.2 | 13.38x | 89f21e3ddf15 | 89f21e3ddf15 | 12.30 | 4.91 |
| rnasim20k | 2 | 230.3 | 17.2 | 13.39x | 89f21e3ddf15 | 89f21e3ddf15 | 12.61 | 5.19 |
| rnasim20k | 3 | 235.3 | 17.6 | 13.37x | 89f21e3ddf15 | 89f21e3ddf15 | 12.93 | 4.76 |
| 16S.B.ALL | 1 | 505.0 | 17.5 | 28.86x | 8a066cfcbfbe | 8a066cfcbfbe | 3.83 | 4.41 |
| 16S.B.ALL | 2 | 503.0 | 17.8 | 28.26x | 8a066cfcbfbe | 8a066cfcbfbe | 3.84 | 4.46 |
| 16S.B.ALL | 3 | 507.4 | 17.5 | 28.99x | 8a066cfcbfbe | 8a066cfcbfbe | 3.84 | 4.41 |
| homfam_sdr | 1 | 267.9 | 11.9 | 22.51x | 4d0ab395b00f | 4d0ab395b00f | 9.96 | 10.26 |
| homfam_sdr | 2 | 269.3 | 11.9 | 22.63x | 4d0ab395b00f | 4d0ab395b00f | 9.96 | 10.23 |
| homfam_sdr | 3 | 270.1 | 12.6 | 21.44x | 4d0ab395b00f | 4d0ab395b00f | 9.96 | 10.25 |
| homfam_Acetyltransf | 1 | 173.1 | 8.9 | 19.45x | f9d7882de1bf | f9d7882de1bf | 8.47 | 8.65 |
| homfam_Acetyltransf | 2 | 172.8 | 9.0 | 19.20x | f9d7882de1bf | f9d7882de1bf | 8.47 | 8.66 |
| homfam_Acetyltransf | 3 | 172.9 | 9.1 | 19.00x | f9d7882de1bf | f9d7882de1bf | 8.47 | 8.66 |
| rnasim25k | 1 | 1183.5 | 27.8 | 42.57x | ce1f97d60307 | ce1f97d60307 | 4.00 | 4.51 |
| rnasim25k | 2 | 1186.9 | 27.0 | 43.96x | ce1f97d60307 | ce1f97d60307 | 4.00 | 4.46 |
| rnasim25k | 3 | 1197.6 | 26.6 | 45.02x | ce1f97d60307 | ce1f97d60307 | 4.00 | 4.48 |
| homfam_rvp | 1 | 878.5 | 26.1 | 33.66x | a6774caadc06 | a6774caadc06 | 33.00 | 33.17 |
| homfam_rvp | 2 | 879.2 | 21.7 | 40.52x | a6774caadc06 | a6774caadc06 | 32.99 | 33.19 |
| homfam_rvp | 3 | 920.0 | 27.1 | 33.95x | a6774caadc06 | a6774caadc06 | 32.97 | 33.20 |
| rnasim50k | 1 | 2512.2 | 70.6 | 35.58x | 47bb47d38946 | 47bb47d38946 | 12.59 | 12.94 |
| rnasim50k | 2 | 2520.7 | 74.3 | 33.93x | 47bb47d38946 | 47bb47d38946 | 12.59 | 12.47 |
| rnasim50k | 3 | 2519.6 | 71.3 | 35.34x | 47bb47d38946 | 47bb47d38946 | 12.55 | 12.48 |
| sarscov2_500 | 1 | 28.7 | 2.5 | 11.48x | e8be67d631bc | e8be67d631bc | 0.78 | 1.55 |
| sarscov2_500 | 2 | 28.6 | 2.5 | 11.44x | e8be67d631bc | e8be67d631bc | 0.80 | 1.57 |
| sarscov2_500 | 3 | 28.7 | 2.4 | 11.96x | e8be67d631bc | e8be67d631bc | 0.71 | 1.55 |
| sarscov2_1000 | 1 | 48.6 | 3.9 | 12.46x | 3b000a61e12a | 3b000a61e12a | 1.13 | 1.94 |
| sarscov2_1000 | 2 | 48.6 | 4.2 | 11.57x | 3b000a61e12a | 3b000a61e12a | 1.15 | 2.00 |
| sarscov2_1000 | 3 | 48.6 | 4.3 | 11.30x | 3b000a61e12a | 3b000a61e12a | 1.37 | 1.97 |
| sarscov2_2000 | 1 | 3800.2 | 18.9 | 201.07x | e6734daa2634 | 300dbffd1556 | 3.00 | 2.78 |
| sarscov2_2000 | 2 | 3810.9 | 18.6 | 204.89x | e6734daa2634 | e6734daa2634 | 2.97 | 2.62 |
| sarscov2_2000 | 3 | 3820.2 | 18.2 | 209.90x | e6734daa2634 | e6734daa2634 | 2.87 | 2.71 |
| edge | 1 | 79.3 | 5.4 | 14.69x | cc0dbbde5397 | cc0dbbde5397 | 6.11 | 5.22 |
| edge | 2 | 79.9 | 5.2 | 15.37x | cc0dbbde5397 | cc0dbbde5397 | 5.84 | 5.49 |
| edge | 3 | 80.1 | 5.5 | 14.56x | cc0dbbde5397 | cc0dbbde5397 | 5.98 | 5.64 |
| long100x30k | 1 | 328.1 | 8.8 | 37.28x | 337e5ebd176c | 337e5ebd176c | 0.39 | 0.78 |
| long100x30k | 2 | 325.6 | 8.8 | 37.00x | 337e5ebd176c | 337e5ebd176c | 0.41 | 0.75 |
| long100x30k | 3 | 325.1 | 9.3 | 34.96x | 337e5ebd176c | 337e5ebd176c | 0.44 | 0.76 |
| long1000x30k | 1 | 937.2 | 22.5 | 41.65x | 25999048f8b2 | 25999048f8b2 | 0.90 | 1.49 |
| long1000x30k | 2 | 941.6 | 23.1 | 40.76x | 25999048f8b2 | 25999048f8b2 | 0.96 | 1.45 |
| long1000x30k | 3 | 940.4 | 23.6 | 39.85x | 25999048f8b2 | 25999048f8b2 | 0.93 | 1.33 |
| long300x50k | 1 | 2056.4 | 26.2 | 78.49x | edb16ecc771b | edb16ecc771b | 0.62 | 1.49 |
| long300x50k | 2 | 2074.3 | 27.1 | 76.54x | edb16ecc771b | edb16ecc771b | 0.68 | 1.47 |
| long300x50k | 3 | 2050.9 | 32.3 | 63.50x | edb16ecc771b | edb16ecc771b | 0.57 | 1.50 |
| long500x30k | 1 | 1635.1 | 37.5 | 43.60x | 2c2158b4d42d | 2c2158b4d42d | 0.77 | 1.47 |
| long500x30k | 2 | 1639.8 | 39.7 | 41.30x | 2c2158b4d42d | 2c2158b4d42d | 0.81 | 1.48 |
| long500x30k | 3 | 1624.3 | 42.1 | 38.58x | 2c2158b4d42d | 2c2158b4d42d | 0.78 | 1.41 |

Rerun after the md5 difference:

| dataset | rep | baseline md5, run 1 / run 2 | fast-mafft md5, run 1 / run 2 |
|---|---:|---|---|
| sarscov2_2000 | 1 | e6734daa2634 / e6734daa2634 | 300dbffd1556 / e6734daa2634 |

## 4. The v9 evaluation (2026-10-13), for comparison

The previous version, measured with [`carc/carc_v9_final.sh`](carc/carc_v9_final.sh) (logs `carc/logs/v9final-*.log`)
on the same 12 datasets. v10 is faster on every one of them (§3.1 against §4.1).

### 4.1 Summary, `--thread 32`

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

### 4.2 Raw runs

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

### 4.3 Identity suite (REP = 1, job 13072136)

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
| timed runs (§4.1/4.2) | **36 of 36** v9 outputs have the same md5 as the baseline's for the same input |

## 5. Earlier CARC measurements of v9 (orientation only)

These come from development jobs recorded in the v8s and v9p development notes (not published). They are
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
