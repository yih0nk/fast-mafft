# How byte-identical output was verified

Every count below comes from the lever notes (`notes/*.md`) or from the final packaging checks run
for this deliverable (§3). "M5" is the development Mac (Apple M5); "CARC" is USC CARC (AMD EPYC
7513/7542, gcc 13.3).

## 1. What "identical" means here

- **The criterion.** The md5 of stdout of `mafft --thread N in.fa` (default FFT-NS-2) must equal
  the md5 from MAFFT v7.526 built from `baseline-src` (commit `ee97999`) with the same compiler and
  flags, run at the same `N`.
- **Why the compiler matters.** Floating-point contraction depends on the compiler and flags, so
  there are two md5 "classes" for the same input. Both baseline and fast-mafft fall into the same
  class for a given build:
  - **FMA-contracting builds:** Apple clang on arm64, and clang with `-march=x86-64-v3` or `native`.
  - **Non-contracting builds:** gcc with any flags (`-std=c99` gives `-ffp-contract=off`), and clang
    targeting plain x86-64.
- **On CARC.** The reference is v7.526 built with gcc 13.3 and the Makefile's default flags. v9 built
  with the default flags and v9 built with `-O3 -march=native -flto` must both match it. v7.526's
  own aggressive build gives the same md5 as its default build.
- **Two places where v7.526 itself is not deterministic,** and how they were handled:
  - **The memsave switch.** When an alignment exceeds 30,000 columns with at most 20,000 sequences
    and more than one thread, the result depends on the schedule. Here fast-mafft must produce
    outputs from the same set (CHANGES.md §3.6). At `--thread 1`, and above 20,000 sequences, it
    must match exactly.
  - **The iterative strategies at `--thread` > 1.** FFT-NS-i, L-INS-i and `--auto` on small inputs
    vary from run to run in v7.526 itself, so they are compared at `--thread 1` only.

## 2. Datasets

| group | sets | sequences | what they exercise |
|---|---|---|---|
| synthetic DNA (`bench/simulate.py`, tree-structured) | sim1k, sim2k, sim5k, sim10k | 1,000 to 10,000 × 1.5 kb | the DNA paths at small and medium N |
| | sim3k_long | 3,000 × 5 kb | long sequences; the memsave switch fires in pass 2 (step 2995 of 2999) |
| synthetic protein (`bench/simulate_protein.py`) | prot2k | 2,000 × 369–446 aa | protein kernels: 6-mer list, inverted index and kind 4; msadist column path; general DP rows |
| edge cases | edge.fa, edge_noshort.fa | sim1k + long unaligned repeats | 59–66 k alignment columns: the memsave switch (`alg = 'M'`) |
| | edge2.fa | 400 | 6-mer counts ≥ 255 (saturation correction), poly-A/(ac)n/g runs, short sequences; 8,559 columns |
| RNASim | rnasim579 | 579 | TWILIGHT's GitHub sample (real RNASim) |
| | rnasim1k, rnasim2k, rnasim5k, rnasim10k | 1,000 to 10,000 | subsets of RNASim 10K (seeds 1, 2, 5) and RNASim 10K itself, from the Recursive MAGUS archive |
| | rnasim20k, rnasim25k, rnasim50k | 20,000 to 50,000 | RNASim 20K and 50K from the same archive; rnasim25k = 25,000 of rnasim50k (seed 25). Above 20,000: the threaded-pass rewrite |
| 16S rRNA (same archive) | 16S.T, 16S.3, 16S.B.ALL | 7,350, 6,323, 27,643 | real rRNA with IUPAC codes: msadist exceptions, general DP rows; 16S.B.ALL is above 20,000 |
| | subsets 16ST600, 16SB800, 16sT800, 16sB1000 | 600 to 1,000 | fast IUPAC checks, including under Rosetta |
| HomFam protein families (same archive) | homfam_aldosered, homfam_p450, homfam_sdr, homfam_Acetyltransf, homfam_rvp | 13,277 to 93,681 | large-N protein: inverted index, kind 4 masks, column path with dense columns, v9 guide tree |
| | subsets p450_1k, rvp2k, ald1500, homfam_rvp_sub5k/20k/30k, homfam_sdr_sub3k/20k, homfam_Acetyltransf_sub20k | 1,000 to 30,000 | the same paths at sizes that run locally or quickly |
| BAliBASE 3.0 core (R1–5) | 386 `.tfa` inputs, RV11–RV50 | small (about 29 on average per reference case) | diverse small real protein inputs: short sequences, `X`, extensions, insertions |

## 3. Checks on the final packaged tree (this deliverable, 2026-10-13)

- **The patch applies cleanly.**
  - Base: `git archive ee97999` of `baseline-src`, which is clean v7.526 without the build artifacts
    present in `baseline-src`'s working tree.
  - `patch -p1 < fast-mafft-v9.patch` applied with no fuzz, offsets or rejects.
  - `git apply` applied too. It warns about trailing whitespace on 9 added lines, which keep MAFFT's
    own style.
- **Both results equal `fast-mafft-src/` exactly.** `diff -r` reports no difference over all 179
  files, and the permission bits are equal as well.
- **`fast-mafft-src/` equals `work/v9p/src`.** It holds every tracked file plus the 7 new sources,
  with no build artifacts. Its `disttbfast.c` md5 is `5a2011f1…` and `mltaln9.c` is `6214cdc5…`, the
  files the final CARC jobs built (`notes/v9p-lever.md`).
- **Builds.** The patched tree builds with Apple clang 15 (`make CC=clang PREFIX=…`, default
  `CFLAGS=-O3`). It also builds with the literal stock form `make -C core && make -C core install
  PREFIX=…`, whose launcher points at `/usr/local/libexec/mafft` (BUILD.md §2). Neither build adds a
  new compiler warning category to any `.c` file compared with the baseline's build log.
- **md5 at `--thread 4`, against `baseline-install`** (Apple clang, the reference build):

  | set | baseline-install | patched tree, clang | stock-make build |
  |---|---|---|---|
  | sim1k | 3e271a550159 | 3e271a550159 | 3e271a550159 |
  | rnasim1k | d3d860c0602f | d3d860c0602f | d3d860c0602f |
  | prot2k | 087f57deb475 | 087f57deb475 | 087f57deb475 |

- **Self-check modes.** `MAFFT_SEXTET=check MAFFT_MSADIST=check MAFFT_UPG=check` on rnasim1k
  (6-mer kind 1, msadist bit-plane path) and prot2k (kind 4, column path) at `--thread 4` passed,
  with the md5 unchanged.
- **Final CARC suite.** REP=1 of `carc/carc_v9_final.sh` runs 42 md5 checks, the self-check modes
  and BAliBASE at 1 and 8 threads. Its result goes into RESULTS.md.

## 4. Checks on the final v9 source (v9p lever)

**CARC, job 13066113** (node a03-17). Reference: v7.526 built with gcc 13.3 default flags. Candidates:
v9 built with default flags and with `-O3 -march=native -flto`. **88 md5 checks, 0 differ.**

| sets | threads |
|---|---|
| sim1k, rnasim579, rnasim1k, prot2k, rnasim2k, homfam_sdr_sub3k, homfam_rvp_sub5k | 1, 8, 32 |
| edge, edge_noshort (each against the original's own md5 at that thread count) | 1, 8 |
| 16S.3 | 8 |
| 16S.T, homfam_aldosered, homfam_sdr_sub20k, homfam_rvp_sub20k, homfam_Acetyltransf_sub20k | 32 |
| full homfam_rvp (a6774caadc06), homfam_sdr (4d0ab395b00f), homfam_p450 (a8da25a0a301), homfam_Acetyltransf (f9d7882de1bf) | 32 |
| `MAFFT_SIMD=scalar\|sse2\|avx2` on rvp_sub5k, sdr_sub3k, prot2k | 4, 32 |

Further checks on CARC:

- **Every timed run's md5 equalled the baseline's.**
- **Self-checks at 32 threads.** `MAFFT_SEXTET=check MAFFT_MSADIST=check,forceskip MAFFT_UPG=check`
  on sdr_sub20k (kind 4), rvp_sub20k (kind 4), Acetyltransf_sub20k (kind 3) and rnasim1k (kind 1).
  `MAFFT_UPG=check` on the full sdr and Acetyltransf. All passed.
- **Job 13065296.** `MAFFT_MSADIST=check,forceskip` at `MAFFT_SIMD=avx2|sse2|scalar` on sdr_sub3k,
  rvp_sub5k, rnasim1k and 16S.T, md5 equal to the baseline's. `MAFFT_UPG=check` on rnasim20k.
- **Job 13063938** (the build just before the thread-pool rework): the same suite, 86 checks, 0
  differ.
- **BAliBASE**, 386 inputs at 1 and 8 threads (default-flags build): **772 runs, 0 differ.**

**Locally (M5), final source:**

- **Identity suite: 29 checks, 0 differ.**
  - rvp_sub20k, rvp_sub30k, sdr_sub20k, Acetyltransf_sub20k, homfam_p450, homfam_aldosered, prot2k,
    rvp_sub5k and sdr_sub3k at 1 and 4 threads.
  - rnasim1k, edge, edge_noshort and 16S.T at 4 threads.
  - The self-check modes.
- **BAliBASE:** 772 runs, 0 differ.
- **Rosetta x86** (`clang -arch x86_64 -O3 -g`, against the same-flag baseline): prot2k, rvp_sub5k,
  sdr_sub3k and rnasim1k at 1 and 4 threads, each at three SIMD settings: `MAFFT_SIMD=scalar`, the
  default (SSE2/POPCNT), and AVX2 (`ROSETTA_ADVERTISE_AVX=1`). **24 runs, 0 differ.** All three
  self-checks also passed at each level.
- **gcc-15 arm64** against the gcc-15 baseline: **10 runs, 0 differ**.
- **Stress builds** with tiny guide-tree parameters (`-DU9_PARMIN=64 -DU9_K=2 -DU9_PD=3
  -DU9_CHUNKMIN=16 -DU9_SYNCF=0.5`, and `-DU9_PARMIN=32 -DU9_K=1`). These make the candidate lists
  fail often and run the parallel phases on small inputs. `MAFFT_UPG=check` on prot2k, sim1k,
  rvp_sub5k, sdr_sub3k and rnasim1k at 1, 3, 5, 8 and 10 threads all passed.

## 5. Checks during development (each lever, against the same-flag baseline)

| lever (notes file) | md5 checks | other checks |
|---|---|---|
| **v3 on CARC** (PLAN.md) | everything run on CARC, with gcc 13.3 `-O3 -g` at `--thread 32`, BAliBASE and the edge sets included: identical | |
| **dist** (distance-lever) | sim1k, sim2k, sim5k, prot2k at 1/4/10 (12 cells); sim10k at 10 | `MAFFT_MSADIST=check` on sim1k, sim2k, sim5k (12.5 M pairs) and prot2k, both paths. sim5k matrix dump identical to the original code's (99,980,000 bytes). Random stress alignments (L from 1 to 4,097, 5 alphabets, 5–97% gaps), including x86 under Rosetta: 0 mismatches |
| **prog** (progressive-lever) | sim1k, sim2k, sim5k, prot2k at 1/4/10; sim10k at 10; sim1k at 2 and 16 | `MAFFT_DPCHECK`: 2.4 M A__align + 0.86 M G__align11 rows on sim1k, 11.8 M + 4.3 M on sim5k. Options `--ep 0.123`, `--op 2.0 --ep 0.05`, `--retree 1`, `--nomemsave`, `--thread 7`: identical |
| **sextet** (sextet-lever) | 8 sets × 1/4/10 = 24 runs; sim3k_long and sim10k at 10 | BAliBASE 772/0. `MAFFT_SEXTET=check` (all modes) on sim1k, sim5k, sim3k_long, sim10k, prot2k, rnasim5k, edge, edge2. Matrix dumps identical |
| **overhead** (overhead-lever) | 5 sets × 1/4/10 = 15; sim3k_long and sim10k at 10 | 72 option/locale combinations (CRLF, Latin-1/UTF-8 names, three locales, stdin, `--auto`, `--localpair`, `--globalpair`, `--maxiterate 2`, `--clustalout`, `--add`, `--treeout`, …). `f2cl` fuzzing on 1,000 + 500 inputs. Guide-tree harness: trees bitwise equal on sim1k, prot2k, sim5k and sim10k dumps at 1, 2, 3, 4, 7, 10 threads and on tie-heavy synthetic matrices. BAliBASE at `--thread 1` and 4, plus the `--auto --ep 0.123 --anysymbol` flags at `--thread 1`: 386/386 each |
| **prog2** (progressive-lever-2) | 8 sets × 1/4/10 = 24; sim3k_long, sim10k, rnasim10k at 10; edge sets at 1/4 (and at 10 as a set of outputs) | BAliBASE 772/0. `MAFFT_TILECHECK` (forced to 16-column bands: 3,727 DPs per run), `DPCHECK`, `PROFCHECK`, `ARCHECK`: no mismatch |
| **gccpar** (gccpar-lever) | clang 18/18 and gcc-15 18/18 (6 sets × 1/4/10); rnasim10k at 10, both compilers | BAliBASE 772/772 with each compiler. Non-default penalties 4/4. gcc-15 `-O3 -mcpu=native -flto`: baseline, v3 and gccpar equal. Rosetta x86 `-O3 -g` 4/4 and `-march=x86-64-v3` 4/4. `DPCHECK`, and `TILECHECK` with forced tiling (1,476 tiled DPs, 3.9 G cells per run), with both compilers |
| **x86simd** (x86simd-lever) | Rosetta x86 matrix: 95 candidate runs at every dispatch level, full 16S.T included, 0 differ. Native 18 runs | Per-pair `MAFFT_SEXTET=check` + `MAFFT_MSADIST=check` at the scalar, SSE2/POPCNT and AVX2 levels and in the x86-64-v3 build: 0 failures. BAliBASE 772/0. Full 16S.T, p450 and aldosered. gcc-15 matrix: 0 differ |
| **prog3** (progressive-lever-3) | 6 sets × 1/4/10; rnasim10k, 16S.3, rnasim20k, aldosered at 10; 16S.T at 4/10; p450, 16S.B.ALL at 4; edge sets at 1/4; Rosetta x86 (two builds) and gcc-15 (two flag sets) on rnasim1k and prot2k | BAliBASE 772/0. `MAFFT_VECCHECK` (2.6–3.9 G cells on rnasim1k, both lane widths, also under Rosetta), `DPCHECK`, forced `TILECHECK`, `LAZYCHECK` (171 deferred tables on rnasim5k, 172 tables and 9 materializations on full 16S.T), `PSEGCHECK` |
| **largen** (largen-lever) | above 20,000: p450, 16S.B.ALL, rnasim25k at 4 and 10 (6/6); at or below 20,000: 14/14 | BAliBASE 772/0. Memsave phases forced at any size (`MAFFT_LARGEN_ORIGSERIAL=1`) on edge, edge_noshort, sim3k_long at 4/10: 6/6. Heavy-step budget stressed: 5/5 (p450 with 5,433 put-backs) |
| **v7 merge** (v7-merge) | 18 + 10 + 3 standard runs; above 20,000 at 4 and 10 (6); `--thread 32` oversubscribed on 8 sets plus the edge sets; gcc-15 8; gcc-15 aggressive 6 + 6; Rosetta x86 18 + 6 | BAliBASE 772/0. Memsave phases with lazy gaps 4/4 (+2 on sim3k_long). Heavy-step stress 4 + 3. Every self-check on rnasim1k and 16S.T (11 runs). Rosetta per-pair checks at three levels (6) |
| **v8s** (v8s-lever) | local 22; above 20,000: p450 at 4/10, 16S.B.ALL and rnasim20k at 10; gcc-15 8. **CARC job 13030789: 42 checks, 0 differ** (1/8/32 threads, both build flavours) | BAliBASE 772/0. All timing-build checks at once on 16S.T at 10 threads (`V8CHECK` covered 684 anchor-profile builds and 1,398 in-place segments). Forced tiling with `TILECHECK`. Stress: rnasim1k at 10 threads, 50 + 60 repeated runs. On CARC: every self-check at 32 threads (14.7 G vector cells verified on 16S.T) |

The 386-input BAliBASE suite was run by the sextet, overhead, prog2, gccpar (once per compiler),
x86simd, prog3, largen, v7, v8s and v9p levers, and on CARC. Each run covered two thread counts, and
none found a difference.

## 6. The self-check modes

Every new kernel can be checked against the code it replaces, bitwise, inside a real run. A mismatch
prints the location and exits non-zero. BUILD.md §4 has the switches.

| mode | build | compares | largest coverage recorded |
|---|---|---|---|
| `MAFFT_SEXTET=check` | normal | every 6-mer pair (stored double) against `commonsextet_p` plus the original conversion | kinds 1–4, every SIMD level, x86 included |
| `MAFFT_MSADIST=check[,forceskip]` | normal | every msa-distance pair against `naivepairscorefast` | 12.5 M pairs on sim5k; bit-plane, exception, column and generic paths; scalar/SSE2/POPCNT/AVX2 |
| `MAFFT_UPG=check` | normal | each guide tree (topol, len, dep, final matrix) against the `_mt` version | full homfam_sdr and Acetyltransf on CARC; stress builds at 1–10 threads |
| `MAFFT_DPCHECK` | timing | every fast DP row against MAFFT's row loop | 11.8 M A__align + 4.3 M G__align11 rows (sim5k) |
| `MAFFT_VECCHECK` | timing | every vectorized cell range against the scalar kernel | 14.7 G cells on 16S.T at 32 threads (CARC) |
| `MAFFT_TILECHECK` | timing | every tiled DP against the serial loop, traceback included | 3,727 forced-tiled DPs per run; 1,476 DPs (3.9 G cells) per run in the gccpar runs |
| `MAFFT_LAZYCHECK` | timing | deferred gap insertion against eagerly updated shadow rows | 1,541 deferred tables on homfam_p450 (above 20,000) |
| `MAFFT_PROFCHECK`, `MAFFT_ARCHECK`, `MAFFT_V8CHECK`, `MAFFT_PSEGCHECK` | timing | profiles, anchor search, parallel Falign pieces, parallel segments | 684 anchor-profile builds and 1,398 in-place segments (16S.T) |

## 7. Compilers and SIMD levels

| build | where | reference | used by |
|---|---|---|---|
| Apple clang 15 arm64, `-O3 -g` | M5 | `baseline-install` | every lever; this package |
| gcc-15.2 arm64, `-O3 -g` | M5 | v7.526 built the same way | gccpar, x86simd, prog3, v7, v8s, v9p |
| gcc-15 arm64, `-O3 -mcpu=native -flto` | M5 | the same flags (= the default build's md5) | gccpar, prog3, v7 |
| clang `-arch x86_64 -O3 -g` under Rosetta, at the default dispatch (SSE2/POPCNT, 2-lane DP), AVX2 (`ROSETTA_ADVERTISE_AVX=1`) and `MAFFT_SIMD=scalar` / `MAFFT_ADPVEC=0` | M5 | `base-x86` | dist, x86simd, prog3, gccpar, v7, v9p |
| clang `-arch x86_64 -march=x86-64-v3 -O3` under Rosetta (FMA class), with and without AVX2 advertised | M5 | `base-x86v3` | x86simd, prog3, gccpar, v7 |
| gcc 13.3 x86-64, default `-O3` and `-O3 -march=native -flto` | CARC | v7.526 gcc 13.3 default | v3, v7 (timing), v8s, v9p, final runs |

Rosetta was used only for correctness. It emulates AVX2 slowly, so its timings mean nothing.

**Thread counts covered:**

- M5: 1, 2, 3, 4, 5, 7, 8, 10, 16, and 32 (oversubscribed).
- CARC: 1, 4, 8 and 32.

## 8. Edge cases

- **The memsave switch (alignments over 30,000 columns).**
  - edge.fa and edge_noshort.fa (59–66 k columns) at 1, 4, 8, 10 and 32 threads, and sim3k_long.
  - At `--thread 1` the output is identical.
  - At more threads it falls within v7.526's own set of outputs. At 10 threads edge_noshort varies
    between two md5s in both programs. At 32 threads edge_noshort gives `6cfea584d675` in both.
  - Above 20,000 the switch phases were forced onto these inputs and matched the serial result.
- **The memsave threshold is alignment width, not sequence length.** 16S.T and 16S.3 contain
  sequences up to 4,066 nt but stay below 30,000 columns (16S.3's alignment has 13,667), so they take
  the normal path.
- **6-mer saturation.** edge2.fa has counts of 255 or more in 5 rows, exercising the overflow
  correction in every SIMD block. It matched at 1, 4 and 10 threads on ARM and at every x86 level.
- **Ambiguity codes.** 16S.T has 59,935 exception columns in 3,443 of 7,350 rows. Sequences shorter
  than 6 kept letters and all-gap or gapless rows were covered by the stress alignments.
- **Above 20,000 sequences.**
  - homfam_p450 (21,013), rnasim25k, 16S.B.ALL (27,643), and on CARC homfam_sdr (50,157),
    homfam_Acetyltransf (46,285) and homfam_rvp (93,681).
  - These run the threaded passes against v7.526's serial passes, whose output does not depend on
    `--thread`. rnasim20k (exactly 20,000) stays on the original threaded path.
- **Repeated runs for races.**
  - rnasim1k at 10 threads, 50 + 60 runs. The race fixed in v8 had shown up once in about 110 runs.
  - The heavy-step budget was stressed with thousands of put-backs.
- **Input handling.**
  - CR and CRLF line ends, Latin-1 and UTF-8 bytes in names, spaces in sequences, stdin input.
  - Three locales (`C`, `C.UTF-8`, `UTF-8`).
  - Names over 254 characters, `=`/`<`/`>`/NUL in sequences, a missing final newline.
  - One fuzz input hangs both the original and fast-mafft, as expected.

## 9. Reference md5s (first 12 hex digits; the same at every thread count unless noted)

For orientation only. The reference in any new environment is v7.526 built there with the same
compiler and flags.

| set | FMA class (Apple clang arm64; clang x86-64-v3) | non-FMA class (gcc; clang x86-64 default) |
|---|---|---|
| sim1k | 3e271a550159 | 85a82397cf39 |
| rnasim579 | 862a63a64dce | 26b86f22a0b1 |
| rnasim1k | d3d860c0602f | 4ec80fa72d83 |
| prot2k | 087f57deb475 | 54420c021453 |
| edge (t1/t4) | 0625e47cebde | cc0dbbde5397 |
| edge_noshort (t1/t4) | 4f0a0d8947bf | 776c3c8d33f6 (t1/t4 on the M5; CARC t8: 4b1480a524ef) |
| edge2 | 96330a6b326d | 8e9d0174d885 |
| 16S.T | 52a159c53210 | f487c1a56702 (M5 gcc-15 and Rosetta x86) |
| sim2k / sim5k / sim10k / sim3k_long | 5ccea931f540 / 2cbec8db3ad2 / fd7d4060b4a7 / 71bca344a3c8 | |
| rnasim5k / rnasim10k | d8f2622e2e80 / 03a55e5dd081 | rnasim10k 5e7bb28c8e8e (gcc-15 arm64) |
| 16S.3 / rnasim20k / homfam_aldosered | 62d15cbda10e / 635823a8a756 / 7fbfa160b9ff | |
| homfam_p450 / 16S.B.ALL / rnasim25k | ebfd9fdd73b0 / b1e342a10946 / 94df3c98061a | p450 a8da25a0a301 (CARC gcc 13.3) |
| homfam_rvp / homfam_sdr / homfam_Acetyltransf | | a6774caadc06 / 4d0ab395b00f / f9d7882de1bf (CARC gcc 13.3) |
