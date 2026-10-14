# fast-mafft v9

fast-mafft is a faster drop-in replacement for **MAFFT v7.526** in its default mode (FFT-NS-2,
`mafft --thread N in.fa`). For the same compiler and flags, its stdout is **byte-identical** to
v7.526's. It was built for the USC CSCI 270 bonus project and is evaluated on USC CARC: 32 cores,
64 GB, large datasets, against v7.526 built with aggressive optimization.

## Contents

| path | what it is |
|---|---|
| `fast-mafft-src/` | Complete MAFFT v7.526 source tree with the changes, and no build artifacts. Build it with the stock Makefile (see BUILD.md). |
| `fast-mafft-v9.patch` | Unified diff against pristine v7.526 (GitLab `sysimm/mafft`, tag v7.526, commit `ee97999`). It covers 19 files under `core/`: 12 changed, 7 new. Apply with `patch -p1` or `git apply` from the top of a clean v7.526 tree. |
| [BUILD.md](BUILD.md) | Build on Linux/gcc (CARC: `module load gcc/13.3.0`) with default and aggressive flags, glibc notes, how to run, and every environment switch. |
| [CHANGES.md](CHANGES.md) | Every optimization, stage by stage: what MAFFT did, what v9 does, why the output stays identical, and the measured effect. Also the bugs found and fixed. |
| [VERIFICATION.md](VERIFICATION.md) | How identity was verified: datasets, thread counts, compilers and SIMD levels, self-check modes, edge cases, reference md5s. |
| [RESULTS.md](RESULTS.md) | Final 3-run CARC evaluation on 12 datasets (`--thread 32`, 32 cores, 64 GB): **9.7x–36.2x mean speedup, minimum single run 8.3x, byte-identical output in every run**. Also methodology, datasets and raw runs. |
| `carc/` | The CARC job script and timing harness (copies of `bench/carc_v9_final.sh` and `bench/carc_run.py`), a log summarizer, the three final job logs (`carc/logs/`), and a README on re-running the evaluation. |

## Quick start (Linux)

```
module load gcc/13.3.0
make -C fast-mafft-src/core clean
make -C fast-mafft-src/core PREFIX=$HOME/fast-mafft -j32
make -C fast-mafft-src/core PREFIX=$HOME/fast-mafft install     # same PREFIX on both commands
$HOME/fast-mafft/bin/mafft --thread 32 in.fa > out.fa
```

## In one paragraph

v9 changes every stage of FFT-NS-2:

- **6-mer distances:** a byte-SAD identity (|a − b| summed with `psadbw`/`uaba`), inverted indexes
  and bit masks.
- **Guide trees:** the same greedy clustering, with candidate lists that certify most
  nearest-neighbour rescans.
- **Progressive passes:**
  - branch-free, fused and AVX2-vectorized Gotoh DP rows, with the row carry computed by a
    prefix scan;
  - critical-path list scheduling;
  - idle threads tiling large DPs;
  - in-place rows and deferred gap insertion;
  - threaded passes above 20,000 sequences, where MAFFT is serial.
- **Distances from the MSA:** bit-parallel popcount kernels.
- **I/O:** much less overhead.

Every floating-point value is produced by the same operations in the same order, integer sums are
exact in any order, tie rules are kept, and the one order-dependent piece of MAFFT (the memsave
switch) is reproduced explicitly. Each new kernel has a bitwise self-check against the original
code.

On CARC development runs, v9 measured 10–11x on rnasim10k and 16S.T, and 15–50x on the large HomFam
sets, against the aggressive baseline. The final numbers will go in RESULTS.md.

## Verification of this package (2026-10-13)

- **The patch applies cleanly.** It applies to a clean `git archive` of v7.526 with `patch -p1` (no
  fuzz or offsets) and with `git apply`. Both results are identical to `fast-mafft-src/` in content
  and permission bits.
- **It builds.** Apple clang 15 builds the patched tree.
- **The output is identical.** md5 at `--thread 4` matches `baseline-install` on sim1k
  (3e271a550159), rnasim1k (d3d860c0602f) and prot2k (087f57deb475).
- **The self-checks pass.** `MAFFT_SEXTET=check MAFFT_MSADIST=check MAFFT_UPG=check` passes on
  rnasim1k and prot2k.

Details are in VERIFICATION.md §3.
