# fast-mafft v9

fast-mafft is a faster drop-in replacement for **MAFFT v7.526** in its default mode (FFT-NS-2,
`mafft --thread N in.fa`). For the same compiler and flags, its stdout is **byte-identical** to
v7.526's. It was built for the USC CSCI 270 bonus project and is evaluated on USC CARC: 32 cores,
64 GB, large datasets, against v7.526 built with aggressive optimization.

## How to evaluate (USC CARC)

```
git clone https://github.com/yih0nk/fast-mafft.git && cd fast-mafft
bash bench/get_data.sh                     # optional: our 32 benchmark inputs (93 MB), checksum-verified
sbatch --account=<account> carc/run_eval.sh                               # all of bench/data, --thread 32
sbatch --account=<account> carc/run_eval.sh --data /path/to/fastas --reps 3   # your own inputs
```

`carc/run_eval.sh` builds **both** MAFFT v7.526 (from this repository's first commit, checked against
known checksums) and fast-mafft on the node that runs the job, with the same compiler (gcc 13.3) and
the same flags (default `-O3 -march=native -flto`). It then runs the two programs back to back on
every input at `--thread 32`, checks that their outputs are **byte-identical**, and prints a table of
wall time, speedup and peak memory to `fast-mafft-eval-<jobid>.log`. Options:

- `--check-default` also compares the outputs of builds with the default Makefile flags;
- `--selfcheck` also runs fast-mafft's internal bitwise self-checks;
- `--threads`, `--reps`, `--flags` and `--work` are listed by `bash carc/run_eval.sh --help`.

Run `bench/get_data.sh` on a login node, since compute nodes may have no internet access.

## Contents

| path | what it is |
|---|---|
| `fast-mafft-src/` | Complete MAFFT v7.526 source tree with the changes, and no build artifacts. Build it with the stock Makefile (see BUILD.md). |
| `fast-mafft-v9.patch` | Unified diff against pristine v7.526 (GitLab `sysimm/mafft`, tag v7.526, commit `ee97999`). It covers 19 files under `core/`: 12 changed, 7 new. Apply with `patch -p1` or `git apply` from the top of a clean v7.526 tree. |
| [BUILD.md](BUILD.md) | Build on Linux/gcc (CARC: `module load gcc/13.3.0`) with default and aggressive flags, glibc notes, how to run, and every environment switch. |
| [CHANGES.md](CHANGES.md) | Every optimization, stage by stage: what MAFFT did, what v9 does, why the output stays identical, and the measured effect. Also the bugs found and fixed. |
| [VERIFICATION.md](VERIFICATION.md) | How identity was verified: datasets, thread counts, compilers and SIMD levels, self-check modes, edge cases, reference md5s. |
| [RESULTS.md](RESULTS.md) | Final 3-run CARC evaluation on 12 datasets (`--thread 32`, 32 cores, 64 GB): **9.7x–36.2x mean speedup, minimum single run 8.3x, byte-identical output in every run**. Also methodology, datasets and raw runs. |
| `carc/` | `run_eval.sh`, the portable evaluation job above; the scripts, summarizer and logs of the final evaluation (`carc_v9_final.sh`, `summarize.py`, `logs/`); and a README. |
| `bench/` | `get_data.sh` and `DATA.md` (benchmark inputs, provenance, licenses, checksums in `DATA_MD5.txt`), the timing harness `carc_run.py`, the development harness `run_bench.py`, and the synthetic-data generators. |

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

In the final CARC evaluation (3 runs on 3 nodes, 12 datasets, `--thread 32`), v9 was 9.7x–36.2x faster
on average than MAFFT v7.526 built with `-O3 -march=native -flto`. Its slowest single run was 8.3x,
and the output was byte-identical in every run. See RESULTS.md.

## Verification

See VERIFICATION.md: identity across compilers, SIMD levels, thread counts, BAliBASE and edge
cases, plus the self-check modes.
