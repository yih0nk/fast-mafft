# fast-mafft v10

fast-mafft is a faster drop-in replacement for **MAFFT v7.526** in its default mode (FFT-NS-2,
`mafft --thread N in.fa`). For the same compiler and flags, its stdout is **byte-identical** to
v7.526's (one exception at more than one thread, where MAFFT's own output depends on thread timing:
see "Known limits"). It was built for the USC CSCI 270 bonus project and is evaluated on USC CARC: 32 cores,
64 GB, large datasets, against v7.526 built with aggressive optimization.

## How to evaluate (USC CARC)

```
git clone https://github.com/yih0nk/fast-mafft.git && cd fast-mafft
bash bench/get_data.sh                     # optional: our 32 benchmark inputs (93 MB), checksum-verified
bash bench/get_data.sh --long              # optional: + 12 long-genome inputs, e.g. real SARS-CoV-2 (50 MB)
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
- when the two outputs of a run differ, both programs run again, and the summary says whether MAFFT
  reproduced its own output (see "Known limits" below);
- `--threads`, `--reps`, `--flags` and `--work` are listed by `bash carc/run_eval.sh --help`.

Run `bench/get_data.sh` on a login node, since compute nodes may have no internet access.

## Contents

| path | what it is |
|---|---|
| `fast-mafft-src/` | Complete MAFFT v7.526 source tree with the changes, and no build artifacts. Build it with the stock Makefile (see BUILD.md). |
| `fast-mafft-v10.patch` | Unified diff against pristine v7.526 (GitLab `sysimm/mafft`, tag v7.526, commit `ee97999`). It covers 21 files under `core/`: 14 changed, 7 new. Apply with `patch -p1` or `git apply` from the top of a clean v7.526 tree. |
| [BUILD.md](BUILD.md) | Build on Linux/gcc (CARC: `module load gcc/13.3.0`) with default and aggressive flags, glibc notes, how to run, and every environment switch. |
| [CHANGES.md](CHANGES.md) | Every optimization, stage by stage: what MAFFT did, what fast-mafft does, why the output stays identical, and the measured effect (§8: the v10 changes for long genomes). Also the bugs found and fixed. |
| [VERIFICATION.md](VERIFICATION.md) | How identity was verified: datasets, thread counts, compilers and SIMD levels, self-check modes, edge cases, reference md5s. |
| [RESULTS.md](RESULTS.md) | Final 3-run CARC evaluation (`--thread 32`, 32 cores, 64 GB): on 12 large benchmark datasets **10.1x–43.9x mean speedup, slowest single run 9.98x, byte-identical in every run**; real SARS-CoV-2 collections 11.6x–205x; long genomes 14.9x–72.8x; smaller inputs. Also methodology, datasets, raw runs, and the v9 evaluation. |
| `carc/` | `run_eval.sh`, the portable evaluation job above; the scripts, summarizer and logs of the final evaluations (`carc_v10_final.sh`, `carc_v9_final.sh`, `summarize.py`, `logs/`); and a README. |
| `bench/` | `get_data.sh` and `DATA.md` (benchmark inputs, provenance, licenses, checksums in `DATA_MD5.txt` and `DATA_MD5_LONG.txt`), the timing harness `carc_run.py`, the development harness `run_bench.py`, and the synthetic-data generators. |

## Quick start (Linux)

```
module load gcc/13.3.0
make -C fast-mafft-src/core clean
make -C fast-mafft-src/core PREFIX=$HOME/fast-mafft -j32
make -C fast-mafft-src/core PREFIX=$HOME/fast-mafft install     # same PREFIX on both commands
$HOME/fast-mafft/bin/mafft --thread 32 in.fa > out.fa
```

## In one paragraph

v9, the first complete version, changes every stage of FFT-NS-2:

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

v10, the current version, adds the paths v9 left at MAFFT's speed:

- **Alignments wider than 30,000 columns** (long genomes), where MAFFT switches to its memsave
  aligners: a branch-free, AVX2, parallel Myers-Miller; an anchor search that builds its profiles
  once instead of 500 times per merge and scores the lags in parallel.
- **Collections of near-identical genomes** (e.g. SARS-CoV-2): parallel Falign segments, an O(n²)
  `blockAlign2`, cache-blocked long FFTs, and a shared pool for the large per-thread buffers.

Every floating-point value is produced by the same operations in the same order, integer sums are
exact in any order, tie rules are kept, and MAFFT's one order-dependent piece (the memsave switch) is
reproduced explicitly. Each new kernel has a bitwise self-check against the original code.

In the final CARC evaluation (3 runs, `--thread 32`), v10 was 10.1x–43.9x faster on average than MAFFT
v7.526 built with `-O3 -march=native -flto` on 12 large benchmark datasets, with a slowest
single run of 9.98x and byte-identical output in every run; 11.6x–205x on real SARS-CoV-2
collections; 14.9x–72.8x on long genomes. See RESULTS.md.

## Known limits

- **Alignments that cross 30,000 columns in the middle of a progressive pass**, at more than one
  thread (in our tests: 2,000 SARS-CoV-2 genomes, and simulated virus collections). There MAFFT
  switches modes while other steps are in flight, and its own output depends on thread timing; it
  usually gives one output on a given machine, but not always. fast-mafft keeps the same mechanism,
  but its timing is different, so it does not always give MAFFT's usual output (11 of 16 runs on 2,000
  SARS-CoV-2 genomes at 32 threads). At `--thread 1` it is identical. Inputs whose alignment stays
  under 30,000 columns, or crosses it in the first steps (most long-genome sets), were identical in
  every run. Details: VERIFICATION.md §3, CHANGES.md §8.4.
- **Small inputs:** about 5x at 1,000 sequences and 7x at 2,000 (MAFFT takes 6–20 s there);
  10x or more on average from 5,000 sequences.

## Verification

See VERIFICATION.md: identity across compilers, SIMD levels, thread counts, BAliBASE and edge
cases, plus the self-check modes.

## License

BSD, as MAFFT's own code: see [LICENSE](LICENSE). MAFFT is by Kazutaka Katoh; his notice is kept in
`fast-mafft-src/license`, and the unchanged `extensions/` keep their own licenses
(`fast-mafft-src/license.extensions`).
