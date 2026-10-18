# Building and running fast-mafft v10

fast-mafft v10 is MAFFT v7.526 (GitLab `sysimm/mafft`, tag `v7.526`, commit `ee97999`) with faster
code in the default strategy, FFT-NS-2. The command line is unchanged. With the same compiler and
flags, stdout is byte-identical to v7.526's. What changed is in [CHANGES.md](CHANGES.md); how
identity was checked is in [VERIFICATION.md](VERIFICATION.md).

## 1. Get the source

Use either of these:

- **`fast-mafft-src/`**: the complete v7.526 source tree with the changes applied. It contains no
  build artifacts.
- **A clean v7.526 tree plus the patch**:

  ```
  git clone https://gitlab.com/sysimm/mafft.git mafft-v10
  git -C mafft-v10 checkout v7.526             # commit ee9799916df6a5d5103d46d54933f8eb6d28e244
  cd mafft-v10 && patch -p1 < /path/to/fast-mafft-v10.patch && cd ..
  ```

  `git apply fast-mafft-v10.patch` also works. It prints harmless trailing-whitespace warnings for 43
  added lines, which keep MAFFT's own style. The patched tree is identical, file for file, to
  `fast-mafft-src/`.

The patch touches 21 files, all under `core/`. Fourteen are changed and seven are new: `adpvec_body.h`,
`crtolf.c`, `dphelp.h`, `gdpvec_body.h`, `lazygap.h`, `mrowblk_body.h` and `progtiming.h`. The new
`crtolf` helper program is added to the Makefile. No other build change is needed.

## 2. Build on Linux with gcc (USC CARC)

```
module load gcc/13.3.0
```

Follow these rules for the baseline and for fast-mafft alike:

1. **Pass the same `PREFIX` to `make` and to `make install`.** The launcher `bin/mafft` is generated
   during `make` from `mafft.tmpl`, and the path to `libexec/mafft` is written into it at that point.
   If you run `make -C core` and then `make -C core install PREFIX=X`, the installed launcher still
   points at `/usr/local/libexec/mafft`. Stock MAFFT behaves the same way. (Checked here: that
   launcher stops with "There is a problem in the configuration of your shell" unless
   `MAFFT_BINARIES` is set.)
2. **Run `make clean` first** in any tree that has been built before. A copied tree can hold a
   generated `core/mafft` from another prefix. `make` does not regenerate it, so `bin/mafft` then
   silently runs the other install's binaries.
3. **Check the launcher after installing:** `grep -m1 '^	prefix=/' $PREFIX/bin/mafft` must print
   `$PREFIX/libexec/mafft`.

### Default build (the identity reference)

This uses the Makefile's defaults: `CC = gcc`, `CFLAGS = -O3`, and `-std=c99` added by the Makefile.

```
build() {  # build SRC PREFIX [CFLAGS]
  make -C "$1/core" clean > /dev/null
  make -C "$1/core" PREFIX="$2" ${3:+CFLAGS="$3"} -j32 > "$2.build.log" 2>&1 &&
  make -C "$1/core" PREFIX="$2" ${3:+CFLAGS="$3"} install >> "$2.build.log" 2>&1
  grep -m1 '^	prefix=/' "$2/bin/mafft"
}
build $PWD/fast-mafft-src $PWD/fm-default
build $PWD/mafft-7.526    $PWD/base-default        # a clean v7.526 tree, for comparison
```

### Aggressive build (the evaluation's "most aggressive compile optimization")

```
build $PWD/fast-mafft-src $PWD/fm-aggr   "-O3 -march=native -flto"
build $PWD/mafft-7.526    $PWD/base-aggr "-O3 -march=native -flto"
```

Notes on flags:

- **With gcc, `-march=native` and `-flto` do not change the output.** MAFFT compiles in ISO C mode
  (`-std=c99`). In that mode gcc defaults to `-ffp-contract=off` and never fuses `a*b + c` into an
  FMA. So `-O3 -march=native -flto` only changes instruction selection, and v7.526 built that way
  gives the same md5 as the default build. This was checked with gcc-15 on the Mac and with gcc 13.3
  on CARC. Use the default build as the reference anyway.
- **With clang, flags can change the output.** Clang contracts `a*b + c` within a single statement
  whenever the target has FMA. A clang `-march=native` or `-march=x86-64-v3` build therefore prints a
  different alignment than a clang default x86-64 build. This is true for v7.526 and fast-mafft
  alike. Always compare against a baseline built with the same compiler and flags.
- **Do not use `-Ofast`, `-ffast-math` or `-ffp-contract=fast`.** They change v7.526's own output.
- **`-mavx2` is not needed.** The SIMD kernels choose SSE2, POPCNT or AVX2 at run time with
  `__builtin_cpu_supports`. Under `-march=native` the same code is compiled for AVX2 directly.
- **`-g` does not change the generated code.**

### glibc notes

MAFFT compiles with `-std=c99` (strict ISO). Under that flag glibc hides POSIX and GNU declarations
unless a feature macro is defined. The baseline relies on the bare `#define _XOPEN_SOURCE` in
`mltaln.h`. fast-mafft needs a few more declarations. This is how the final tree handles them
(checked with `grep -n _GNU_SOURCE core/*.c`):

- **`disttbfast.c` and `Salignmm.c` begin with `#define _GNU_SOURCE`**, guarded by `#ifndef`, before
  any `#include`.
  - `disttbfast.c` needs it for `pthread_rwlock_t` and for the anonymous `mmap` of the distance
    matrix (`MAP_ANONYMOUS`, `madvise(MADV_HUGEPAGE)`).
  - `Salignmm.c` needs it for `posix_memalign` and `nanosleep`.

  Without the macro, gcc 13 warns "implicit declaration of function 'posix_memalign'", and gcc 14 or
  later rejects the file. v7 had this bug; v8 fixed it. Nothing needs to be added to `CFLAGS`.
- **The other changed files need no macro.** They call only functions that glibc declares
  unconditionally or under MAFFT's own `_XOPEN_SOURCE`:
  - `sched_yield`, `sysconf` and `gettimeofday` in `mltaln9.c`;
  - `fstat`, `fileno`, `isatty` and `setvbuf` in `f2cl.c`.

  `io.c` declares `flockfile`, `getc_unlocked` and `funlockfile` itself.
- **Instrumented timing builds need `-D_GNU_SOURCE`** in `CFLAGS` (`-DPROGTIMING`; see section 4),
  for `clock_gettime` and `strdup`. Normal builds do not.
- **Tested so far:** the final source was built on CARC with gcc 13.3, with default flags and with
  `-O3 -march=native -flto` (job 13066113). Locally it builds with Apple clang 15 and with gcc-15
  (arm64). gcc 14 or later on glibc has not been tried.
- **Transparent huge pages.** fast-mafft holds the N(N+1)/2 distance half-matrix in one anonymous mapping.
  When `/sys/kernel/mm/transparent_hugepage/enabled` is `madvise`, it requests huge pages with
  `madvise(MADV_HUGEPAGE)`. Under `always` (CARC's setting) it gets them without asking.
  `MAFFT_V9SLAB=0` restores MAFFT's per-row `calloc`. The output does not change either way.
- **Memory.** If peak RSS is a concern, measure it with `/usr/bin/time -v` (maximum resident set
  size). `MALLOC_ARENA_MAX=4` limits glibc's per-thread arenas and does not change the output.

## 3. Run

```
fm-default/bin/mafft --thread 32 in.fa > out.fa
```

The command line is MAFFT v7.526's. `mafft in.fa` with no strategy option runs FFT-NS-2, and the
alignment on stdout is byte-identical to v7.526's at the same `--thread` value.

- **`--thread N`.** The 6-mer distances, the guide trees and the msa distances use N threads.
  - v7.526 caps the progressive passes at 16 threads when N ≥ 16. The launcher passes
    `disttbfast -C 32-16`.
  - fast-mafft passes `-C N-16-N`. It uses all N threads wherever the progressive result cannot depend on
    the order of the steps or on the thread count: the ready queue, and every pass above 20,000
    sequences.
  - It keeps v7.526's 16 step-claiming threads in the one order-dependent case: MAFFT's in-order
    dispatch, used for a pass in which the memsave switch fires (an alignment grows past 30,000
    columns, with at most 20,000 sequences). CHANGES.md §3.6 explains it. v10 starts the other
    N − 16 threads there too, but they only help inside steps (CHANGES.md §8.3).
- **Above 20,000 sequences.** v7.526 runs both progressive passes on one core whatever `--thread`
  says. fast-mafft runs them threaded with the same result. `--threadtb 1` restores the serial passes.
- **Compare runs at the same `--thread` value.** `mafft` without `--thread` takes MAFFT's
  unthreaded code path, and that can give a different alignment than `--thread N` (v7.526 does the
  same).
- **stderr.** Progress messages and stage markers keep v7.526's wording. From 4,000 sequences on,
  the guide-tree progress counter is printed every 100 merges (every 1000 from 40,000 sequences)
  instead of every 10, because each write cost about 120 µs when stderr is a pipe. Identity is about
  stdout only.
- **Other strategies.** `--auto`, L-INS-i, FFT-NS-i, `--add` and so on share some of the changed
  code: input reading, `f2cl`, the launcher, and the A__align/G__align11 kernels. During development
  they were spot-checked for identical output, mostly on earlier builds (VERIFICATION.md §6). They
  were not tuned or timed. In v7.526 itself the iterative modes are not deterministic at
  `--thread` > 1.

Identity check: md5 of stdout, same thread count, both binaries built with the same compiler and
flags.

```
h() { "$1" --thread "$2" "$3" 2>/dev/null | md5sum | cut -c1-12; }
for d in sim1k rnasim1k prot2k 16S.T; do for t in 1 8 32; do
  echo "$d t=$t base=$(h base-default/bin/mafft $t $d.fa) fm=$(h fm-default/bin/mafft $t $d.fa)"
done; done
```

## 4. Self-check and diagnostic switches (environment variables)

These were collected by grepping every `getenv` in the final `core/`. A failed check prints the
mismatch and exits non-zero. The self-checks keep the output identical but make runs slower, some
by a lot.

### Normal build (no rebuild needed)

| variable | effect |
|---|---|
| `MAFFT_SEXTET=check` | Every 6-mer pair is also computed with v7.526's `makecompositiontable_p` + `commonsextet_p` plus the original conversion, and the stored double is compared bitwise. Prints `sextetfast check passed: kind K, simd level S, ...`. |
| `MAFFT_SEXTET=...` (other words, combinable with `,`) | `ref` uses the original code for every pair. `time` prints phase times, the kernel kind, the SIMD level and block counts. `scalar` skips the SIMD block path. `list` (DNA) uses the list path. `inv` / `noinv` force or forbid the protein inverted index (kind 3). `dense` / `nodense` force or forbid the protein element masks (kind 4). `counts` stores counts and lets the original conversion loop run. `nolftab` computes `lenfac` for every pair. |
| `MAFFT_SEXTET_DUMP=/abs/path` | Write the finished 6-mer distance matrix. |
| `MAFFT_MSADIST=check` | Every pair of the "distance matrix from msa" stage is also scored by v7.526's `naivepairscorefast` and compared bitwise. Silent on success; prints `msadist check FAILED at i=.. j=..` and exits on a mismatch. |
| `MAFFT_MSADIST=...` (combinable) | `stat` prints the path taken: bit-plane or column, avx2/popcnt/neon/portable, exception kinds, dense columns. `ref` uses the original scorer. `time` prints thread-seconds per part. `generic`, `ctz`, `noexc`, `nolane`, `nodense` select older or simpler paths. `forceskip` / `noskip` force or forbid skipping all-gap words in the gap-run loops. `nont` turns off the x86 non-temporal matrix stores. |
| `MAFFT_MSADIST_DUMP=/abs/path`, `MAFFT_MSADIST_SAVEMSA=/abs/path` | Write the finished matrix, or the aligned rows plus `amino_dis`, `penalty_dist` and self-scores (input for the standalone kernel benchmark). |
| `MAFFT_UPG=check` | Build each guide tree with the v9 routine and, from a copy of the matrix, with the earlier `_mt` routine. That routine was checked against v7.526's own function. `topol`, `len`, `dep` and the matrix left behind are compared bitwise. Prints `UPGMA check passed: ...`. Needs memory for a second matrix. |
| `MAFFT_UPG=old` | Use the `_mt` guide-tree routine. |
| `MAFFT_UPGTHREADS=n` | At most n guide-tree threads. The tree does not depend on n. |
| `MAFFT_UPGSTAT=1` | Guide-tree timers and counts on stderr. |
| `MAFFT_SIMD=scalar\|sse2\|popcnt\|avx2` (`neon` on ARM) | Lower the SIMD level of the 6-mer and msa-distance kernels. A level above the CPU's prints a warning and uses the CPU's level. Use it with the two `check` modes to test each level. |
| `MAFFT_MSCHECK=1` | (v10) Every Myers-Miller call (`MSalignmm`, past 30,000 columns) is also run with the original recursion; the score, both gap tables, every row and the fused profiles are compared. Prints `MSCHECK: N MSalignmm calls ... no mismatch`. About MAFFT's own time on long genomes. |
| `MAFFT_MSKCHECK=1` | (v10) Every recursion node's forward and backward sweeps are recomputed by the original loops; the midpoint, jumps and split are compared. Slow. |
| `MAFFT_ARLCHECK=1` | (v10) The anchor search: every lag's segments against the original `alignableReagion`, the lag lists against the original `getKouho`, the FFT input vectors against the original loops, and every lag scored by a helper. |
| `MAFFT_BA2CHECK=1` | (v10) Every `blockAlign2` call also runs the full original loops; the cut points are compared. |
| `MAFFT_FFTCHECK=1` | (v10) Every transform of 16,384 points or more is also computed in the original order and compared bitwise. |
| `MAFFT_MSFAST=0`, `MAFFT_UDLV10=0` | (v10, A/B) The original `MSalignmm` recursion; the original `Falign_udpari_long`. |
| `MAFFT_PSEG=0`, `MAFFT_LAG2=0` | (v10, A/B) Align Falign's segments one after another; score the FFT lags one at a time. |
| `MAFFT_MSFORCEPAR=1`, `MAFFT_MSBANDW=n` | (v10, tests) Use the Myers-Miller task stack and bands even without idle helpers; band width in columns (default 1,024). |
| `MAFFT_ADPVEC=0\|2\|4` | Width of the progressive DP row kernels. 4 = AVX2 lanes, 2 = SSE2/NEON lanes, 0 = scalar. Default: 4 on x86 with AVX2, 2 on other x86, 0 on ARM. In v10 it also sets the Myers-Miller cell kernels. |
| `MAFFT_V9SLAB=0\|1\|2` | Distance half-matrix: `0` uses v7.526's per-row `calloc`, `1` uses one mapping and never calls madvise, `2` always calls madvise(MADV_HUGEPAGE). By default it calls madvise only when the THP mode is `madvise`. |
| `MAFFT_V9PTIME=1` | Stage timers for the distance stages, and huge-page use. |
| `MAFFT_BINARIES=dir` | v7.526's own variable: run the programs in `dir` instead of the installed `libexec/mafft`. |

Recommended self-check run (it passes on every set tried; see VERIFICATION.md):

```
MAFFT_SEXTET=check MAFFT_MSADIST=check MAFFT_UPG=check fm-default/bin/mafft --thread 32 in.fa | md5sum
MAFFT_ARLCHECK=1 MAFFT_BA2CHECK=1 MAFFT_FFTCHECK=1 MAFFT_MSCHECK=1 fm-default/bin/mafft --thread 32 long.fa | md5sum
for s in scalar sse2 avx2; do MAFFT_SIMD=$s MAFFT_SEXTET=check MAFFT_MSADIST=check,forceskip fm-default/bin/mafft --thread 8 in.fa | md5sum; done
```

### Timing build (`-DPROGTIMING`): progressive-pass self-checks, logs and A/B switches

Only `disttbfast` needs rebuilding. Point an installed launcher at it with `MAFFT_BINARIES`:

```
cp -R fast-mafft-src fm-pt && make -C fm-pt/core clean > /dev/null
make -C fm-pt/core CFLAGS="-O3 -DPROGTIMING -DGCCPAR_TBTEST -D_GNU_SOURCE" disttbfast -j32
mkdir -p fm-pt/libexec && cp fm-pt/core/disttbfast fm-pt/libexec/
for f in fm-default/libexec/mafft/*; do [ "$(basename $f)" = disttbfast ] || ln -sf $PWD/$f fm-pt/libexec/; done
MAFFT_BINARIES=$PWD/fm-pt/libexec MAFFT_VECCHECK=1 MAFFT_DPCHECK=1 MAFFT_LAZYCHECK=1 \
  fm-default/bin/mafft --thread 8 rnasim1k.fa | md5sum
```

| variable (timing build) | checks or does |
|---|---|
| `MAFFT_VECCHECK=1` | Every vectorized DP cell range (A__align 4-letter and general-letter kernels, G__align11) is redone by the scalar kernel on shadow copies, and m, mp, the row, the traceback cells, mi, mpi and wm are compared bitwise. Combine with `MAFFT_ADPVEC=2` or `4` to check a given width. |
| `MAFFT_DPCHECK=1` | Every fast DP row is compared with v7.526's original row loop. |
| `MAFFT_TILECHECK=1` | Every tiled (multi-thread) DP is redone serially from the saved state and compared, traceback codes included. To force tiling: `MAFFT_TILEFORCE=1 MAFFT_TILEMINCELLS=0 MAFFT_TILEWIDTH=16`. |
| `MAFFT_LAZYCHECK=1` | Deferred gap insertion: shadow rows receive every gap table immediately and must equal the real rows whenever those are current. |
| `MAFFT_PROFCHECK=1`, `MAFFT_ARCHECK=1` | One-pass profiles and the dense first row and column; the hoisted anchor search. Each is compared with the original functions. |
| `MAFFT_V8CHECK=1` | Parallel FFT input vectors and transform pairs, anchor-profile builds, parallel row copies and in-place Falign segments, each redone by the original code and compared. |
| `MAFFT_PSEGCHECK=1` (`MAFFT_PSEG=1`, `MAFFT_PSEGFORCE=1`) | Parallel Falign segments are also run sequentially and compared. (In v10 they are on by default for calls with at least 16 segments; `MAFFT_PSEG=1` also uses them for smaller calls, `MAFFT_PSEGFORCE=1` for every call.) |
| `MAFFT_PROGLOG=/abs/file`, `MAFFT_V8STAT=1` | Per-step log of the progressive passes; row-kernel and tiling statistics. |
| A/B switches | `MAFFT_NOFASTDP` (original DP loops), `MAFFT_INORDER` (original in-order dispatch), `MAFFT_READYQ1` (ready queue at one thread), `MAFFT_NOHELP`, `MAFFT_NOHELPFOR`, `MAFFT_NOTILE`, `MAFFT_NOGENERAL`, `MAFFT_NOLAZY`, `MAFFT_LZTWOPHASE=0`, `MAFFT_COPYROWS` (row copies under `aseqlock`; this restores the old deadlock-prone path), `MAFFT_SCHEDPRIO=0` (old priorities), `MAFFT_CRITHELP=1`, `MAFFT_FAINPLACE=0`, `MAFFT_LARGEN_TB1` |
| Tiling parameters | `MAFFT_TILEWIDTH`, `MAFFT_TILEMINCELLS`, `MAFFT_TILEFORCE`, `MAFFT_TILEOLD`, `MAFFT_TILENBFAC`, `MAFFT_TILEALIGN`, `MAFFT_TILEBATCH`, `MAFFT_TILECLAIM`, `MAFFT_TILEBANDMAJOR`, `MAFFT_TILELOCALSETUP`, `MAFFT_CRITNB`, `MAFFT_CRITMINCELLS`. `MAFFT_TILEIJPTEST` is for diagnosis only: it makes the output wrong on purpose. |

### Test build (`-DGCCPAR_TBTEST`)

| variable | effect |
|---|---|
| `MAFFT_SERIALTB=1` | Use the serial `treebase()` progressive path, as v7.526 does above 20,000 sequences. |
| `MAFFT_LARGEN_ORIGSERIAL=1` | Apply the above-20,000 logic at any size, for example to test the memsave-switch phases on `edge.fa`. |
| `MAFFT_LARGEN_HEAVY=f`, `MAFFT_LARGEN_HEAVYMIN=bytes`, `MAFFT_LARGEN_KEEP=bytes` | Memory budget for concurrent large DPs, and arena trimming. `HEAVYMIN=1 HEAVY=0.000001` stresses the budget. |
