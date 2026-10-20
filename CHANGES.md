# fast-mafft v10: what changed, and why the output is still byte-identical

This document explains every optimization in fast-mafft, stage by stage. §1–§7 are v9, which covers
every stage of MAFFT's default mode; §8 is v10, which adds alignments wider than 30,000 columns (long
genomes) and collections of near-identical genomes. For each change it covers what MAFFT v7.526 does,
what fast-mafft does instead, why the alignment printed on stdout cannot change, and
what the change measured. All numbers are copied from the project's development notes (one notes
file per lever, kept privately and not published in this repository) or from the CARC logs in
`carc/logs/`, with the machine stated: "M5" is the development Mac (Apple M5, 4 performance + 6 efficiency cores,
clang); "CARC" is a USC CARC node (AMD EPYC 7513, gcc 13.3, `--thread 32`). The final CARC
evaluation goes in [RESULTS.md](RESULTS.md).

## 0. Background

### 0.1 What MAFFT's default mode does

`mafft in.fa` (FFT-NS-2, with `--thread N`) runs one program, `disttbfast`, which does this:

1. **6-mer distances.** Count shared 6-mers for all N(N−1)/2 pairs and convert the counts to
   distances in a double half-matrix.
2. **Guide tree 1.** A UPGMA-like agglomerative clustering of that matrix.
3. **Progressive alignment 1/2.** One profile-profile alignment per tree node, leaves to root, using
   FFT anchors and a Gotoh-style dynamic program (DP).
4. **Distance matrix from msa.** Re-score all pairs from the pass-1 alignment.
5. **Guide tree 2**, built from those scores.
6. **Progressive alignment 2/2** on the new tree.

The shell launcher then reformats the result with `f2cl`.

### 0.2 Where the baseline spends its time

On the M5 (sim5k, 10 threads, 51.3 s), a profile of the baseline splits the run as follows:

| stage | share | hot function |
|---|---:|---|
| progressive passes | 46% | `A__align` (41.6% of all wall time) |
| distance from msa | 39% | `naivepairscorefast` (36.7%) |
| 6-mer distances | 10% | `commonsextet_p` (10.4%) |
| serial work (script, parsing, `f2cl`, trees) | 5.5% | |

On CARC (rnasim10k, 32 threads, 99.3 s) the split was:

- progressive passes: 37.8 + 25.9 s;
- distance from msa: 22.2 s;
- 6-mer distances: 11.1 s;
- everything else: about 2.3 s.

Above 20,000 sequences MAFFT runs both progressive passes on one core. On the largest protein sets
the two guide trees also become dominant (§2, §6).

### 0.3 How v9 was built up

v9 was built in rounds. Separate "levers" worked in parallel copies of the source and were then
merged:

| lever | work |
|---|---|
| dist | distance from msa |
| prog, prog2, prog3 | progressive passes |
| sextet | 6-mer distances |
| overhead | I/O, script, trees |
| gccpar | gcc code generation |
| x86simd | x86 kernels |
| largen | more than 20,000 sequences |
| v8s | IUPAC rRNA at 32 threads |
| v9p | huge protein N |

Each lever's development notes record its measurements (not published here). v3 was the first build tested on CARC, v7 the
merge of the round-3 levers, v8 = v7 + v8s, and **v9 = v8 + v9p**.

### 0.4 The exactness rules every change follows

"Byte-identical" means the same md5 of stdout as v7.526 built with the same compiler and flags, at
the same `--thread`. Every change keeps to six rules:

1. **Integer quantities may be computed in any order.**
   - Shared 6-mer counts, the msa-distance score and self-scores are sums of integers.
   - MAFFT sometimes accumulates them in `double`, but every partial sum is an integer far below
     2^53, so the double result equals the exact integer sum.
   - Such sums can be reordered, vectorized, computed with popcounts or split across threads freely.
2. **Every floating-point value is produced by the same operations on the same operands in the same
   order.**
   - Sums over group members or alphabet letters keep their original order.
   - Each `a*b + c` stays one C statement. Whether it becomes an FMA is then decided by the compiler
     and flags exactly as in the baseline:
     - clang contracts within a statement when the target has FMA;
     - gcc in MAFFT's `-std=c99` mode never contracts (`-ffp-contract=off`).
   - New SIMD code is compiled with `target("avx2")`, never `"fma"`.
   - The x86 kernels of the distance stages carry a target attribute only on integer code. All their
     double arithmetic stays in default-target code.
3. **Padding terms must be exact no-ops.** A kernel may add terms MAFFT skipped only when they are
   exact no-ops:
   - `acc + s*0.0` and `acc + 0.0*w` return `acc` bit for bit whenever `acc` is not −0.0;
   - MAFFT's match accumulators start at +0.0, add non-negative weights, and can never become −0.0.
4. **Tie rules are kept.**
   - The DP keeps its exact mix of `>` and `>=`.
   - Minimum searches keep "strict `<`, first index wins".
   - Where a sequential selection is parallelized, the parallel combination is proven to return the
     same element on ties (§2, §3.1).
5. **Schedule independence.**
   - Each progressive step is a pure function of its two children's rows, so any topological order
     of the tree gives the same alignment.
   - The one piece of order-dependent state in MAFFT, the "memsave" switch at 30,000 columns, is
     reproduced explicitly (§3.6, §6).
6. **Every new kernel has a bitwise self-check.** The original code recomputes the same value and the
   two are compared, either in the normal build (`MAFFT_SEXTET=check`, `MAFFT_MSADIST=check`,
   `MAFFT_UPG=check`) or in an instrumented build (`MAFFT_DPCHECK`, `VECCHECK`, `TILECHECK`,
   `LAZYCHECK`, `PROFCHECK`, `ARCHECK`, `V8CHECK`). See [VERIFICATION.md](VERIFICATION.md).

---

## 1. The 6-mer distance matrix

### What MAFFT does

- **Encoding.** Each sequence becomes a list of rolling 6-mer codes. DNA has 4^6 = 4096 codes;
  ambiguity letters are dropped. Protein has 6^6 = 46656 codes over 6 residue groups.
- **Pair loop.** For each row i, `distancematrixthread` allocates a 4096-int count table, then calls
  `commonsextet_p(table_i, pointt[j])` for every j ≥ i. That returns the integer Σ_p min(c_i(p),
  c_j(p)) by walking j's list with a thread-local memo table: random accesses, latency-bound, about
  **2.2 µs per pair** (sim5k).
- **Conversion.** A serial O(N²) pass turns the counts into distances:
  `(1 − cnt/min(self_i, self_j)) · lenfac · 2`, where `lenfac` depends on the two lengths.

### What v9 does

- **DNA, "kind 1": dense saturated byte tables and a sum-of-absolute-differences kernel.**
  - Each row gets a 4096-byte table t_i(p) = min(c_i(p), 255).
  - The pair count uses the identity Σ min(a, b) = (Σa + Σb − Σ|a − b|) / 2. A *min*-reduction is
    thereby turned into |a − b|, which the hardware sums directly:
    - NEON `uaba` on ARM;
    - SSE2 `psadbw` / AVX2 `vpsadbw` on x86, which sum 8 absolute byte differences per 64-bit lane.
  - The kernel processes 4×4 blocks of rows × columns. Threads claim "super-rows" of 16 rows, so each
    16 KB block of column tables is reused by four row blocks.
  - **Saturation correction.** The byte cap only matters when both counts are ≥ 255. Those tuples are
    kept in sorted overflow lists, and their correction min(c_i, c_j) − 255 is added by merging the
    two lists.
- **Protein: three exact kernels, chosen per input by an operation-count estimate.**
  - *Kind 2, list path.* A dense u16 table of row i, summed over j's list of distinct tuples.
  - *Kind 3, inverted index.* For every tuple p, the ascending list of rows holding it. Row i adds
    min(c_i(p), c_j(p)) into `acc[j]` for each later row j in p's list. The work becomes
    Σ_p n_p²/2 instead of N²/2 × (tuples per row): 14x to 71x less work on homfam_aldosered, p450,
    sdr and Acetyltransf. This is a sparse join. On homfam_rvp the saving is only 1.3x, so rvp kept
    the list path until kind 4.
  - *Kind 4, element bit masks + index (v9).*
    - Write the count c(p) as the set of "elements" (p, 1), …, (p, c(p)). Then min(c_i, c_j) is the
      size of the intersection of the two element sets.
    - The ≤ 1024 most frequent elements become bits of a per-row mask, so their contribution is
      `AND` + `POPCNT` per pair. This is bit-parallel set intersection.
    - What is left comes from an inverted index as in kind 3.
    - On homfam_rvp, 128 elements leave 0.09 index steps per pair, against 67 for kind 3.
- **The conversion is fused into the pair threads.**
  - Each count is stored as a double, then converted by the original statements in the same order.
  - `lenfac` is read from a table filled by a verbatim copy of MAFFT's `lenfac` statement whenever the
    lengths span fewer than 512 values.
  - This removes the serial O(N²) pass and three of four divisions per pair.
- **The distance half-matrix lives in one anonymous `mmap` (v9).**
  - It was one `calloc` per row, zeroed serially by glibc.
  - With transparent huge pages, the kernel zero-fills it on first touch inside the threaded pair
    loop, and the guide tree's strided reads hit 2 MB pages instead of one 4 KB TLB entry per value.

### Why the output is identical

- **Counts.** Every path returns exactly the integer Σ_p min(c_i(p), c_j(p)), including MAFFT's
  quirks:
  - 6-mers span dropped ambiguity letters;
  - sequences shorter than 6 kept letters give count 0;
  - the diagonal is the tuple count.
- **Distances.** The same double expression is applied to the same doubles. `otool` shows the same
  instruction sequence as the original loop: fdiv, fadd, fdiv, fmadd, fadd, fdiv, fdiv, fsub, fmul,
  fadd.
- **Order.** Each matrix cell is written once by one thread.
- **Check.** `MAFFT_SEXTET=check` reruns `commonsextet_p` plus the original conversion for every pair
  and compares the stored 8 bytes. It passed on every set tried, at every SIMD level, x86 included.

### Measured

| what | before | after |
|---|---|---|
| DNA pair kernel (M5, NEON) | 2214 ns/pair | 20 ns/pair |
| sim5k stage, 1 thread (M5) | 29–33 s | 0.27–0.31 s |
| homfam_p450 stage, 4 threads (M5) | 48.4 s (baseline), 16.4 s (v3) | 1.29 s |
| homfam_rvp pair loop (CARC) | 7.9–8.7 s (list path, v8) | 0.94–1.09 s (kind 4) |
| homfam_rvp half-matrix allocation (CARC) | 3.7 s | 0.02 s |

- On x86, v3's portable fallback took 1.7 s for rnasim10k at 32 threads on CARC, about 3300 cycles
  per pair. llvm-mca puts the AVX2 kernel at 96–128 cycles per pair (Zen 2/3).

---

## 2. The UPGMA guide trees

### What MAFFT does

`fixed_musclesupg_double_realloc_nobk_halfmtx_memsave` is greedy agglomerative clustering with
nearest-neighbour pointers. Every node i keeps `nearest[i]` and `mindisfrom[i]`. Each of the N−1
steps:

1. scans the active nodes for the smallest `mindisfrom` (strict `<`, lowest index wins, the last
   list node excluded by a loop quirk);
2. merges that pair (im, jm);
3. computes each active node's new distance to the merged cluster with
   `cluster_mix_double`: 0.9·min(d1, d2) + 0.05·(d1 + d2), one fmadd;
4. runs a full `setnearest()` rescan for every node whose nearest neighbour was im or jm.

The work is O(N²) per tree in the typical case, and nearly every value is a strided read from a
different row of the half-matrix. That is one cache line per value, and on 4 KB pages one TLB miss
per value.

### What v9 does

- **First rewrite (overhead lever, the `_mt` version).**
  - The active nodes live in an ascending array instead of a linked list, with the same visiting
    order.
  - The update loop and rescans run on a small thread pool. Chunks of the scan order are claimed
    through one atomic word.
  - Each chunk keeps MAFFT's running minimum: strict `<`, starting from 999.9. The chunk minima are
    then combined in scan order with strict `<`, which is exactly the serial running minimum and its
    first-index-wins tie rule.
- **v9 rewrite (`_v9`).**
  - *setnearest is a lexicographic minimum.* MAFFT's rescan of node i walks the row part (j > i,
    ascending) and then the column part (j < i, ascending) with strict `<`. That returns the
    lexicographic minimum of (d(i, j), rank_i(j)), where rank_i(j) = j for j > i and j + N for j < i.
    Everything below relies on this identity.
  - *Candidate lists as certificates.*
    - Every node keeps its 16 lexicographically smallest entries from its last full scan, plus a bound
      θ: the 17th entry.
    - Distances change only in the update loop. When node i's distance to the merged node drops
      below θ_i, that node joins i's list; when a list overflows it is dropped.
    - A rescan reads only the listed distances. If their minimum is below θ it is provably the minimum
      over all active nodes. Otherwise a full threaded scan runs and rebuilds the list.
    - On homfam_rvp, 95–97% of about 670,000–700,000 rescans per tree were settled by the lists,
      with 16 reads each.
  - *Initialization in two sequential sweeps.* A sweep over column bands, then one over rows, instead
    of N column-strided scans.
  - *Update loop on a thread pool.* Software prefetch runs 12 nodes ahead. A phase uses
    √(active/50) threads, because synchronizing 32 cores costs about 1 µs per thread per step.
  - *The next pair to merge* comes from per-chunk minima combined in index order.
  - *Lazy deletion*: no memmove per merge.
  - *Progress output* is throttled.

### Why the output is identical

- Every distance comes from the same `cluster_mix_double` call on the same two values in the same
  loop body.
- `nearest`/`mindisfrom` receive the values the original statements give: the list shortcut returns
  exactly setnearest's lexicographic minimum, and parallel minima are combined in scan order with
  the original tie rule.
- `MAFFT_UPG=check` builds each tree a second time with the `_mt` version, which was itself checked
  bit for bit against v7.526's function by the overhead lever's harness. It compares `topol`, `len`,
  `dep` and the matrix left behind.
- Stress builds with tiny parameters made the lists fail often and the parallel phases kick in on
  small inputs (1 to 10 threads). The check passed in every configuration.

### Measured

| what | before | after |
|---|---|---|
| sim5k, both trees, 10 threads (M5) | 0.36 s | 0.15 s (`_mt`) |
| sim10k, both trees, 10 threads (M5) | 1.62 s | 0.63 s (`_mt`) |
| homfam_rvp, each tree (CARC) | 122 s and 119 s (v8: 90% of its 268 s run) | 7.2–9.4 s |

- The `_mt` version stops scaling at 3–4 threads because it is memory-bound.

---

## 3. Progressive alignment

The two progressive passes do one profile-profile alignment per guide-tree node. For DNA most steps
go through `Falign`:

1. FFT cross-correlation of the two profiles;
2. anchor search;
3. a Gotoh DP (`A__align`, or `G__align11` for one sequence against one) on each segment between
   anchors.

Steps whose groups have more sequences than the sequences are long skip the FFT and run one
full-size DP over the whole alignment width. On sim5k the largest such step aligns 9,887 × 8,225
columns.

### 3.1 The DP kernels

**What MAFFT does.** `A__align` fills the DP one row at a time.

- **`match_calc` per row:** 26 × 26 fmadds to form per-letter row scores, a `calloc`/`free` per row,
  then a sparse sum for each column.
- **Per cell:** four fmadds for the gap candidates and four comparisons that pick the gap states,
  the move and its traceback code.
- **Why it is slow:** the row-recurrence variables `mi`, `mpi` and the pointer `m` are `static`
  thread-local, so every cell loads and stores them through TLS memory, and every comparison is a
  data-dependent branch.
- **Cost:** 6.0 ns per cell on the M5, which is 86% of the progressive time.

**What v9 does.**

- **Fused row kernel (`Adp_row_fused4`).**
  - One pass per row computes the match score, the column recurrence (`m[]`, `mp[]`), the row
    recurrence (`mi`, `mpi`, now in registers) and the traceback code, using branch-free selects with
    the same comparisons in the same order.
  - *Match scores.* For nucleotides, the match chain runs over a fixed list of ≤ 4 letters in
    ascending order, padded with exact `0.0` terms (rule 3).
  - *The `+ fpenalty_ex` trick.* MAFFT computes `mi = max(...) + fpenalty_ex` with
    `fpenalty_ex == 0.0`. Adding +0.0 changes nothing except −0.0, which `mi` can no longer be after
    column 1. So `g2 + 0.0` is computed off the loop-carried chain and the select uses it directly.
    This leaves only compare + select on the serial dependency.
- **Fused `G__align11`.** The 1-vs-1 DP gets the same treatment, with match scores read inline.
- **Cheaper setup.**
  - Row scores are summed over the letters actually present.
  - Profiles and gap counts are built in one pass over the rows (`profile_fused`).
  - The first row and column are formed densely.
  - Gaps are inserted in place.
- **gcc code generation (gccpar lever).** gcc compiled the branch-free C as branchy code. Its
  if-conversion kept integer selects whose arms are expressions (`i-1`, `j-mpi`) as conditional
  blocks, which mispredict about once per cell.
  - The integer selects are now written as bit arithmetic, `ISEL(C, A, B) = B ^ ((B ^ A) & -C)`,
    which gcc folds into a straight-line select (`csel`/`cmov`).
  - The "horizontal gap wins" block became a double select plus an integer mask.
  - On x86, every `c >= x ? c : x` is spelled `x > c ? x : c`. That has exactly the semantics of one
    `maxsd` (return the second operand unless the first is greater), and it selects the same operand
    for every non-NaN input, signed zeros included.
- **SIMD across cells with a prefix scan (prog3 lever).**
  - Everything in a cell except the row carry (`mi`, `mpi`) is independent across columns, so x86
    runs 4 cells per AVX2 vector (2 per SSE2 vector). Each lane runs the scalar statements of its own
    cell.
  - The carry update `if (g >= mi) { mi = g; mpi = j-1; }` is a *selection*: "keep the latest element
    of maximal value". That selection is **associative**: combine(a, b) = b if b.v ≥ a.v else a. Both
    groupings of (a, b, c) return the latest maximal element, because `>=` is a total preorder on
    non-NaN doubles, and ties, including +0.0 versus −0.0, go to the later element.
  - So an in-register Hillis–Steele prefix scan (shifts by 1 and 2), combined with the incoming carry,
    gives every lane exactly the `(mi, mpi)` the scalar loop has at that cell. Values are only
    selected, never recomputed, so their bits match.
  - Indices travel as doubles, which are exact below 2^53.
  - ARM keeps the scalar fused loop, because NEON lanes were not faster on the M5.
- **Rows with more than 4 letters** (IUPAC rRNA, protein).
  - The match scores of a block of columns are formed first, as one fmadd chain per column over the
    letters present in that block, in ascending order and with exact zero terms.
  - Four blocks run with independent accumulators, because a single chain is limited by add latency.
  - The cell loop then reads the scores back.
  - On CARC these rows cost only 3–6% more per cell than 4-letter rows.

**Measured.**

| what | before | after |
|---|---|---|
| A__align DP fill, sim5k pass 1, 1 thread (M5, clang) | 40.4 s (6.0 ns/cell) | 10.0 s (1.5 ns/cell) |
| DP fill under gcc (rnasim2k) | 3.8–3.9 ns/cell (v3's kernel) | 1.7 ns/cell |
| G__align11 under gcc | 1.9 ns/cell | 0.9 ns/cell |
| rnasim5k progressive passes vs the gcc baseline, 1 thread | 1.90x (v3) | 3.34x (gccpar) |
| AVX2 4-letter row, llvm-mca Zen 3 | 12.4–13.1 cycles/cell (scalar) | 5.5–6.2 cycles/cell |

- The remaining gap between gcc and clang is the FMA contraction that gcc legitimately does not do.
  The gcc baseline pays it as well.

### 3.2 The scheduler

**What MAFFT does.** Threads take guide-tree steps strictly in merge order.

- A thread that takes step l waits, still holding it, until both children are done. Independent
   subtrees can only overlap within a window of N consecutive step indices: head-of-line blocking.
- Every completion wakes every waiting thread.

**What v9 does.** It uses list scheduling by critical path.

- A step enters a max-heap when its last child finishes, and idle threads pop the step with the
  largest estimated remaining critical path.
- A step's priority is its own cost estimate summed over itself and all its ancestors up to the root:
  - originally size^0.75;
  - since v8, size × Lmax + 3 × L0 × L1, using only sequence lengths known before the pass.
- No thread ever waits while holding a step.

**Why identical.** Every step is a pure function of its children (rule 5).

**Measured.**

- A replay of measured step times gives the ready queue 15–25% at 10 threads (M5).
- At 32 threads, the v8 priority brings 16S.T pass 1 from 1.83–1.86 s to 1.51–1.54 s. A perfect
  oracle would reach 1.46–1.51 s.
- The baseline profile already showed that order is not the whole story. With W the total work, CP
  the critical path and P the thread count, T ≥ max(W/P, CP), and the average parallelism W/CP of
  these trees is only about 7 to 15. At 32 threads the passes are bound by the critical path, which is what
  §3.3 and §3.4 attack.

### 3.3 Idle threads help inside large steps

- **Helper pool (`dphelp.h`).** A thread with no ready step joins a "posted job" of a running step.
  - Jobs include DP tiles, profile columns, gap insertion by rows, FFT input vectors, anchor profiles
    and row copies.
  - Every job is split so that each value is computed by exactly the statement the one-thread code
    uses, just on another thread.
- **Tiled wavefront DP.**
  - A large DP (≥ 2 M cells, with idle threads around) is cut into up to 32 column bands. Each band
    walks down the rows and needs, from its left neighbour, the row carry after the neighbour's last
    column and that column's value from the previous row. The neighbour publishes these before
    advancing.
  - This is a pipelined wavefront over anti-diagonals: a divide-and-conquer of the DP table by
    columns that respects its dependencies.
  - Every cell is computed by the same macro with the same operands, in the same two alternating row
    buffers. Even the unused stale value MAFFT keeps in the last column evolves as before.
- **v8 tuning for 32 Zen cores spread over sockets.** The first version ran at about 18.5 ns per cell
  per thread on CARC, against 2.8 ns serially. The fixes:
  - **Sticky bands.** 53% of claims used to change band; now a thread follows its band and hands
    rows over in batches.
  - **No shared cache lines between bands.**
  - **Band-major traceback.** Each band writes its traceback codes to its own contiguous region
    instead of a slice of every row. Those scattered writes were about 2/3 of the tiled thread time.
  - **Letter selections in parallel**, and narrower bands.
  - **Robust waiting:** spin, then yield, then sleep.
- **Parallel O(rows × columns) passes inside a Falign step:** FFT input vectors, anchor profiles, the
  two forward FFTs at once, and row copies.
- **In-place single segment.** 9 of 10 Falign calls on 16S.T's critical chain end with one segment
  covering both groups. That segment is now aligned in place instead of being copied out, aligned
  and appended back. The same bytes reach the same A__align call.

**Measured.**

- On the M5 (10 threads), the tail of each pass, with 1–2 steps in flight, went from 1.4–1.9 s to
  about 0.5 s.
- On CARC, 16S.T's tiled row work fell from 40.7 to 21.1 thread-seconds. The progressive passes went
  from 4.0–4.1 + 3.7 s (v7) to 2.0–2.1 + 1.9 s (v8), and the whole run from 8.9 s to 5.1 s.
- The baseline built with aggressive flags took 57.4 s on 16S.T.

### 3.4 Memory traffic, copies and lazy gap insertion

- **Allocation churn.**
  - The traceback matrix (up to (1.3 L + 110)² ints) was freed and re-zeroed whenever a larger
    alignment came along.
  - The Falign buffers and `mseq` were `calloc`ed on every call.
  - All of these now live in per-thread arenas, and big tracebacks go back into a small shared pool.
  - Page faults per pass fell from 480 K to about 230 K (M5, 10 threads).
- **No row copies.**
  - MAFFT copied every member row of a step into private buffers and back, because another thread's
    "Reallocating.." could move every row under the global lock.
  - Now no thread moves a row it does not own. A step grows only its own rows (to `alloclen + 11`
    bytes, the size MAFFT's global realloc gives) and aligns them in place.
  - Concurrent steps never share a row, and a parent starts only after its children are marked done
    under the mutex.
- **Deferred ("lazy") gap insertion.**
  - Near the root, a step that adds one sequence to a group of 9,999 had to insert its gaps into
    9,999 rows of about 10,400 columns: 104 MB moved, about 45% of the step.
  - A group with a cached profile does not need its rows until a later step reads them. So A__align
    returns the gap table instead of applying it.
  - Pending tables are composed per row, which is pure position bookkeeping, and applied once, in
    parallel, when the rows are next needed. Since v8 the gaps are inserted on the way into each
    row's grown buffer in one pass.
  - Inserting gaps only places characters, so the rows end up exactly as if every table had been
    applied immediately, and the DP only ever sees the cached profile and the lengths.
- **Small fixes.**
  - Falign appended segments with `strcat`, which is quadratic. It now appends at a known offset.
  - An unused label string stops calling `sprintf` once it reaches the 100 characters it can hold.
  - Atracking no longer fills a row and a column it never reads.

### 3.5 Falign's anchor search

- MAFFT rebuilt the per-column profiles of both groups for every candidate FFT lag.
- They are now built once per call, member by member in the original order. The per-lag score sums
  use the same terms in the same order.
- The running window sum that decides segment boundaries is therefore bit-identical
  (`MAFFT_ARCHECK`).

### 3.6 The memsave switch (an order-dependent corner of MAFFT)

When a child alignment grows past 30,000 columns, MAFFT flips a global variable (`alg = 'M'`) so that
every later step uses the memory-saving aligner. Which steps run after the flip depends on the
schedule, so MAFFT's own threaded output is racy here. In v7.526, `edge_noshort` at 10 threads gives
two different md5s from run to run. v9 handles it in three ways:

- **At `--thread 1`** v9 always uses MAFFT's in-order dispatch.
- **At more threads, up to 20,000 sequences:**
  - The ready queue runs optimistically. Until the switch fires, every step uses the same aligner in
    any order.
  - If a step would fire it, the pass is abandoned: rows, `nlen`, profile caches, member lists, done
    flags and the FFT log are restored.
  - The pass is then redone with MAFFT's own in-order dispatch and its original thread count (16 at
    `--thread ≥ 16`). v9 therefore produces the same set of possible outputs as MAFFT.
  - Normal inputs never abandon a pass.
- **Above 20,000 sequences**, see §6. There MAFFT is serial and deterministic, and v9 reproduces its
  exact switch step.

---

## 4. The distance matrix from the MSA

### What MAFFT does

For every pair i < j, `naivepairscorefast`:

1. calls `strlen`;
2. `calloc`s two buffers of the aligned length;
3. copies both rows while dropping the columns that are gaps in both;
4. scans the copies with two data-dependent branches per column;
5. frees the buffers.

On sim5k that is about 10 µs per pair over 10,939 columns, with 86% of the cells gaps. The score is then
turned into a distance with `(1.0 − S/min(self_i, self_j)) · 2.0`.

### What v9 does

- **Closed form.** S = Σ over columns where both rows have a residue of `amino_dis[a][b]`, plus
  `penal` × (number of maximal runs of (gap, residue) columns + runs of (residue, gap) columns),
  counted after the gap-gap columns are removed.
- **Gap runs by carry propagation (bit-parallelism).**
  - Take row i's gap bitmask G and seeds X = G & ¬G_j: the columns inside i's gaps where j has a
    residue.
  - In the multiword integer sum G + X, every gap block containing a seed carries out into the first
    column after the block, which is 0 in G. So popcount((G + X) & ¬G) counts exactly the runs.
  - The two directions can never end at the same column, so one popcount of the OR counts both.
  - Padding bits catch a carry off the end.
- **DNA: bit planes.**
  - The four bases get a 2-bit code such that `amino_dis` depends only on code(a) XOR code(b). For
    MAFFT's DNA matrix that is a=0, c=1, g=2, t=3: XOR 0 is a match, 2 a transition, 1 or 3 a
    transversion. All 24 codings are tried at run time.
  - The residue sum then needs 3–4 popcounts per word:

    ```
    S_res = Td[3]*(len - cnt(gor)) + (Td[0]-Td[3])*(len - cnt(gor|x0|x1))
          + (Td[1]-Td[3])*cnt(x0 & ~(gor|x1)) + (Td[2]-Td[3])*cnt(x1 & ~(gor|x0))
    ```
- **IUPAC codes (real rRNA).**
  - The 4 most frequent kinds go into the planes. Rare other letters are "exceptions": a per-row
    sorted list of (column, byte).
  - The pair's two lists are merged, and a correction `amino_dis[a][b] − Td[code(a) ⊕ code(b)]` is
    added from small tables, using the other row's class read from the planes already streamed.
  - On 16S.T, 0% of pairs now leave the fast path; in v3, every pair did.
- **Protein: the column path.**
  - For row i, every column where i has a residue walks the list of later rows with a residue there
    and adds `amino_dis` into `acc[j]`.
  - On AVX2, dense columns (≥ 1/8 of rows occupied) are processed 32 rows at a time with `vpshufb`
    table lookups.
  - Gap runs come from the gap masks as above.
- **x86 kernels.**
  - A POPCNT kernel with an `adc` carry chain.
  - An AVX2 version split into a gap-run pass, 4 pairs per vector, and a residue pass with
    nibble-table popcounts. In the gap-run pass the carry-out of the 3-input add is the top bit of
    (a & ¬b) | (a & ¬s).
- **v9 extras.**
  - Non-temporal stores of the matrix on x86 (N ≥ 4096): the rows are cold by then.
  - All-gap words are skipped in the gap-run loops when enough of them exist. That is exact, because
    such a word passes both carries through unchanged.
- **No per-pair allocation at all.** MAFFT's skip table is no longer built.

### Why the output is identical

- Every summand of MAFFT's double accumulator is an int (`amino_dis` or `penal`), and |S| stays far
  below 2^53. So MAFFT's double sum equals the exact int64 sum, which every path computes.
- The final expression `(1.0 − S/bunbo)·2.0` is unchanged and has no multiply-add to contract.
- Rows of unequal length, or bytes ≥ 0x80, fall back to the original function.
- `MAFFT_MSADIST=check` compares every pair with `naivepairscorefast`, bitwise, on every path and SIMD
  level. A matrix dump on sim5k (99,980,000 bytes) is identical to the original code's.

### Measured

| set | before | after |
|---|---|---|
| sim5k, 1 thread (M5) | 122.2 s | 1.77 s (69x) |
| sim5k, 10 threads (M5) | 19.4 s | 0.37 s (53x) |
| 16S.T, 4 threads (M5) | 119.6 s (baseline), 24.6 s (v3's generic path) | 6.7 s |
| homfam_p450, 4 threads (M5) | 344.8 s | 14.8 s |
| rnasim10k, 32 threads (CARC) | 22.2 s | 1.8 s (v3, portable kernel) |

---

## 5. I/O, `f2cl` and startup

### What MAFFT does

- The script copies the input with `cat | tr`.
- It runs about 80 small processes, 33 of them for option files that default to `/dev/null`.
- `disttbfast` reads the input twice with one locked `getc()` per byte.
- `f2cl` re-wraps the alignment at 60 columns. It reads its input three times with a locked `getc()`
  per byte (165 M locked calls for 55 MB on sim5k), computes CLUSTAL conservation marks that FASTA
  output never prints, and writes one `fprintf` per line.

### What v9 does

- **Readers.** `flockfile` once per call, then `getc_unlocked`: the same bytes in the same order, and
  still thread-safe.
- **`f2cl` in-memory fast path.**
  - It reads a regular file once and splits it at line-start `>` with `memchr`.
  - It fills the same structures exactly as the original readers do: names cut at 254 characters,
    the same DNA/protein guess.
  - Anything unusual falls back to the original code after `rewind()`: a pipe, junk before the first
    `>`, NUL, `=`, `<` or `>` in a sequence, or sequences of 5,000,000 bytes or more.
  - Marks are computed only for CLUSTAL output, and wrapped records are written with one `fwrite`.
- **The launcher script.**
  - A small helper, `crtolf`, replaces `cat | tr "\r" "\n"` for regular all-ASCII files. For pure
    ASCII, `tr` is the identity plus CR→LF in the C, UTF-8, SJIS and eucJP locales, which was checked.
    Any other input goes through `tr` as before.
  - `/dev/null` option files are created empty without running three processes each.
- **Buffers.** `disttbfast` uses 1 MB stdio buffers when stdin or stdout is not a terminal.
- **Progress output.** The guide-tree progress counter is printed less often (every 100 or 1000
  merges for large inputs). With stderr on a pipe, each write cost about 120 µs, which added up to
  1.1–1.3 s per tree on homfam_rvp. stdout is unchanged.

### Why the output is identical

- The bytes parsed and written are the same. The fast paths accept only inputs on which they provably
  match the original, and everything else takes the original code.
- `f2cl` was fuzzed with 1,500 generated inputs and option sets against the original, with
  identical stdout and exit status.
- 72 option and locale combinations of the launcher were identical.

### Measured

| what | before | after |
|---|---|---|
| `f2cl` on sim5k's 55 MB alignment (M5) | 1.73–1.82 s | 0.10 s |
| `disttbfast` start-up, sim5k (M5) | 0.31 s | 0.047 s |
| "other" (start-up + output), sim5k (M5) | 2.5 s | 0.3 s |
| "other", sim10k (M5) | 5.7 s | 0.43 s |

---

## 6. Threading above 20,000 sequences

### What MAFFT does

`disttbfast` contains `if( njob > 20000 ) nthreadtb = 0;`, so above 20,000 sequences both
progressive passes run in the serial `treebase()` on one core, whatever `--thread` says. The 2018
comment gives two reasons:

- stack depth in `topolorder_lessargs`;
- the inefficiency of the local row copies.

Above 20,000 sequences, then, a 32-core run of MAFFT spends most of its time on one core.

### What v9 does

- **The cutoff is gone.**
  - The serial and threaded paths call the same aligners with the same arguments. The only
    differences are inert under `--thread N`:
    - a profile-reuse branch that needs `nthread == 0`;
    - a `<` versus `<=` in a buffer-size test;
    - a gap character that is `-` in both paths.
  - The tree recursion that worried the authors is no longer called by any progressive thread (since
    2018/Mar). Worker stacks are max(8 MB, 200·N) bytes anyway, the same allowance MAFFT gives its
    main thread.
  - Local copies are gone (§3.4).
- **The memsave switch is reproduced exactly as in the serial pass,** with threads, in three phases:
  - *Phase 0.* The ready queue runs. A step that finds a child over 30,000 columns puts itself back
    and lowers a limit to its own index, and no step at or above the limit starts. The set of steps
    with such a child is upward-closed: alignments only grow towards the root, and parents have
    larger indices. So every step below the serial switch step s* still runs with `alg = 'A'`, as in
    MAFFT, and the limit ends at s*.
  - *Phase 1.* Subtrees that had already run past s* are reset and rerun.
  - *Phase 2.* `alg = 'M'`, then all remaining steps, in any order.
- **Memory caps for these passes.**
  - No traceback pool.
  - Per-thread arenas over 16 MB are freed after each step.
  - Falign rows are sized by the step.
  - A budget admits concurrent full-size DPs (tracebacks ≥ 64 MB) only while their tracebacks fit in
    N × alloclen bytes. That is the size of the alignment copy MAFFT's own root step allocates.
- **The 16-thread cap.**
  - MAFFT's launcher caps the progressive threads at 16 for `--thread ≥ 16`.
  - v9's launcher passes `-C N-16-N`. The order-independent passes (the ready queue, the
    above-20,000 phases) use all N threads.
  - MAFFT's 16 are kept for the in-order redo after a memsave switch at ≤ 20,000 sequences, where
    the thread count can change MAFFT's own result.

### Why the output is identical

- Every step computes the same function of its children as in the serial pass.
- The phase logic makes the switch happen at the serial pass's step s*. Three inputs that exceed
  30,000 columns were tested with the above-20,000 logic forced, against MAFFT's 1-thread result.
- MAFFT's own output above 20,000 does not depend on `--thread`.

### Measured

| set | baseline | serial passes (v4) | threaded (largen) |
|---|---|---|---|
| rnasim25k, whole run, 10 threads (M5) | 1358.7 s | | 63.5 s |
| rnasim25k, progressive passes (M5) | | 135.1 s | 40.3 s |
| 16S.B.ALL, progressive passes (M5) | | 68.6 s | 18.9 s |

On CARC with v9, homfam_rvp (93,681 sequences) took 24.3 s against 950–1220 s for the baseline, and
rnasim50k took 87.6 s. RESULTS.md has the final runs.

---

## 7. x86 SIMD dispatch

- **Three levels.**
  - 0: portable C.
  - 1: SSE2 for the 6-mer kernel and POPCNT for msadist on x86-64; NEON on arm64.
  - 2: AVX2, plus POPCNT and BMI1 for msadist.

  The level is detected once with `__builtin_cpu_supports`, before any worker thread starts.
- **Per-target functions.** Each kernel is a separate function with
  `__attribute__((target("avx2")))` (or `"popcnt"`, `"avx2,popcnt,bmi"`). A default build with no `-m`
  flags therefore runs AVX2 on Zen 2/3, and a `-march=native` build compiles the same code for AVX2
  directly.
- **Exactness.**
  - In the distance stages only integer kernels carry a target attribute. The double arithmetic
    (`lenfac`, `(1 − S/bunbo)·2`) stays in default-target code, so contraction and rounding are those
    of the baseline built with the same flags.
  - The DP vector kernels (§3.1) carry `target("avx2")` but never `"fma"`. Each lane runs the scalar
    cell's statements, so whatever the build flags make of `a*b + c` in the baseline, they make of it
    here too.
  - Because gcc's ISO C mode never contracts, a gcc `-O3 -march=native -flto` build gives the default
    build's md5 in both the baseline and v9.
- **Portability details.**
  - `<immintrin.h>` is included between `#pragma push_macro`/`pop_macro` of MAFFT's one-letter macros
    (`B`, `C`, ...), which clash with parameter names in some compiler headers.
  - Only intrinsics available since gcc 4.9–5 are used.
- **Testing knobs.** `MAFFT_SIMD=scalar|sse2|popcnt|avx2` and `MAFFT_ADPVEC=0|2|4` force lower levels.
  Every level was run with the check modes under Rosetta, and with md5 comparisons on CARC
  (VERIFICATION.md).

**Measured (llvm-mca, Zen 2/3).**

| loop | before | after |
|---|---|---|
| 6-mer pair kernel | about 3300 cycles/pair (derived from v3's portable path on CARC) | 96–128 cycles/pair (AVX2) |
| msadist DNA kernel, per 128 columns | 27–35 cycles (portable) | about 8–9 cycles (AVX2) |
| DP row, 4 letters, per cell | 12.4–13.1 cycles (scalar) | 5.5–6.2 cycles (AVX2) |

---

## 8. v10: long genomes and near-identical genome collections

v9 was built and timed on alignments of up to a few thousand columns. Two kinds of input it left at
about MAFFT's own speed:

- **Alignments wider than 30,000 columns** (long genomes, or shorter sequences with long unaligned
  stretches). Past 30,000 columns MAFFT switches every later merge to its memory-saving aligners,
  `Falign_udpari_long` and `MSalignmm`, which v9 never touched. v9 was 1.0–1.2× there.
- **Collections of near-identical genomes**, such as SARS-CoV-2. They stay just under 30,000 columns
  (29,875 for 500 genomes) and use v9's ordinary path, but each merge finds about 200 FFT anchors, and
  the merges form a deep chain (129 steps per pass for 500 genomes). On CARC, v9 was 3.4× on 500
  genomes and 4.6× on 1,000.

v10 = v9 + the changes below, in six files (`MSalignmm.c`, `Falign.c`, `fftFunctions.c`, `fft.c`,
`disttbfast.c`, `functions.h`). Every change keeps the rules of §0.4, and each new kernel has a bitwise
self-check in the normal build (§8.7).

### 8.1 Myers-Miller (`MSalignmm`)

**What MAFFT does.** `MSalignmm` is a Myers-Miller (Hirschberg) divide and conquer: a forward and a
backward DP sweep find the midpoint of the optimal path, and the two halves recurse. Its cells are
branchy (data-dependent `if`s, an `i == imid` test per cell); every row recomputes a 26×26 match table
after a `calloc`; every recursion node allocates and zeroes work arrays of (max length) × 26. It runs on
one thread. On 30 genomes of 30 kb at one thread it is 70% of the run.

**What v10 does.**

- **Match scores** come from a table per call (per row of profile 1, per letter of profile 2) and its
  transpose, built once. Each score is the same chain of `match_calc` terms.
- **Cells.** The forward and backward rows run fused, branch-free kernels: the original's four
  candidates and selects, without the jump bookkeeping that only rows `imid`, `imid−1` and `jumpi`
  need; those rows run the original statements. On x86, 4 cells per AVX2 step (2 per SSE2 step), with
  runtime dispatch; the selects are `maxpd`, and the row carry is a "later element wins ties" prefix
  maximum inside the vector. ARM keeps the scalar kernel.
- **Parallel divide and conquer.** The two halves of a node are independent; idle helper threads take
  them from a per-call task stack. The gap tables are joined in order at the end, and the score is
  summed in the original order.
- **Parallel sweeps.** In large nodes the forward and backward sweeps are split into up to 16 bands of
  at least 1,024 columns that run as a wavefront: band b computes row i after band b−1 has published
  that row's last value and row carry.
- **Output.** The recursion now produces only the two gap tables; the rows are made from them at the
  end, in parallel by rows (`gapireru` of each input row, which is what the recursion writes piece by
  piece).
- The full-DP leaf (`MSalignmm_tanni`, nodes under 100 rows or columns) and the profile set-up
  (`cpmx_calc_new`, gap counts) are ported the same way.

### 8.2 The anchor search of `Falign_udpari_long`

**What MAFFT does.** For each of up to 500 candidate FFT lags (`NKOUHO_LONG`), `alignableReagion`
rebuilds the per-column profiles of both groups, O((c1 + c2) × L) per lag. The early exit on a lag
without anchors is commented out there, so all 500 lags run. On low-divergence long genomes this is
70% of the run. `getKouho` picks the 500 best lags with 500 passes over the correlation array.

**What v10 does.**

- The profiles are built once per call, and the nonzero entries of each column are kept as compressed
  rows. The 500 lags are scored by idle helpers in parallel, each into its own segment buffer, and the
  segments are taken over in lag order, exactly as the original loop finds them (with the same
  `TOO MANY SEGMENTS` check).
- `getKouho_fast` picks the 500 best lags in one pass with a bounded heap. Its order (value
  descending, index ascending) and its padding when fewer than 500 values qualify are those of the
  original's 500 rounds.
- No row copies for the FFT input (the rows were only read); the input vectors are built by column
  ranges, and the two forward transforms run at the same time.
- A single segment covering both groups is aligned in place. Otherwise the segment results are
  appended at their offsets (not with `strcat`, which is quadratic) in a per-thread arena, instead of
  four `clus × alloclen` matrices allocated per call.
- The original stays, as `Falign_udpari_long_orig`, for constraints and scoring-matrix sets.

### 8.3 Threads in MAFFT's in-order dispatch

Once a pass switches to the memsave mode it runs in MAFFT's in-order dispatch (§3.6). There, a thread
that finds every step handed out exits, so the largest merges, near the root, ran on one thread.

- Such threads now stay as helpers until the root step is done. They help inside the steps still
  running (§8.1, §8.2).
- At `--thread` above 16, the in-order dispatch also starts helper-only threads (MAFFT's dispatch uses
  16 threads there). They never take a step.
- Which thread takes which step, in which order, and when each step checks the switch, is MAFFT's.

### 8.4 When the switch fires: a race in MAFFT, and the fix

In MAFFT's in-order dispatch a thread claims the next step under the mutex, checks there whether a
child is wider than 30,000 columns (if so it sets `alg = 'M'`), copies the member rows, unlocks, and
only then reads `alg` to choose `Falign` or `Falign_udpari_long`. When the switch fires at one of the
first steps, whether the steps already claimed read 'A' or 'M' depends on timing. MAFFT's row copies
inside the mutex serialize the claims, so one outcome dominates; v9 moved the row copies out of the
mutex, which made the others common. On the development Mac, at `--thread 10`, `long500x30k` then gave
an md5 different from MAFFT's (tallies of the first 8 steps of pass 1 over repeated runs: MAFFT
`AAAMMMMM` in 14 of 22 runs, v10 before the fix `MMMMMMMM` in 8 of 17).

**The fix.** Each step reads `alg` under the mutex, right after the switch check. Because `alg` only
goes from 'A' to 'M', and only under this mutex, that value is one MAFFT's own step can read; it is the
outcome MAFFT's serialized claims make dominant (after the fix: `AAAMMMMM` in 5 of 5 runs). And a pass
with an input sequence longer than 30,000 starts in the in-order dispatch at once, as MAFFT does,
instead of first trying the ready queue (the switch is certain there).

**What remains: a switch in the middle of a pass.** When the switch fires late in a pass, while
other steps are in flight (the `viral100` set, or 2,000 SARS-CoV-2 genomes), which of those steps read
'M' at the top level, and which of their Falign segments (each segment reads `alg` again) use the long
aligner, depends on timing. MAFFT itself then varies: on CARC, `viral100` at 32 threads gave two
different outputs in four runs. Usually, though, one outcome wins on a given machine, because MAFFT's
steps near the switch take about the same time from run to run (`viral100` at 8 threads: 8 of 8 runs
alike; `sarscov2_2000` at 32: 5 of 5). fast-mafft keeps MAFFT's dispatch and its per-segment reads, so
the same mechanism decides its output, and every run switches at the same step as MAFFT. But it runs
those steps 20–200 times faster, by different factors for different steps, so the race has different
odds: `sarscov2_2000` at 32 threads gave MAFFT's usual output in 11 of 16 runs, `viral100` at 8 threads
in 6 of 8. Reproducing MAFFT's usual outcome every time would mean reproducing MAFFT's relative step
timings, which a faster program cannot do. At `--thread 1` there is no race and the output is
identical. VERIFICATION.md §3 has the measurements.

### 8.5 Near-identical genomes: the ordinary Falign path

One step of the critical chain on 500 SARS-CoV-2 genomes took about 38 ms on the Mac at 10 threads:
three FFTs of 32,768 points (~5 ms), `blockAlign2` (~6 ms) and about 200 small segment DPs run one
after another (~21 ms).

- **`blockAlign2`** (the anchor chaining DP). `permit()` always returns 0, so in its two inner k-loops
  every k other than 0 is skipped unless j (or i) is the last cut point. The loops now stop after
  k = 0 in that case: the same k, in the same order, with the same comparisons. O(ncut²) instead of
  O(ncut³).
- **The long FFT** (`fft.c`, 16,384 points or more). After the bit reversal, each stage's butterflies
  run block by block, with the stage's twiddles copied into two contiguous arrays, instead of twiddle
  by twiddle with a stride of 2k through the whole vector. Same butterflies, same operands; the
  butterflies of a stage touch disjoint elements, so their order inside a stage does not matter.
  32,768 points: 0.80 → 0.26 ms on the Mac. Below 16,384 points the original order is kept.
- **Parallel segments.** A call with at least 16 anchors aligns its segments in parallel by default,
  largest DP first. Each segment's inputs are columns of the original rows and the gap states at its
  cut points, so it can be aligned anywhere, any time, into its own buffer. The outputs are
  concatenated in segment order straight into the rows, and the score is summed in segment order.
  (v9 had this code but off by default, because of its memory; it now uses the per-thread row arena.)
- **The first two candidate lags at once.** Falign scores its FFT lags one by one and stops at the first
  that yields no segment; on these genomes that is the second. With an idle helper the first two are
  scored at the same time, the second by the helper, and the loop takes their segments over where it
  would have scored them.
- `getKouho_fast` (§8.2) also replaces Falign's `getKouho`.

### 8.6 Memory: a pool for the large per-thread arenas

Each thread kept its Falign row arena and its member-row block at the size of the largest step it had
run. On the critical chain the big steps run one after another on whichever thread is free, so every
thread ended up holding arenas for the largest groups: peak RSS on CARC was 2.4 / 3.3 GB on 500 / 1,000
SARS-CoV-2 genomes, against MAFFT's 0.8 / 1.3 GB. Now a thread gives an arena of 4 MB or more back to a
small shared pool when its step ends, and the next step, on any thread, takes the best-fitting pooled
block. At most two blocks stay pooled; the pool is emptied when a pass ends.

### 8.7 Why the output is identical, and the self-checks

- **Match scores** (§8.1): the table entries keep the original statements and orders. Terms left out
  (letters of zero weight) or added (padding to 4 letters) are `x * 0.0` added to a sum that starts at
  +0.0 and cannot become −0.0, so they are exact no-ops (§0.4, rule 3).
- **Cells**: the same four candidates, each `a + b * c` one statement, compared in the same order with
  the same strictness. `maxpd` returns its first operand only when it is strictly larger, which is the
  scalar select lane by lane. The row carry is a running "later wins ties" maximum, which is
  associative on non-NaN values, so the in-vector scan gives every lane the scalar loop's value bit for
  bit (values are selected, never recomputed).
- **Bands and halves**: a band reads only its own columns plus its neighbour's published value and
  carry; no column is written by two bands. A node's halves depend only on its own sweeps; the gap
  tables and the score are joined in the original order.
- **Anchors, lags, blockAlign2, FFT, segments**: the same terms in the same order; the parallel pieces
  are consumed in the original order (§8.2, §8.5). The skipped `blockAlign2` iterations did nothing.
- **Scheduling** (§8.3, §8.4): the claims, their order and the switch check are MAFFT's; extra threads
  only help inside steps.
- **Pool** (§8.6): only where buffers live changes.

Self-checks in the normal build (each exits with a `MISMATCH` message on any difference):

| variable | what it recomputes the original way |
|---|---|
| `MAFFT_MSCHECK=1` | every `MSalignmm` call (score, both gap tables, every row) and the fused profiles |
| `MAFFT_MSKCHECK=1` | every recursion node's sweeps (midpoint, jumps, split) |
| `MAFFT_ARLCHECK=1` | every lag's segments, the lag lists, the FFT input vectors, the helper-scored lags |
| `MAFFT_BA2CHECK=1` | every `blockAlign2` call |
| `MAFFT_FFTCHECK=1` | every transform of 16,384 points or more |

Switches for A/B tests: `MAFFT_MSFAST=0` and `MAFFT_UDLV10=0` (the original memsave code),
`MAFFT_PSEG=0` (sequential segments), `MAFFT_LAG2=0`.

### 8.8 Measured

Final CARC evaluation (RESULTS.md §3), `--thread 32`, mean of 3 runs, against MAFFT v7.526 built with
`-O3 -march=native -flto`:

| input | MAFFT (s) | v9 | v10 (s) | v10 speedup |
|---|---:|---|---:|---:|
| sarscov2_500 (500 real genomes) | 28.7 | 3.4x (v9's path, job 13159523) | 2.5 | 11.6x |
| sarscov2_1000 | 48.6 | 4.6x (v9's path, job 13159523) | 4.1 | 11.8x |
| sarscov2_2000 (switch in pass 1) | 3,810 | not run | 18.6 | 205x |
| edge (59–66 k columns) | 79.8 | 1.2x (Mac, 10 threads) | 5.4 | 14.9x |
| long100x30k | 326 | 1.1x (Mac, 10 threads) | 9.0 | 36.4x |
| long1000x30k | 940 | | 23.1 | 40.8x |
| long500x30k | 1,633 | | 39.8 | 41.2x |
| long300x50k | 2,061 | | 28.5 | 72.8x |

The changes of §8.5 also apply to ordinary inputs: on the 12 datasets of the v9 evaluation, v10's
mean speedups are 10.1x–43.9x against v9's 9.7x–36.2x, and its slowest run 9.98x against 8.3x.
Single steps, on the Mac at 10 threads: a late chain step of sarscov2_500 from ~38 ms to ~5.3 ms; the
critical path of its pass 1 from 4.91 s to 1.16 s.

---

## 9. Bugs found and fixed along the way

1. **The memsave switch under the ready queue** (found by the sextet lever, fixed in prog2 and, above
   20,000, in largen).
   - The first ready-queue scheduler ran steps out of index order, so the `alg = 'M'` switch at
     30,000 columns took effect at different steps than in MAFFT.
   - On `edge.fa` / `edge_noshort.fa` (59–66 k columns) that build differed from the baseline even at
     `--thread 1`.
   - The fixes are in §3.6 and §6.
2. **The v3 deadlock** (found by the gccpar lever). It was a lock-order inversion:
   - A step that reallocated all rows held the global mutex and waited for the `aseqlock` write lock.
   - Another step's helper-assisted row copy held `aseqlock` for reading and waited for the global
     mutex inside `dphelp_for`.
   - 2 of 5 runs of unfixed code on rnasim10k and 16S.T at 10 threads hung.
   - First fix: copy rows on the calling thread. Later the row copies and the global reallocation
     were removed altogether (§3.4), so the default build no longer takes these locks.
3. **The `dphelp_finish` race** (found by the v8s lever, present since prog2).
   - The poster read the job's `slots` count without the mutex.
   - A helper decrements `slots` before it increments `refs`, so the poster could see 0 and 0,
     return, and free the job's stack frame while the helper was about to call `job->run` from it.
   - It was seen once (SIGBUS) in about 110 runs.
   - Fix: a posted job is always withdrawn under the mutex.
4. **`posix_memalign` on glibc** (found on CARC by the v8s lever).
   - Under `-std=c99`, glibc declares `posix_memalign` and `nanosleep` only with a feature macro.
     gcc 13 warned about the implicit declaration; gcc 14 and later reject it.
   - `Salignmm.c` now defines `_GNU_SOURCE` before its includes, as `disttbfast.c` already did.
5. **A thread-local variable in the tiled DP's band-0 row setup** (found by `MAFFT_TILECHECK`, prog3
   lever). Fixed before merging.
6. **The aligner read after the unlock** (found by the v10 lever on the Mac, `--thread 10`).
   - v9's shorter critical section changed which outcome of MAFFT's own race (the steps claimed just
     before the memsave switch read 'A' or 'M') was common, so `long500x30k` gave an md5 that MAFFT
     rarely gives. CARC's 32-thread runs happened to match.
   - Fix: the step's aligner is read under the mutex (§8.4).

Behaviour that fast-mafft deliberately keeps from MAFFT:

- The racy memsave switch at ≤ 20,000 sequences with several threads (§3.6), when it fires in the
  middle of a pass (§8.4).
- The nondeterminism of MAFFT's iterative modes (`--auto` on small inputs, `--maxiterate`) at
  `--thread` > 1. v7.526 itself gives different md5s from run to run there.

## 10. Known limitations

- **Alignments that cross 30,000 columns in the middle of a pass**, at more than one thread: MAFFT's
  own output depends on thread timing there, and fast-mafft's race has different odds, so it does not
  always give the output MAFFT usually gives on that machine (§8.4; 11 of 16 runs on 2,000 SARS-CoV-2
  genomes at 32 threads). At `--thread 1` it is identical.
- **Small inputs.** At 1,000 sequences fast-mafft is 5.0–5.6x faster at 32 threads (MAFFT takes
  6–10 s), and 2.4x on prot2k, where MAFFT takes 1.5 s; at 2,000 sequences 7.0–7.7x (RESULTS.md §3.2).
  The fixed costs (launcher, input, the serial parts of each stage) dominate there.
- **Peak memory** is close to, but on several sets above, the baseline's on CARC (v10, RESULTS.md §3).
  - Large sets: homfam_rvp 33.20 vs 33.00 GB, homfam_sdr 10.26 vs 9.96 GB, homfam_Acetyltransf 8.66
    vs 8.47 GB, 16S.T 4.02 vs 3.69 GB; lower on rnasim10k (3.68 vs 5.81) and rnasim20k (5.19 vs 12.93).
  - Long genomes: 0.8–1.5 GB against 0.4–1.0 GB on the simulated sets; sarscov2_500 1.57 vs 0.80 GB.
  - Small sets: 0.8–2.5 GB against 0.3–1.8 GB (per-thread buffers at 32 threads).
  - All are well inside the 64 GB of an evaluation node. The excess comes from several large steps
    being in flight at once, and from per-thread buffers.
- **Tuning.** Everything x86-specific was tuned on CARC's EPYC 7513 (Zen 3) nodes and estimated for
  Zen 2 with llvm-mca. The ARM build keeps the scalar DP kernels.
- **Other strategies.** Only FFT-NS-2 (the default) was optimized and timed.
