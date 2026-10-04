/* work/prog3: deferred ("lazy") gap insertion for the progressive pass.
 *
 * When a guide-tree step aligns two groups with A__align, the gaps it opens in group s are
 * normally inserted into every row of that group right away (O(rows x columns) of memory traffic).
 * A group whose profile is cached (cpmxhist, made by A__align for every result with more than 20
 * sequences) does not need its rows again until a later step reads them (FFT path, a memsave
 * aligner, or no cached profile) or the pass ends.  So treebasethread may ask A__align to *return*
 * the gap table of such a group instead of applying it; the table is attached to the group's
 * tree node, and the rows are brought up to date ("materialized") only when needed, by applying
 * all deferred tables of a subtree at once (composed into one table per row, see lazygap_* in
 * disttbfast.c).  Inserting gaps is pure character placement, so the rows end up exactly as if
 * every table had been applied when it was made; the DP only ever sees the cached profile and the
 * lengths, which are the same either way.
 *
 * Protocol: the caller fills A_lz (on = 1) right before A__align; A__align reads and clears it.
 *   stale[s]  rows of group s are not current (only the logical length len[s] is): A__align must
 *             not read or write them (the caller guarantees a cached profile for that side)
 *   defer[s]  return group s's gap table in outgt[s] (malloc'ed, NULL if no gap was inserted)
 *             instead of inserting it (implied by stale[s])
 *   outlen    length of the alignment (the new logical length of both groups)
 */
#ifndef LAZYGAP_H
#define LAZYGAP_H
typedef struct
{
	int on;
	int stale[2];
	int len[2];
	int defer[2];
	char *outgt[2];
	int outlen;
	int used;      // set by A__align when it honoured the request
} a_lazy_t;
extern TLS a_lazy_t A_lz;
#endif
