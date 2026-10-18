#!/usr/bin/env python3
# diag_mix.py TREE : patches a MAFFT or fast-mafft source tree in place with stderr diagnostics for a
# memsave switch in the middle of a pass (stdout unchanged):
#   "SWITCH l=<step>"            when a step sets alg = 'M' (treebasethread)
#   "TOP l=<step> A|M"           the step's top-level aligner (Falign / Falign_udpari_long)
#   "MIXED l=<step> nA=.. nM=.." a Falign call whose segments were aligned partly with 'A' and partly with 'M'
#                                (the per-segment read of alg met the switch made by another thread)
# Works on the baseline tree and on fast-mafft-src (its parallel-segment path counts as all 'A').
import sys
t = sys.argv[1] + "/core/"
def rd(f): return open(t + f, encoding="latin-1").read()
def wr(f, s): open(t + f, "w", encoding="latin-1").write(s)

s = rd("Falign.c")
s = s.replace('#include "mltaln.h"\n', '#include "mltaln.h"\nTLS int diag_curstep = -1; // DIAG\n', 1)
a = s.index("double Falign( int **whichmtx, double ***scoringmatrices, double **n_dynamicmtx,")
b = s.index("\n{\n", a) + 3
s = s[:b] + "\tint diag_na = 0, diag_nm = 0; // DIAG\n" + s[b:]
c = s.index("\t\tswitch( alg )\n", b)
s = s[:c] + "\t\tif( alg == 'M' ) diag_nm++; else diag_na++; // DIAG\n" + s[c:]
e = s.index("\treturn( totalscore );\n}", c)
s = s[:e] + '\tif( diag_na && diag_nm ) reporterr( "\\nMIXED l=%d nA=%d nM=%d\\n", diag_curstep, diag_na, diag_nm ); // DIAG\n' + s[e:]
if "totalscore = Falign_psegs(" in s:
    p = s.index("totalscore = Falign_psegs(", b)
    q = s.index("\n", p) + 1
    s = s[:q] + "\t\tdiag_na += count-1; // DIAG\n" + s[q:]
wr("Falign.c", s)

s = rd("disttbfast.c")
s = s.replace('#include "mltaln.h"\n', '#include "mltaln.h"\nextern TLS int diag_curstep; // DIAG\n', 1)
a = s.index("static void *treebasethread(")
b = s.index("\n}\n", a)
body = s[a:b]
o1 = "pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );"
o2 = "pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );"
assert body.count(o1) == 1 and body.count(o2) == 1, (body.count(o1), body.count(o2))
body = body.replace(o1, 'diag_curstep = l; reporterr( "\\nTOP l=%d M\\n", l ); ' + o1)
body = body.replace(o2, 'diag_curstep = l; reporterr( "\\nTOP l=%d A\\n", l ); ' + o2)
n = body.count("alg = 'M';")
assert n >= 1, n
body = body.replace("alg = 'M';", "alg = 'M'; reporterr( \"\\nSWITCH l=%d\\n\", l ); // DIAG")
s = s[:a] + body + s[b:]
wr("disttbfast.c", s)
print("ok", n)
