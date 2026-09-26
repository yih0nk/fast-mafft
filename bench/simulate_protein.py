# Simulate tree-structured protein sequences (substitutions + short indels) for MAFFT timing.
# Protein counterpart of simulate.py: same random-recursive-tree scheme, 20-letter alphabet,
# residues drawn from approximate natural amino-acid background frequencies.
# Usage: python3 simulate_protein.py N LENGTH [SEED] > file.fa
import random, sys
n, L, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1
rng = random.Random(seed)
AA = "ARNDCQEGHILKMFPSTWYV"
# Background frequencies (percent), roughly the Robinson & Robinson / BLOSUM62 composition.
FREQ = [7.8, 5.1, 4.3, 5.3, 1.9, 4.3, 6.3, 7.4, 2.2, 5.1, 9.1, 5.7, 2.2, 3.9, 5.1, 7.1, 5.8, 1.3, 3.2, 6.4]
def draw(k=1):
    return "".join(rng.choices(AA, weights=FREQ, k=k))
def mutate(s, sub=0.04, indel=0.004):
    out = []
    i = 0
    while i < len(s):
        r = rng.random()
        if r < indel / 2:                      # deletion
            i += rng.randint(1, 5); continue
        if r < indel:                          # insertion
            out.append(draw(rng.randint(1, 5)))
        c = s[i]
        if rng.random() < sub:
            c = draw()
        out.append(c); i += 1
    return "".join(out)
seqs = [draw(L)]
while len(seqs) < n:
    seqs.append(mutate(rng.choice(seqs)))
seqs = [mutate(s, 0.02, 0.002) for s in seqs]
for k, s in enumerate(seqs):
    print(f">p{k}\n{s}")
