# Simulate tree-structured DNA sequences (substitutions + short indels) for MAFFT timing.
import random, sys
n, L, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1
rng = random.Random(seed)
B = "ACGT"
def mutate(s, sub=0.02, indel=0.003):
    out = []
    i = 0
    while i < len(s):
        r = rng.random()
        if r < indel / 2:                      # deletion
            i += rng.randint(1, 5); continue
        if r < indel:                          # insertion
            out.append("".join(rng.choice(B) for _ in range(rng.randint(1, 5))))
        c = s[i]
        if rng.random() < sub:
            c = rng.choice(B)
        out.append(c); i += 1
    return "".join(out)
seqs = ["".join(rng.choice(B) for _ in range(L))]
while len(seqs) < n:
    seqs.append(mutate(rng.choice(seqs)))
seqs = [mutate(s, 0.01, 0.002) for s in seqs]
for k, s in enumerate(seqs):
    print(f">s{k}\n{s}")
