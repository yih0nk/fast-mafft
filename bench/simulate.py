# Simulate tree-structured DNA sequences (substitutions + short indels) for MAFFT timing.
#   python3 simulate.py n L seed [sub indel final_sub final_indel]
# The rates are per generation (default 0.02 0.003) and for the final pass over every sequence
# (default 0.01 0.002); with the defaults the output is the same as before the rates were arguments.
import random, sys
n, L, seed = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]) if len(sys.argv) > 3 else 1
SUB, INDEL, FSUB, FINDEL = [float(x) for x in sys.argv[4:8]] if len(sys.argv) > 7 else (0.02, 0.003, 0.01, 0.002)
rng = random.Random(seed)
B = "ACGT"
def mutate(s, sub=SUB, indel=INDEL):
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
seqs = [mutate(s, FSUB, FINDEL) for s in seqs]
for k, s in enumerate(seqs):
    print(f">s{k}\n{s}")
