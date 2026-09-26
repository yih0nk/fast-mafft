# Usage: MAFFT=/path/to/bin/mafft python3 profile.py IN.fa THREADS
# Run MAFFT (default FFT-NS-2) and time each stage from its stderr progress messages.
import subprocess, sys, time, re, hashlib, os
inp, threads = sys.argv[1], int(sys.argv[2])
t0 = time.time()
p = subprocess.Popen([os.environ.get("MAFFT", "mafft"), "--thread", str(threads), inp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
import threading
out_chunks = []
threading.Thread(target=lambda: out_chunks.append(p.stdout.read()), daemon=True).start()
events, buf = [], b""
while True:
    ch = p.stderr.read(1)
    if not ch: break
    if ch in (b"\n", b"\r"):
        line = buf.decode(errors="replace").strip(); buf = b""
        if line and not line.startswith("STEP") and not re.match(r"^\s*\d+\s*/\s*\d+", line):
            events.append((time.time() - t0, line))
    else:
        buf += ch
p.wait(); total = time.time() - t0
time.sleep(0.2)
out = b"".join(out_chunks)
print(f"threads={threads} total={total:.1f}s md5={hashlib.md5(out).hexdigest()[:12]}")
keys = ("Making a distance matrix", "Constructing a UPGMA", "Progressive alignment", "Strategy", "done")
prev = 0.0
for t, line in events:
    if any(k in line for k in keys):
        print(f"  {t:8.1f}s  (+{t-prev:6.1f})  {line[:70]}"); prev = t
os.makedirs("out", exist_ok=True)
open(f"out/{os.path.basename(inp)}.t{threads}.aln", "wb").write(out)
