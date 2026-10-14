#!/usr/bin/env python3
# Python 3.6-compatible (CARC). Usage: carc_run.py MAFFT_BIN THREADS INPUT LABEL
# Runs `MAFFT_BIN --thread THREADS INPUT`, prints one result line:
#   LABEL dataset t=N wall=S maxrss_gb=G md5=XXXXXXXXXXXX  stages=6mer:a,prog1:b,distmsa:c,prog2:d
# Stage times come from MAFFT's stderr progress markers; peak RSS is the largest child (disttbfast).
import subprocess, sys, time, hashlib, resource, os, re, threading
mafft, threads, inp, label = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]
t0 = time.time()
p = subprocess.Popen([mafft, "--thread", threads, inp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
h = hashlib.md5()
def pump():
    for chunk in iter(lambda: p.stdout.read(1 << 20), b""):
        h.update(chunk)
th = threading.Thread(target=pump); th.start()
marks, buf = {}, b""
keys = [("6mer", b"Making a distance matrix .."), ("prog1", b"Progressive alignment 1/2"),
        ("distmsa", b"Making a distance matrix from msa"), ("prog2", b"Progressive alignment 2/2"),
        ("end", b"Strategy:")]
while True:
    ch = p.stderr.read(4096)
    if not ch:
        break
    buf += ch
    now = time.time() - t0
    for k, pat in keys:
        if k not in marks and pat in buf:
            marks[k] = now
    if len(buf) > 65536:
        buf = buf[-4096:]
p.wait(); th.join()
wall = time.time() - t0
rss_gb = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024.0 / 1024.0  # Linux: KiB
order = ["6mer", "prog1", "distmsa", "prog2", "end"]
st = []
for a, b in zip(order, order[1:]):
    if a in marks and b in marks:
        st.append("%s:%.1f" % (a, marks[b] - marks[a]))
print("%s %s t=%s wall=%.1f maxrss_gb=%.2f md5=%s exit=%d stages=%s" % (
    label, os.path.basename(inp).replace(".fa", ""), threads, wall, rss_gb, h.hexdigest()[:12],
    p.returncode, ",".join(st)), flush=True)
