#!/bin/bash
# Download the exact benchmark inputs used in RESULTS.md into bench/data/ and verify their checksums.
#
#   bash bench/get_data.sh              # the 32 FASTA inputs (93 MB download, ~325 MB unpacked)
#   bash bench/get_data.sh --long       # also the 12 long-genome inputs (real SARS-CoV-2 collections and
#                                       # simulated 30-50 kb sets; 50 MB download, ~210 MB) into bench/data/long/
#   bash bench/get_data.sh --balibase   # also BAliBASE 3.0 (386 *.tfa cases) into bench/data/balibase/
#
# On CARC run this on a login node (compute nodes may not have internet access).
# Provenance and licenses of every file: bench/DATA.md. Checksums: bench/DATA_MD5.txt.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
URL=https://github.com/yih0nk/fast-mafft/releases/download/data-v1/fast-mafft-bench-data-v1.tar.gz
LONG_URL=https://github.com/yih0nk/fast-mafft/releases/download/data-long-v1/fast-mafft-bench-long-v1.tar.gz
BALIBASE_URL=https://bigest-icube.fr/balibase/BalibaseDownload/BAliBASE_R1-5.tar.gz
cd "$HERE"

if [ -f data/rnasim50k.fa ] || [ -f data/large/rnasim50k.fa ]; then
  echo "bench/data already present; verifying only"
else
  tmp=$(mktemp -d)
  echo "downloading $URL"
  curl -fL --retry 3 -# -o "$tmp/data.tar.gz" "$URL"
  tar -xzf "$tmp/data.tar.gz" -C "$HERE"
  rm -rf "$tmp"
fi

LONG=0; BALI=0
for a in "$@"; do
  case $a in --long) LONG=1;; --balibase) BALI=1;; *) echo "unknown option $a"; exit 2;; esac
done
if [ $LONG = 1 ]; then
  if [ -f data/long/sarscov2_2000.fa ]; then
    echo "bench/data/long already present; verifying only"
  else
    tmp=$(mktemp -d)
    echo "downloading $LONG_URL"
    curl -fL --retry 3 -# -o "$tmp/long.tar.gz" "$LONG_URL"
    tar -xzf "$tmp/long.tar.gz" -C "$HERE"
    rm -rf "$tmp"
  fi
fi

python3 - "$HERE" $LONG <<'EOF'
import hashlib, os, sys
here = sys.argv[1]
bad = 0; n = 0
lists = ["DATA_MD5.txt"] + (["DATA_MD5_LONG.txt"] if sys.argv[2] == "1" else [])
for line in (l for f in lists for l in open(os.path.join(here, f))):
    want, path = line.split()
    p = os.path.join(here, path)
    n += 1
    if not os.path.exists(p):
        print("MISSING", path); bad += 1; continue
    h = hashlib.md5()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    if h.hexdigest() != want:
        print("BAD MD5", path); bad += 1
print("%d files checked, %d problems" % (n, bad))
sys.exit(1 if bad else 0)
EOF

if [ $BALI = 1 ]; then
  if [ -d data/balibase/bb3_release ]; then
    echo "BAliBASE already present"
  else
    tmp=$(mktemp -d)
    echo "downloading BAliBASE 3.0 from its official site"
    curl -fL --retry 3 -# -o "$tmp/bb.tar.gz" "$BALIBASE_URL"
    mkdir -p data/balibase && tar -xzf "$tmp/bb.tar.gz" -C data/balibase
    rm -rf "$tmp"
  fi
  echo "BAliBASE cases: $(find data/balibase -name '*.tfa' | wc -l)"
fi
