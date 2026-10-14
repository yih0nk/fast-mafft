#!/bin/bash
# Download the exact benchmark inputs used in RESULTS.md into bench/data/ and verify their checksums.
#
#   bash bench/get_data.sh              # the 32 FASTA inputs (93 MB download, ~325 MB unpacked)
#   bash bench/get_data.sh --balibase   # also BAliBASE 3.0 (386 *.tfa cases) into bench/data/balibase/
#
# On CARC run this on a login node (compute nodes may not have internet access).
# Provenance and licenses of every file: bench/DATA.md. Checksums: bench/DATA_MD5.txt.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
URL=https://github.com/yih0nk/fast-mafft/releases/download/data-v1/fast-mafft-bench-data-v1.tar.gz
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

python3 - "$HERE" <<'EOF'
import hashlib, os, sys
here = sys.argv[1]
bad = 0; n = 0
for line in open(os.path.join(here, "DATA_MD5.txt")):
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

if [ "${1:-}" = "--balibase" ]; then
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
