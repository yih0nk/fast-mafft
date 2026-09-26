#!/bin/bash
# Run a command while holding an exclusive timing lock (mkdir is atomic; macOS has no flock).
LOCK="$(cd "$(dirname "$0")" && pwd)/.timing.lock"
until mkdir "$LOCK" 2>/dev/null; do sleep 3; done
trap 'rmdir "$LOCK"' EXIT
"$@"
