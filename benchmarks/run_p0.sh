#!/usr/bin/env bash
set -euo pipefail

godot_bin="${GODOT_BIN:-godot}"
output_path="${1:-p0-baseline.json}"

"$godot_bin" --headless --path . --script benchmarks/p0_baseline.gd -- \
  --warmup=10 \
  --ticks=30 \
  --populations=1,16,64,256,1000 \
  "--output=$output_path"
