#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

mkdir -p build

clang++ -std=c++17 -O2 engine.cpp -o build/pinball \
  $(pkg-config --cflags --libs sdl2 SDL2_image SDL2_mixer)

echo "Built: build/pinball"

