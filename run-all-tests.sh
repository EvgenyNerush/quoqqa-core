#!/usr/bin/env bash
set -u

cargo build --bins || exit 1

for src in src/bin/*.rs; do
    [[ -e "$src" ]] || { echo "No binaries found in src/bin"; exit 1; }
    name=$(basename "$src" .rs)
    echo -e "\033[1m=== Running $name ===\033[0m"
    cargo run --bin "$name"
done

