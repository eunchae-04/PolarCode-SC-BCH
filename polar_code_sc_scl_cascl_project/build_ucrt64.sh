#!/usr/bin/env bash
set -euo pipefail

gcc -std=c11 -O2 -Wall -Wextra -Wpedantic \
    -o compare.exe polar_code_sc_scl_cascl_compare.c -lm

echo
echo "Build successful: compare.exe"
