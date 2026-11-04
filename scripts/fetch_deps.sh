#!/usr/bin/env sh
# Fetches the pinned dependencies into deps/ and builds param_report.
#   t-evolve-cryptoalg     release v1.1 (commit 1395c3f; protocol code = release v1.0, 54a4e7e):
#                          logged parameter chain used as the baseline
#   t-evolve-national-scale  commit d0885a3: param_report used by the tool (identical to v1.0 for
#                          E <= 2183 blocks, i.e. up to about 1.02 * 10^6 voters, where the tail
#                          exponent stays 187 ln 2 and two sampler levels suffice; adds the tail-exponent rule for
#                          larger trees and a third sampler level above about 2^27.5)
#   lattice-estimator      commit 53da598: hardness estimates (needs SageMath or passagemath)
set -eu
cd "$(dirname "$0")/.."
mkdir -p deps

fetch() {   # fetch <dir> <url> <ref> <expected short hash>
    if [ ! -d "deps/$1" ]; then git clone -q "$2" "deps/$1"; fi
    git -C "deps/$1" fetch -q --tags origin
    git -C "deps/$1" checkout -q "$3"
    test "$(git -C "deps/$1" rev-parse --short=7 HEAD)" = "$4"
    echo "$1: $4"
}

fetch t-evolve-cryptoalg https://github.com/futures1mple/t-evolve-cryptoalg.git v1.1 1395c3f
fetch t-evolve-national-scale https://github.com/futures1mple/t-evolve-national-scale.git d0885a3 d0885a3
fetch lattice-estimator https://github.com/malb/lattice-estimator.git 53da598 53da598

make -s -C deps/t-evolve-national-scale build/param_report
echo "param_report: deps/t-evolve-national-scale/build/param_report"
