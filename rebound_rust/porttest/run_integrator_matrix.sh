#!/usr/bin/env bash
# run_integrator_matrix.sh — run all 63 integrator configurations in both
# the clang-built C REBOUND and the pure-Rust port, and compare the raw-bit
# dumps. Prints one line per configuration and a final tally.
#
# This is the macOS/POSIX twin of run_integrator_matrix.ps1 (the Windows
# version). Usage:
#     bash run_integrator_matrix.sh [nsteps]
# Default nsteps = 500, which is what the master document records.
#
# Before running it you need both executables built:
#     cd .. && cargo build --release --example integrators_test
#     clang -I"../../rebound/rebound/src" -D_GNU_SOURCE -O2 -ffp-contract=off \
#           integrators_test.c ../../rebound/rebound/src/librebound_static.a \
#           -lm -o integrators_test
#
# NOTE: every harness in this folder writes its dump to the same two
# filenames, state_c_final.txt and state_rust_final.txt. This script
# consumes them (renaming each pair to matrix_c.txt / matrix_rust.txt
# before comparing), so after it runs, the shearing-sheet dumps are gone.
# Re-create them by re-running that pair:
#     ./problem_test 400
#     ../target/release/examples/shearing_sheet_test 400
#
# Part of the rebound_rs verification suite, GPL-3.0-or-later.
set -u
nsteps="${1:-500}"
here="$(cd "$(dirname "$0")" && pwd)"
rust="$here/../target/release/examples/integrators_test"
cexe="$here/integrators_test"

for p in "$rust" "$cexe"; do
    [ -x "$p" ] || { echo "missing executable: $p" >&2; exit 1; }
done

# The 63 configurations. The second field (the leapfrog order) is only
# read by leapfrog; the others encode their settings in the name (see
# integrators_test.c for the full mapping).
configs=(
    "none 2" "ias15 2"
    "leapfrog 2" "leapfrog 4" "leapfrog 6" "leapfrog 8"
    "whfast 2" "whfast-c11 2" "whfast-c17 2" "whfast-dh 2" "whfast-whds 2"
    "whfast-bary 2" "whfast-mk 2" "whfast-comp 2" "whfast-lazy 2"
    "whfast-usafe 2"
    "saba 2" "saba-1 2" "saba-2 2" "saba-3 2" "saba-4 2" "saba-cm2 2"
    "saba-cl2 2" "saba-104 2" "saba-864 2" "saba-h844 2" "saba-h864 2"
    "saba-h1064 2" "saba-usafe 2"
    "janus 2" "janus-2 2" "janus-4 2" "janus-8 2" "janus-10 2"
    "eos 2"
    "eos-0-0 2" "eos-1-1 2" "eos-2-2 2" "eos-3-3 2" "eos-4-4 2" "eos-5-5 2"
    "eos-6-6 2" "eos-7-7 2" "eos-8-8 2"
    "eos-2-7 2" "eos-5-8 2" "eos-usafe 2"
    "mercurius 2" "mercurius-usafe 2" "mercurius-c4 2" "mercurius-c5 2"
    "mercurius-inf 2" "mercurius-hill01 2"
    "bs 2" "bs-tight 2" "bs-loose 2" "bs-maxdt 2"
    "trace 2" "trace-pbs 2" "trace-ias15 2" "trace-hill1 2"
    "trace-perinone 2" "trace-eta001 2"
)

echo "Running ${#configs[@]} configurations at $nsteps steps each."

identical=0
failed=()

cd "$here"
for c in "${configs[@]}"; do
    name="${c% *}"; order="${c#* }"
    if [ "$name" = "leapfrog" ]; then label="leapfrog(order $order)"; else label="$name"; fi

    rm -f state_c_final.txt state_rust_final.txt
    "$cexe" "$name" "$order" "$nsteps" >/dev/null 2>&1
    if [ ! -f state_c_final.txt ]; then failed+=("$label (C produced no output)"); continue; fi
    mv -f state_c_final.txt matrix_c.txt

    "$rust" "$name" "$order" "$nsteps" >/dev/null 2>&1
    if [ ! -f state_rust_final.txt ]; then failed+=("$label (Rust produced no output)"); continue; fi
    mv -f state_rust_final.txt matrix_rust.txt

    if cmp -s matrix_c.txt matrix_rust.txt; then
        identical=$((identical+1))
        printf '  %-24s identical\n' "$label"
    else
        n=$(diff matrix_c.txt matrix_rust.txt | grep -c '^[<>]')
        failed+=("$label")
        printf '  %-24s MISMATCH (%s differing lines)\n' "$label" "$n"
    fi
done

echo ""
echo "$identical of ${#configs[@]} configurations bit-identical."
if [ "${#failed[@]}" -gt 0 ]; then
    echo "FAILURES:"
    for f in "${failed[@]}"; do echo "  $f"; done
    exit 1
fi
echo "ALL CONFIGURATIONS BIT-IDENTICAL"
