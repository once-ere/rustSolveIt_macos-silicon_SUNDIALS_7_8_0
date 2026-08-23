#!/usr/bin/env bash
# macos_verify_physics.sh — the macOS/Apple-Silicon byte-identity gate for
# the physics.
#
# Re-runs, on macOS (arm64), exactly the three suites whose Linux outputs
# are recorded in evidence/port-7.8.0/, in the same concatenated-log
# formats, and diffs the results:
#
#   examples-7.8.0.log            6 self-checking physics examples
#   collision-scripts-7.8.0.log   12 collision scripts
#   dynamic-notebooks-7.8.0.log   59 dynamic notebooks
#
# Normalisations (the same ones the originals carry):
#   * the OS-assigned scene port  ->  <port>
#   * the current directory       ->  <cwd>
# Line endings are compared with --strip-trailing-cr so the gate is
# indifferent to how the evidence files were checked out.
#
# Run from the repository root:
#   bash tools/macos_verify_physics.sh
#
# Writes evidence/macos/{examples,collision-scripts,dynamic-notebooks}-macos.log
# and prints IDENTICAL or DIFFERS per suite. Exit 0 only if every suite is
# identical or matches a pinned, documented divergence diff byte-for-byte.
set -u
cd "$(dirname "$0")/.."
ROOT="$PWD"
BIN="$ROOT/target/release/posim"
OUT="$ROOT/evidence/macos"
REF="$ROOT/evidence/port-7.8.0"
mkdir -p "$OUT"
export POSIM_NO_BROWSER=1

[ -x "$BIN" ] || { echo "build first: cargo build --release -p posim" >&2; exit 2; }

normalize () {
    # scene ports vary per run; cwd is machine-specific
    sed -e 's#127\.0\.0\.1:[0-9][0-9]*#127.0.0.1:<port>#g' \
        -e "s#$(printf '%s' "$ROOT" | sed 's/[.[\*^$/]/\\&/g')#<cwd>#g"
}

fail=0

# The reference logs were concatenated on a glibc host whose locale
# collation ordered the section globs; rather than imitate that locale,
# take the section order straight out of each reference log.
sections () { grep '^#####* ' "$1" | tr -d '\r' | sed 's/^#* //'; }

# ---- 1. the six examples --------------------------------------------------
: > "$OUT/examples-macos.log"
sections "$REF/examples-7.8.0.log" | while read -r ex; do
    echo "############ $ex" >> "$OUT/examples-macos.log"
    "$ROOT/target/release/examples/$ex" >> "$OUT/examples-macos.log" 2>&1
    echo "exit=$?" >> "$OUT/examples-macos.log"
done

# ---- 2. the twelve collision scripts -------------------------------------
: > "$OUT/collision-scripts-macos.log"
for f in scripts/collisions/*.posim; do
    echo "##### $f" >> "$OUT/collision-scripts-macos.log"
    "$BIN" --script "$f" >> "$OUT/collision-scripts-macos.log" 2>&1 \
        || { echo "script $f FAILED" >&2; fail=1; }
done

# ---- 3. the fifty-nine dynamic notebooks ----------------------------------
{ echo "<cwd>"; } > "$OUT/dynamic-notebooks-macos.log"
sections "$REF/dynamic-notebooks-7.8.0.log" | while read -r f; do
    echo "##### $f" >> "$OUT/dynamic-notebooks-macos.log"
    "$BIN" --script "$f" 2>&1 | normalize >> "$OUT/dynamic-notebooks-macos.log"
    echo "rc=${PIPESTATUS[0]}" >> "$OUT/dynamic-notebooks-macos.log"
done

# ---- compare --------------------------------------------------------------
# A divergence from the Linux evidence is ACCEPTED only when it is pinned
# byte-for-byte in evidence/macos/accepted-divergences-*.diff, with the
# cause documented in VERIFICATION_MACOS.md. Anything beyond a pinned diff
# is a regression and fails the gate.
for pair in "examples-macos.log:examples-7.8.0.log:accepted-divergences-examples.diff" \
            "collision-scripts-macos.log:collision-scripts-7.8.0.log:accepted-divergences-collisions.diff" \
            "dynamic-notebooks-macos.log:dynamic-notebooks-7.8.0.log:accepted-divergences-dynamic.diff"; do
    IFS=: read -r got want accepted <<< "$pair"
    d="$(diff --strip-trailing-cr "$REF/$want" "$OUT/$got")"
    if [ -z "$d" ]; then
        echo "IDENTICAL  $got == $want"
    elif [ -n "$accepted" ] && [ -f "$OUT/$accepted" ] \
         && [ "$d" = "$(cat "$OUT/$accepted" | tr -d '\r')" ]; then
        echo "IDENTICAL-MODULO-DOCUMENTED  $got vs $want (pinned: $accepted)"
    else
        echo "DIFFERS    $got vs $want (beyond the pinned divergences)"
        printf '%s\n' "$d" | head -20
        fail=1
    fi
done
exit $fail
