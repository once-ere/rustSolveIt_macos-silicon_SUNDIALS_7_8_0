# run_integrator_matrix.ps1 — run all 63 integrator configurations in both
# the MSVC-built C REBOUND and the pure-Rust port, and compare the raw-bit
# dumps. Prints one line per configuration and a final tally.
#
# Usage:   .\run_integrator_matrix.ps1  [nsteps]
# Default nsteps = 500, which is what section 15.1 of rebound_rust.md records.
#
# NOTE: every harness in this folder writes its dump to the same two
# filenames, state_c_final.txt and state_rust_final.txt. This script
# consumes them (renaming each pair to matrix_c.txt / matrix_rust.txt
# before comparing), so after it runs, the shearing-sheet dumps are gone.
# Re-create them by re-running that pair:
#     .\problem_test.exe
#     ..\target\release\examples\shearing_sheet_test.exe
#
# Part of the rebound_rs verification suite, GPL-3.0-or-later.

param([int]$nsteps = 500)

# Deliberately NOT "Stop": some integrator configurations legitimately print
# warnings to the error stream (for example bs-maxdt reports a clamped step),
# and PowerShell would otherwise treat a native program's stderr output as a
# terminating error and abandon the sweep partway through.
$ErrorActionPreference = "Continue"
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$rust = Join-Path $here "..\target\release\examples\integrators_test.exe"
$cexe = Join-Path $here "integrators_test.exe"

foreach ($p in @($rust, $cexe)) {
    if (-not (Test-Path $p)) { Write-Host "missing executable: $p" -ForegroundColor Red; exit 1 }
}

# The 63 configurations. `order` is only read by leapfrog; the others encode
# their settings in the name (see integrators_test.c for the full mapping).
$configs = @()
$configs += ,@("none", 2)
$configs += ,@("ias15", 2)
foreach ($o in 2, 4, 6, 8) { $configs += ,@("leapfrog", $o) }
foreach ($n in "whfast", "whfast-c11", "whfast-c17", "whfast-dh", "whfast-whds",
                "whfast-bary", "whfast-mk", "whfast-comp", "whfast-lazy",
                "whfast-usafe") { $configs += ,@($n, 2) }
foreach ($n in "saba", "saba-1", "saba-2", "saba-3", "saba-4", "saba-cm2",
                "saba-cl2", "saba-104", "saba-864", "saba-h844", "saba-h864",
                "saba-h1064", "saba-usafe") { $configs += ,@($n, 2) }
foreach ($n in "janus", "janus-2", "janus-4", "janus-8", "janus-10") { $configs += ,@($n, 2) }
$configs += ,@("eos", 2)
foreach ($d in 0..8) { $configs += ,@("eos-$d-$d", 2) }
foreach ($n in "eos-2-7", "eos-5-8", "eos-usafe") { $configs += ,@($n, 2) }
foreach ($n in "mercurius", "mercurius-usafe", "mercurius-c4", "mercurius-c5",
                "mercurius-inf", "mercurius-hill01") { $configs += ,@($n, 2) }
foreach ($n in "bs", "bs-tight", "bs-loose", "bs-maxdt") { $configs += ,@($n, 2) }
foreach ($n in "trace", "trace-pbs", "trace-ias15", "trace-hill1",
                "trace-perinone", "trace-eta001") { $configs += ,@($n, 2) }

Write-Host ("Running {0} configurations at {1} steps each." -f $configs.Count, $nsteps)

$identical = 0
$failed = @()

foreach ($c in $configs) {
    $name = $c[0]; $order = $c[1]
    $label = if ($name -eq "leapfrog") { "leapfrog(order $order)" } else { $name }

    Push-Location $here
    & $cexe $name $order $nsteps 2>$null | Out-Null
    if (-not (Test-Path "state_c_final.txt")) { Pop-Location; $failed += "$label (C produced no output)"; continue }
    Move-Item -Force "state_c_final.txt" "matrix_c.txt"

    & $rust $name $order $nsteps 2>$null | Out-Null
    if (-not (Test-Path "state_rust_final.txt")) { Pop-Location; $failed += "$label (Rust produced no output)"; continue }
    Move-Item -Force "state_rust_final.txt" "matrix_rust.txt"
    Pop-Location

    $diff = Compare-Object (Get-Content (Join-Path $here "matrix_c.txt")) `
                           (Get-Content (Join-Path $here "matrix_rust.txt"))
    if ($diff) {
        $failed += $label
        Write-Host ("  {0,-24} MISMATCH ({1} differing lines)" -f $label, $diff.Count) -ForegroundColor Red
    } else {
        $identical++
        Write-Host ("  {0,-24} identical" -f $label)
    }
}

Write-Host ""
Write-Host ("{0} of {1} configurations bit-identical." -f $identical, $configs.Count)
if ($failed.Count -gt 0) {
    Write-Host "FAILURES:" -ForegroundColor Red
    $failed | ForEach-Object { Write-Host "  $_" -ForegroundColor Red }
    exit 1
}
Write-Host "ALL CONFIGURATIONS BIT-IDENTICAL" -ForegroundColor Green
