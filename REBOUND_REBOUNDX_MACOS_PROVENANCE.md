# REBOUND_REBOUNDX_MACOS_PROVENANCE — the N-body ports, verified on this Mac

This page records how the two pure-Rust astronomy libraries in this
repository — `rebound_rust/` and `reboundx_rust/` — got here, and every
command used to build, test, verify and run them on **macOS / Apple
Silicon**, with the measured result of each. It is self-contained: every
command you need appears on this page, run from the repository root.

## 1. What these two folders are

- **`rebound_rust/`** (crate `rebound_rs` 5.1.1) — a pure-Rust, line-for-line
  translation of [REBOUND](https://github.com/hannorein/rebound), the
  open-source N-body gravity simulator by Hanno Rein and collaborators
  (version 5.1.1, commit `dad5f97806ecbb408dcaff728851c64e67f9f6eb`).
  29 source modules, 13 example programs, 394 automated tests, and one
  Jupyter notebook per example in `rebound_rust/notebooks/`.
- **`reboundx_rust/`** (crate `reboundx_rs` 5.1.0) — a pure-Rust translation
  of [REBOUNDx](https://github.com/dtamayo/reboundx) by Dan Tamayo, Hanno
  Rein and collaborators, which adds extra physics (general relativity,
  tides, spin, migration, radiation forces…) to a REBOUND simulation.
  33 source modules, 4 example programs, 137 automated tests; its
  notebooks live beside the REBOUND ones.

Both are **zero `unsafe`, zero external dependencies (std only), zero
build warnings**, keep the exact C function names, and are
**GPL-3.0-or-later** (each folder carries its own LICENSE — note this
differs from the BSD-3-Clause license of the physics simulator in the
rest of this repository; the folders are independent works that live
together here).

They were first produced and verified bit-for-bit on Windows 11, as the
physics acceptance test of the rustSolveIt Windows port. This repository
carries them **verified the same way on macOS / Apple Silicon**: the same
experiments, re-run on this machine against the original C compiled with
Apple clang on this machine.

The master document — written for a reader who has never programmed, with
the complete story and every command — is
`rebound_rust/rebound_rust.md`, also typeset as
`rebound_rust/rebound_rust.pdf` (36 pages). You do not need it to
reproduce the results below; everything needed is on this page.

## 2. Machine and toolchain (as measured, 2026-08-27)

| item | value |
|---|---|
| machine | Apple **M5 Max**, 128 GB RAM, arm64 |
| OS | macOS Tahoe 26 (26.6.1, build 25G76) |
| C compiler | Apple clang 21.0.0 (Xcode Command Line Tools) |
| Rust | rustc / cargo 1.94.0, `aarch64-apple-darwin` |
| Python (notebooks only) | 3.12 via `uv`-managed venv |

## 3. What changed in the port to macOS

No Rust code changed, with one commented exception class: **doc-comment
text only** (rustc 1.94's documentation checker rejects bare
`<placeholder>` and `particles[0]`-style text under the crates'
`deny(warnings)`; those spans are now backtick-quoted in 2 REBOUND and 5
REBOUNDx modules — no executable line moved, and both test suites were
re-run green afterwards). Everything else:

1. `rebound_rust/porttest/macos_shim/rand_r_glibc.c` — **new**. Upstream
   REBOUND vendors glibc's `rand_r` random generator only under
   `#ifdef _WIN32`, so on macOS the C reference silently used Apple's
   different `rand_r` and every random initial condition diverged (the
   shearing sheet placed 1,441 particles instead of 1,482). The shim
   supplies the same glibc algorithm at link time; upstream source is
   untouched. The Rust port implements the glibc algorithm everywhere
   and needs no shim.
2. `rebound_rust/porttest/run_integrator_matrix.sh` — **new**: the
   macOS/POSIX twin of the Windows PowerShell sweep over all 63
   integrator configurations.
3. `rebound_rust/notebooks/make_notebooks.py` — no longer hard-codes the
   Windows `.exe` suffix (five sites), so the generated notebooks run on
   both platforms; `tides_spin_pseudo.ipynb` is now generated (the
   Windows folder had defined but never carried it) — **every one of the
   17 examples has exactly one companion notebook**.
4. Documentation — `rebound_rust.md`/`.tex`/`.pdf`, both READMEs and the
   three acceptance-test write-ups re-platformed to macOS with this
   machine's measured results; the Windows stories (the MSVC `pow`
   difference, the MSVC variable-length-array shim) are kept and
   labelled as history.
5. The upstream C trees are **not** committed: `.gitignore` reserves
   `/rebound/` and `/reboundx/` for local clones used as the
   verification reference (commands below).

## 4. Build and test (commands + measured results)

```bash
cd rebound_rust
cargo build --release --all-targets
cargo test  --release
cargo clippy --release --all-targets
cargo doc   --no-deps
cd ../reboundx_rust
cargo build --release --all-targets
cargo test  --release
cargo clippy --release --all-targets
cargo doc   --no-deps
```

Measured: rebound_rs **394 passed, 0 failed**; reboundx_rs **137 passed,
0 failed**; **zero** compiler, clippy and rustdoc warnings in all of it.

## 5. The C references, built on this machine

```bash
git clone https://github.com/hannorein/rebound.git rebound/rebound
(cd rebound/rebound && git checkout dad5f97806ecbb408dcaff728851c64e67f9f6eb)
git clone https://github.com/dtamayo/reboundx.git reboundx

cd rebound/rebound/src
clang -c -DBUILDINGLIBREBOUND -D_GNU_SOURCE -DSERVER \
      -DGITHASH=dad5f97806ecbb408dcaff728851c64e67f9f6eb \
      -O2 -ffp-contract=off *.c
ar rcs librebound_static.a *.o
cd ../../../reboundx/src
clang -c -I../../rebound/rebound/src -I. -D_GNU_SOURCE -DLIBREBOUNDX \
      -O2 -ffp-contract=off *.c
ar rcs libreboundx.a *.o
cd ../..
```

`-ffp-contract=off` stops clang fusing multiply-adds into FMA
instructions, which would change rounding relative to the Rust build —
the clang equivalent of the MSVC `/fp:precise` the Windows reference
used. All 64 C files compile; the MSVC variable-length-array shim is not
needed under clang.

Build the sixteen comparison harnesses:

```bash
cd rebound_rust/porttest
clang -c -O2 -ffp-contract=off macos_shim/rand_r_glibc.c -o macos_shim/rand_r_glibc.o
for f in addfmt_test archive_test bs_pow_diff derivatives_test \
         frequency_test integrators_test kepler_rectilinear_c libm_diff \
         movetocom_var_c movetocom_var_test problem_test; do
  clang -I../../rebound/rebound/src -D_GNU_SOURCE -O2 -ffp-contract=off \
        "$f.c" macos_shim/rand_r_glibc.o \
        ../../rebound/rebound/src/librebound_static.a -lm -o "$f"
done
cd ../../reboundx_rust/porttest
for f in tides_spin_pseudo_c tides_spin_kozai_c tides_spin_migration_c \
         rebx_binary_roundtrip_c rebx_binary_read_c; do
  clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
        -O2 -ffp-contract=off "$f.c" \
        ../../rebound_rust/porttest/macos_shim/rand_r_glibc.o \
        ../../reboundx/src/libreboundx.a \
        ../../rebound/rebound/src/librebound_static.a -lm -o "$f"
done
cd ../..
```

## 6. The verification record (every result measured here)

| check | command (from the repo root) | measured result |
|---|---|---|
| maths library: sin, cos, tan, atan2, pow, sqrt, fmod, exp, log × 200,000 samples each, + an exp/log differential | `cd rebound_rust/porttest && ./libm_diff && ../target/release/examples/libm_diff && cmp libm_c.txt libm_rust.txt` | **bit-identical — `pow` included** (the Windows port's one known `pow` difference does not exist here: both languages call Apple's libm) |
| BS step-size `pow` shapes | `./bs_pow_diff && ../target/release/examples/bs_pow_diff && cmp bs_pow_c.txt bs_pow_rust.txt` | **bit-identical** |
| integrator matrix, 63 configurations × 500 steps | `bash rebound_rust/porttest/run_integrator_matrix.sh 500` | **"63 of 63 configurations bit-identical. ALL CONFIGURATIONS BIT-IDENTICAL"** (run twice: before and after the shim) |
| ias15, 1000-step adaptive run | `./integrators_test ias15 2 1000` both sides, `cmp` | **bit-identical** |
| shearing sheet: seed 42, 1,482 particles, 400 steps, 102,478 collisions | `./problem_test 400 && ../target/release/examples/shearing_sheet_test 400 && shasum -a 256 state_c_final.txt state_rust_final.txt` | **byte-identical**, both SHA-256 `418c864dd1a610cbe8ea6d81ecafa1e4ce6d36837494177d9875ee820ef0766f` |
| 65 orbital-derivative functions | `./derivatives_test` + Rust twin, `cmp` | **130/130 outputs bit-identical** |
| frequency analysis (MFT, FMFT, FMFT2) | `./frequency_test` + Rust twin, `cmp` | **bit-identical** |
| `add_fmt` + built-in datasets | `./addfmt_test` + Rust twin, `cmp` | **bit-identical** |
| MEGNO / variational move-to-com | `./movetocom_var_c` + Rust twin, `cmp` | **bit-identical** |
| rectilinear Kepler probe | `./kepler_rectilinear_c` + Rust twin | **every hex bit-pattern identical** (one human-readable banner line formats its zero exponent differently — `e+00` vs `e0` — with no numbers involved) |
| Simulationarchive, whfast-usafe and ias15 | `./archive_test <cfg> write` / `continue` alternating languages, `cmp` | **bit-identical in all six directions** |
| web server | Rust `server_test` serves; `curl` fetches `/simulation` (3,448 bytes); C `./archive_test whfast load served.bin`, `cmp` | **bit-identical state** |
| REBOUNDx `tides_spin`: pseudo, kozai, migration — short and full runs | `cd reboundx_rust/porttest && ./tides_spin_<n>_c <t>` + Rust twins, `cmp` (t = 62.83185307179586 / 628.3185307179586; kozai 1000 and default 100,000) | **6/6 bit-identical** — the full kozai means both languages took the identical sequence of thousands of adaptive IAS15 steps |
| REBOUNDx binary files | `../target/release/examples/rebx_binary_roundtrip` and `./rebx_binary_roundtrip_c rebx_c_reference.bin`; `cmp -l a b \| wc -l`; read-backs via `./rebx_binary_read_c` | both files **10,784 bytes**; **26 differing bytes = the git-hash header stamp only**; Rust reads C's file **25/25 checks + full re-serialization**; C's dumps of the Rust file and its own file **identical (28 lines)** |
| the 17 Jupyter notebooks | section 7 below | **"All 17 notebooks executed with no errors"** |

## 7. The notebooks — every example has one

```bash
cd rebound_rust/notebooks
uv venv .venv
uv pip install -p .venv/bin/python nbformat nbclient ipykernel matplotlib
python3 make_notebooks.py            # regenerates all 17
.venv/bin/python run_notebooks.py    # executes all 17, writes outputs back
```

Measured: `17 notebooks written`, then
**`All 17 notebooks executed with no errors.`**

To read them in your browser:

```bash
python3 -m pip install --user jupyterlab
jupyter lab rebound_rust/notebooks/
```

## 8. Run any example directly

```bash
cd rebound_rust
cargo run --release --example integrators_test -- whfast 2 500
cargo run --release --example shearing_sheet_test -- 400
cd ../reboundx_rust
cargo run --release --example tides_spin_pseudo -- 62.83185307179586
```

(13 REBOUND examples, 4 REBOUNDx examples; each has the same name as its
notebook.)

## 9. Credit and licensing

REBOUND is © Hanno Rein and collaborators; REBOUNDx is © Dan Tamayo,
Hanno Rein and collaborators; both GPL-3.0-or-later, as are these
translations. **All of the science belongs to the original authors** —
cite **Rein & Liu 2012 (A&A 537, A128)** for any REBOUND use and
**Tamayo, Rein, Shi & Hernandez 2019 (MNRAS 491, 2885)** for any
REBOUNDx use, plus the per-feature papers listed in
`rebound_rust/README.md`. Nothing from these ports has been or will be
submitted to the upstream projects (their contribution policy asks that
no AI-generated code be submitted, and this port honours it); report
problems with the port here, not upstream.
