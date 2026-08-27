# The REBOUNDx Tests: Did the Rust Version Get the Same Answer?

**Short answer: yes — every single bit, in all six test runs, on this Mac.**

This document explains what REBOUNDx is, which three simulations we used to test
our Rust translation of it, exactly how we ran the test on macOS (Apple
Silicon), and what we found. It is written for someone who has never done this
before. Every command is printed in full, so you never need to open another
document.

> **A note on `~/work` in the commands below.** The `~` character is Terminal
> shorthand for your *home folder* (the one holding your Desktop, Documents and
> so on). `~/work` is a stand-in for whatever folder you put this project in —
> you do not have to create one with that name. If your copy lives in
> `~/astronomy`, then read every `~/work/...` below as `~/astronomy/...`. Only
> that leading part changes; the folder names after it are real and must match.

---

## Table of contents

1. [What is REBOUNDx?](#1-what-is-reboundx)
2. [What are "tides" and "spin", in plain language?](#2-what-are-tides-and-spin-in-plain-language)
3. [The three test simulations, and why these three](#3-the-three-test-simulations-and-why-these-three)
4. [What "bit-for-bit identical" means](#4-what-bit-for-bit-identical-means)
5. [The equipment](#5-the-equipment)
6. [Step 1 — Build the original C REBOUNDx](#6-step-1--build-the-original-c-reboundx)
7. [Step 2 — Build the Rust version](#7-step-2--build-the-rust-version)
8. [Step 3 — Run both and compare](#8-step-3--run-both-and-compare)
9. [The results](#9-the-results)
10. [Speed](#10-speed)
11. [What this proves, and what it does not](#11-what-this-proves-and-what-it-does-not)
12. [Files involved](#12-files-involved)

---

## 1. What is REBOUNDx?

**REBOUND** is a program that calculates how objects move under gravity —
planets around a star, moons around a planet, and so on. It was written in the
C programming language by Hanno Rein and collaborators.

**REBOUNDx** ("REBOUND eXtras") is a companion library by **Dan Tamayo**, Hanno
Rein and collaborators. It adds physics *beyond* plain gravity. You pick the
extra effects you want and switch them on:

| Effect | What it does |
|---|---|
| general relativity | Einstein's correction to Newton's gravity — it makes Mercury's orbit slowly rotate |
| tides | the stretching a planet feels from its star, which slowly circularises orbits |
| spin evolution | how a planet's rotation axis tips and drifts over time |
| radiation pressure | sunlight physically pushing on small dust grains |
| migration | planets spiralling inward or outward through a gas disk |

We translated REBOUNDx from C into **Rust**, and this document is the evidence
that the translation computes exactly the same numbers.

## 2. What are "tides" and "spin", in plain language?

**Tides.** The Moon pulls harder on the side of the Earth facing it than on the
far side. That difference stretches the Earth slightly — that is what makes the
ocean tides. The same thing happens to a planet orbiting close to its star, but
much more strongly. Because the planet is not perfectly rigid, the stretching
lags behind, and that lag steadily drains energy from the orbit.

**Spin.** Planets rotate. Tides push on that rotation too. Over time they tend
to:

- make the orbit **rounder** (circularisation),
- make the spin axis **stand up straight** relative to the orbit (obliquity
  damping),
- lock the rotation rate to the orbit (this is why the Moon always shows us the
  same face).

The `tides_spin` effect tracks all of this at once. It is the most complicated
effect in REBOUNDx, which is exactly why we chose it as the acceptance test.

## 3. The three test simulations, and why these three

These are three of REBOUNDx's own shipped examples. Together they exercise
every awkward part of the machinery.

### Test 1 — `tides_spin_pseudo_synchronization`

A "hot Jupiter": a giant planet orbiting extremely close to its star (its year
lasts about 3 of our days). It starts slightly elliptical, tilted 30°, and
spinning fast. Tides should rapidly round out the orbit, straighten the tilt,
and settle the spin.

**Exercises:** the WHFast integrator, the `tides_spin` force, and REBOUNDx's
ability to evolve spin vectors as extra differential equations alongside the
orbits.

### Test 2 — `tides_spin_kozai`

A planet with a distant stellar companion, undergoing a **Kozai cycle** — the
companion's pull periodically drives the planet's orbit to a very high
eccentricity, then back. During those high-eccentricity phases the planet
swings very close to its star.

**Exercises:** the **adaptive** IAS15 integrator, which chooses its own step
size. This is a much harder test than it sounds: for the results to match
bit-for-bit, the entire *sequence of chosen step sizes* has to match, which
means the extra forces must feed IAS15's internal error estimate identically.
It also runs a **second** effect at the same time (`gr_potential`, general
relativistic precession).

### Test 3 — `tides_spin_migration_driven_obliquity_tides`

Two Earth-sized planets migrating inward through a gas disk while tides evolve
their spins (the setup from Millholland & Laughlin 2019).

**Exercises:** two REBOUNDx forces at once (`tides_spin` +
`modify_orbits_forces`), three bodies instead of two, and — importantly —
**changing a parameter in the middle of the run**: migration is switched off
half-way through by setting the migration timescale to infinity.

### Coverage summary

| | Test 1 | Test 2 | Test 3 |
|---|---|---|---|
| Integrator | WHFast (fixed step) | IAS15 (**adaptive**) | WHFast (fixed step) |
| Bodies | 2 | 3 | 3 |
| REBOUNDx forces | 1 | 2 | 2 |
| Spin ODE evolution | yes | yes | yes |
| Rotation into the invariable plane | yes | yes | yes |
| Parameter changed mid-run | no | no | **yes** |

## 4. What "bit-for-bit identical" means

Computers store decimal numbers as 64 ones-and-zeros (**bits**), in a format
called IEEE-754 double precision.

**Bit-for-bit identical** means *all 64 bits match*, for every number — not
"agrees to 10 decimal places".

This strictness is necessary because these systems are **chaotic**: a difference
in the very last bit grows with every step until the two runs look completely
different. So bit-for-bit agreement is the only test that really proves the
arithmetic is the same.

Both programs therefore print every number as raw bits in hexadecimal. The
number 1.0 prints as `3ff0000000000000`. We then compare the files.

## 5. The equipment

| Thing | Version |
|---|---|
| Computer | Apple M5 Max (Apple Silicon, arm64), 128 GB RAM |
| Operating system | macOS Tahoe 26 (version 26.6.1, build 25G76) |
| C compiler | Apple `clang` 21.0.0 (clang-2100.1.1.101), from the Xcode Command Line Tools |
| Rust | `rustc` 1.94.0, `cargo` 1.94.0 (`aarch64-apple-darwin`) |
| REBOUND (C) | version 5.1.1, commit `dad5f978` |
| REBOUNDx (C) | version 5.1.0 |

No Linux, no virtual machine, no Homebrew compiler — everything native macOS:
Apple's clang, Apple's system maths library, Apple Silicon.

**What you need installed.** Two things, both free:

```bash
# 1. The Xcode Command Line Tools: Apple's C compiler (clang) and its
#    librarian (ar). A dialog pops up; click Install.
xcode-select --install

# 2. Rust, from the official installer at https://rustup.rs — it gives you
#    the rustc compiler and the cargo build tool for Apple Silicon.
```

## 6. Step 1 — Build the original C REBOUNDx

### 6a. Get both sources

```bash
cd ~/work
git clone https://github.com/hannorein/rebound.git rebound/rebound
git clone https://github.com/dtamayo/reboundx.git reboundx

# Pin REBOUND to the exact revision we tested against:
git -C rebound/rebound checkout dad5f97806ecbb408dcaff728851c64e67f9f6eb
```

(The REBOUNDx copy we tested was version 5.1.0.)

### 6b. Build C REBOUND, and make a static library from it

```bash
cd ~/work/rebound/rebound/src
clang -c -DBUILDINGLIBREBOUND -D_GNU_SOURCE -DSERVER \
      -DGITHASH=dad5f97806ecbb408dcaff728851c64e67f9f6eb \
      -O2 -ffp-contract=off *.c
ar rcs librebound_static.a *.o
```

What this does, piece by piece: `clang -c` compiles every C file (`*.c`) into
an **object file** (a `.o` file of machine code) without yet making a runnable
program; `-O2` turns on normal optimisation; the `-D...` flags set the same
build switches REBOUND's own build uses (including burning the source-code
revision into `GITHASH`); and `ar rcs` is the **librarian** — it bundles all 31
object files into one archive, `librebound_static.a`. No errors are expected.

> **Why `-ffp-contract=off`?** Modern Apple Silicon chips have an instruction
> that computes `a*b + c` in **one step** (a "fused multiply-add"), with one
> rounding instead of two. That is slightly *more* accurate — but it produces
> *different last bits* than doing the multiply and the add separately. Rust
> never fuses on its own, so we tell clang not to either; otherwise C and Rust
> could both be correct and still disagree in the last bit. This flag plays the
> same role that `/fp:precise` played for the Windows edition's MSVC compiler.

> **Why a *static* library?** The normal build makes a shared library (a
> `.dylib` on macOS), which only exposes the functions marked as public API.
> Our test harnesses call a few internal ones, so we bundle all the object
> files into a static library instead.

### 6c. Build C REBOUNDx

```bash
cd ~/work/reboundx/src
clang -c -I../../rebound/rebound/src -I. -D_GNU_SOURCE -DLIBREBOUNDX \
      -O2 -ffp-contract=off *.c
ar rcs libreboundx.a *.o
```

All **33** files compile, unmodified, with no errors — and `ar rcs` bundles
them into `libreboundx.a`. (`-I` tells the compiler an extra folder to search
for header files, here so REBOUNDx can find REBOUND's.)

> ### Windows history: on MSVC, two of the 33 files would not compile
>
> On the **Windows 11 edition** of this port, `gr_full.c` and
> `interpolation.c` failed to build, because they use a C99 feature called
> **variable-length arrays** — arrays whose size is decided while the program
> runs:
>
> ```c
> double a_const[N][3];   /* gr_full.c   */
> double u[n];            /* interpolation.c */
> ```
>
> Microsoft's C compiler has never supported these; GCC and Clang do. The
> Windows fix changed only *where the memory comes from* (the heap via
> `malloc` instead of the stack) and touched no arithmetic; the patched copies
> still live in `reboundx_rust/porttest/msvc_shim/` as part of that edition's
> record.
>
> **On macOS none of this is needed.** Apple's clang supports variable-length
> arrays, so both files compile exactly as their authors wrote them. And the
> Rust port never had the problem on any platform: Rust's `Vec` is a growable
> array and needs no workaround.

### 6d. One small macOS extra: the `rand_r` shim

REBOUND's source file `rebound.c` carries its own copy of a simple
random-number generator called `rand_r` (the algorithm from glibc, the GNU C
library) — but it only switches that copy on when building **for Windows**
(`#ifdef _WIN32`). On macOS that switch is off, so a plain C build silently
uses *Apple's* built-in `rand_r` instead — a **different generator**. Same
seed, different sequence of "random" numbers. Our Rust port implements the
glibc algorithm on every platform.

The fix keeps the upstream source untouched: the same glibc algorithm is
compiled once, as its own object file, and linked into **every** C test
harness. The linker then satisfies `rebound.c`'s call to `rand_r` from our
object file instead of from Apple's C library. Build it once:

```bash
cd ~/work/rebound_rust/porttest/macos_shim
clang -c -O2 -ffp-contract=off rand_r_glibc.c
```

That produces `rand_r_glibc.o`, which every harness link line below includes.

For honesty's sake: the three `tides_spin` tests in this document start from
**fixed, hand-written initial conditions** and draw no random numbers, so this
shim does not change their results — we link it everywhere so that every C
reference binary uses one and the same generator. Where it *does* matter is
the shearing-sheet test (1,482 randomly placed particles), whose comparison
died at particle 1 until the shim existed; that discovery is told in full in
the companion write-up `shearing_sheet_port_test.md`.

### 6e. Build the three C test harnesses

The shipped examples print rounded numbers, include `<unistd.h>` (a Unix
header nothing in them actually uses — and which does not exist on Windows,
where these harnesses were first written), and call `system("rm ...")`. Our
harnesses are those examples with exactly three changes, listed at the top of
each file:

1. `<unistd.h>` dropped (nothing from it is used),
2. `system("rm ...")` removed,
3. the text output replaced by a final dump of every state variable as raw bits.

**The physics setup is byte-for-byte the stock example.**

```bash
cd ~/work/reboundx_rust/porttest

clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off \
      tides_spin_pseudo_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.o \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a \
      -lm -o tides_spin_pseudo_c

clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off \
      tides_spin_kozai_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.o \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a \
      -lm -o tides_spin_kozai_c

clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off \
      tides_spin_migration_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.o \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a \
      -lm -o tides_spin_migration_c
```

> **The order of that link line matters.** The linker reads left to right,
> and when it meets a library it keeps only the pieces something *earlier*
> asked for. So: the harness comes first (it asks for REBOUNDx and REBOUND
> functions); then `rand_r_glibc.o` (an object file is always kept, which is
> what guarantees `rand_r` resolves to the glibc algorithm, not Apple's);
> then `libreboundx.a` **before** `librebound_static.a` (REBOUNDx calls into
> REBOUND, so the asker precedes the provider); and finally `-lm`, the system
> maths library.

> **One harmless warning is expected.** Compiling `tides_spin_pseudo_c.c`
> prints a single warning about a `void**` pointer cast. It is in the harness
> only — both libraries compile warning-free — and it does not affect any
> number.

## 7. Step 2 — Build the Rust version

```bash
cd ~/work/reboundx_rust
cargo build --release
cargo build --release --examples
```

That is the whole build. It finishes with **zero warnings**, downloads
**nothing**, and uses **no third-party libraries** — only Rust's standard
library. The code is compiled with `#![forbid(unsafe_code)]`, so the compiler
rejects any use of Rust's escape hatch for unchecked operations.

## 8. Step 3 — Run both and compare

Each program takes the simulation end time as its argument, and each writes a
file of raw bits into the current folder (`state_*_c.txt` from the C program,
`state_*_rust.txt` from the Rust one). The `./` prefix simply means "the
program in this folder"; macOS executables have no `.exe` suffix.

We compare with `cmp`, a built-in macOS tool that checks two files **byte by
byte**. With `-s` ("silent") it prints nothing and just reports success or
failure, so we chain it with `&&`/`||` to print a one-word verdict:
`IDENTICAL` if every byte matches, `DIFFERENT` otherwise.

Each simulation is run twice — once at a short end time, once at a long one:

```bash
cd ~/work/reboundx_rust/porttest

# ---- Test 1: pseudo-synchronisation — 10 orbits, then 100 orbits ----
./tides_spin_pseudo_c 62.83185307179586
../target/release/examples/tides_spin_pseudo 62.83185307179586
cmp -s state_pseudo_c.txt state_pseudo_rust.txt && echo IDENTICAL || echo DIFFERENT

./tides_spin_pseudo_c 628.3185307179586
../target/release/examples/tides_spin_pseudo 628.3185307179586
cmp -s state_pseudo_c.txt state_pseudo_rust.txt && echo IDENTICAL || echo DIFFERENT

# ---- Test 2: Kozai cycle — t = 1,000, then the example's full default t = 100,000 ----
./tides_spin_kozai_c 1000.0
../target/release/examples/tides_spin_kozai 1000.0
cmp -s state_kozai_c.txt state_kozai_rust.txt && echo IDENTICAL || echo DIFFERENT

./tides_spin_kozai_c 100000.0
../target/release/examples/tides_spin_kozai 100000.0
cmp -s state_kozai_c.txt state_kozai_rust.txt && echo IDENTICAL || echo DIFFERENT

# ---- Test 3: migration + obliquity tides — 10 orbits, then 100 orbits ----
./tides_spin_migration_c 62.83185307179586
../target/release/examples/tides_spin_migration 62.83185307179586
cmp -s state_migration_c.txt state_migration_rust.txt && echo IDENTICAL || echo DIFFERENT

./tides_spin_migration_c 628.3185307179586
../target/release/examples/tides_spin_migration 628.3185307179586
cmp -s state_migration_c.txt state_migration_rust.txt && echo IDENTICAL || echo DIFFERENT
```

If you ever see `DIFFERENT` and want to know *where*, run `diff` on the pair —
it prints the lines that disagree:

```bash
diff state_pseudo_c.txt state_pseudo_rust.txt
```

And to record a short "fingerprint" of a result file (a **SHA-256 hash** — a
64-character summary that changes if even one bit of the file changes):

```bash
shasum -a 256 state_pseudo_c.txt state_pseudo_rust.txt
```

Two identical files always print the same fingerprint.

> **Heads-up: the C prints a lot of warnings.** REBOUNDx warns, on *every
> single timestep*, that you are giving a velocity-dependent force to WHFast.
> That is correct and expected behaviour (our Rust prints it too), but for a
> long run it is hundreds of thousands of lines. Send them to nowhere by
> adding `2>/dev/null` — stream 2 is where a program's warnings go, and
> `/dev/null` is macOS's built-in discard-everything file:
>
> ```bash
> ./tides_spin_pseudo_c 628.3185307179586 2>/dev/null
> ```

## 9. The results

All six runs, measured on the machine in §5:

| Test | End time | Result |
|---|---|---|
| pseudo-synchronisation | 10 orbits (t = 62.83) | **BIT-IDENTICAL** |
| pseudo-synchronisation | 100 orbits (t = 628.3) | **BIT-IDENTICAL** |
| Kozai cycle | t = 1,000 | **BIT-IDENTICAL** |
| Kozai cycle | t = 100,000 (the example's full default) | **BIT-IDENTICAL** |
| migration + obliquity tides | 10 orbits (t = 62.83) | **BIT-IDENTICAL** |
| migration + obliquity tides | 100 orbits (t = 628.3) | **BIT-IDENTICAL** |

Every position, every velocity, every mass and **every spin vector** matched to
all 64 bits, in all six runs.

The Kozai result deserves a note. Because IAS15 chooses its own step size, and
because those choices depend on the extra forces, matching bit-for-bit at
t = 100,000 means the two programs took **the identical sequence of thousands
of adaptive steps** — not merely that they ended up in the same place.

For the record: the Windows 11 edition of this port passed the same six runs
bit-identically on its own platform. macOS is the second, independent
compiler-and-maths-library pairing on which the agreement holds.

## 10. Speed

The timings below are **Windows measurements**, kept from the Windows 11
edition of this write-up (same tests, MSVC-built C, same machine for both
columns). We did not re-time on macOS: speed is not the acceptance criterion
here, correctness is — but the lesson these numbers teach is worth keeping.

| Test | C (Windows) | Rust (Windows) |
|---|---|---|
| pseudo-synchronisation (t = 628.3) | 83.6 s | 0.9 s |
| Kozai (t = 100,000) | 1.3 s | 1.6 s |
| migration (t = 628.3) | 163.9 s | 1.7 s |

**Please do not read those first and third rows as "Rust is 90× faster than C".
They are not a fair comparison of the physics.** Nearly all of the C's time in
those two runs is spent *printing the per-timestep warning described above* —
hundreds of thousands of unbuffered writes. Rust buffers its output, so it pays
far less for the same messages.

The Kozai row (1.3 s vs 1.6 s) is the honest comparison, because that run is
short enough that printing does not dominate. The fair summary is: **the two
are in the same league**, with the Rust carrying extra safety guarantees for
free.

## 11. What this proves, and what it does not

**What it proves.** For everything these three simulations compute — the tidal
and spin forces, the spin differential equations, general-relativistic
precession, migration forces, the fixed-step and adaptive integrators, the
rotation into the invariable plane, and mid-run parameter changes — our Rust
REBOUNDx produces results identical to the original C, down to the last bit.

**What it does not prove.**

1. **Not every effect is covered.** REBOUNDx has about 20 effects. These three
   tests exercise `tides_spin`, `gr_potential` and `modify_orbits_forces`
   thoroughly, and the machinery all of them share (parameters, forces,
   operators, the extras state). The others are covered by the Rust test suite
   rather than by a bit-exact comparison against C.
2. **One machine, one compiler.** These runs are macOS on Apple Silicon, with
   Apple's clang and Apple's system maths library. A different platform's C
   library computes some maths functions with different last bits, so it would
   give a different — equally valid — reference. The acceptance criterion is
   that C and Rust agree *within* a platform, and here they do, exactly.

**What about the maths functions themselves?** We measured that too, on this
Mac: for each of **21** library maths functions (`sin`, `cos`, `tan`, `atan2`,
`sqrt`, `fmod`, `exp`, `log`, `cbrt`, `pow`, and more), C and Rust were
compared on 200,000 sample inputs each — and **every one, including `pow`, is
bit-identical**. Both the clang-built C and the Rust resolve every such call
to Apple's system maths library, so there is **no known divergent maths
function on macOS**. The full comparison is written up in `rebound_rust.md`.

> **Windows history: the `pow` caveat.** On the Windows 11 edition, Rust's
> `pow` (raise-to-a-power) and Microsoft's C library disagreed on about 0.03%
> of inputs, by at most 2 ULP (the smallest representable difference) — the
> single known divergence of that platform, and a real caveat there for any
> effect calling `pow` at run time. That caveat **does not exist on macOS**,
> as the measurement above shows.

## 12. Files involved

| File | What it is |
|---|---|
| `reboundx_rust/porttest/tides_spin_pseudo_c.c` | C harness, test 1 |
| `reboundx_rust/porttest/tides_spin_kozai_c.c` | C harness, test 2 |
| `reboundx_rust/porttest/tides_spin_migration_c.c` | C harness, test 3 |
| `reboundx_rust/examples/tides_spin_pseudo.rs` | Rust twin, test 1 |
| `reboundx_rust/examples/tides_spin_kozai.rs` | Rust twin, test 2 |
| `reboundx_rust/examples/tides_spin_migration.rs` | Rust twin, test 3 |
| `rebound_rust/porttest/macos_shim/rand_r_glibc.c` | the glibc `rand_r` shim linked into every C harness (§6d) |
| `reboundx_rust/porttest/msvc_shim/` | Windows history: the two VLA-free files the MSVC build needed (§6c); not used on macOS |
| `reboundx_rust/porttest/state_*_c.txt` | the C programs' raw-bit output |
| `reboundx_rust/porttest/state_*_rust.txt` | the Rust programs' raw-bit output |
| `reboundx_rust/tests/` | the Rust unit/integration test suite |

---

## Credit where it is due

**REBOUNDx** was written by **Dan Tamayo**, Hanno Rein and collaborators, and is
published under the GNU General Public License v3.

If you use it for published work, cite **Tamayo, Rein, Shi & Hernandez 2019**
(*Monthly Notices of the Royal Astronomical Society*, volume 491, page 2885;
[arXiv:1908.05634](https://arxiv.org/abs/1908.05634)), plus the paper for each
effect you switch on. For the `tides_spin` effect used here, that is **Lu,
Hernandez & Rein 2023** (MNRAS 526, 66), building on **Eggleton, Kiseleva & Hut
1998** and **Hut 1981**. REBOUND itself is **Rein & Liu 2012** (A&A 537, A128).

Everything here is a translation of their work. All the science is theirs.
