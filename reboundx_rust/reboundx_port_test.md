# The REBOUNDx Tests: Did the Rust Version Get the Same Answer?

**Short answer: yes — every single bit, in all three test simulations, and
byte-for-byte in the save-file round trip (apart from one deliberately
blanked version-stamp field in the file header, explained in section 10).**

This document explains what REBOUNDx is, which three simulations we used to test
our Rust translation of it, exactly how we ran the test, and what we found. It
is written for someone who has never done this before. Every command is printed
in full, so you never need to open another document. Everything here was run on
a Mac with an Apple Silicon processor; the same test has also been run on
Windows 11, and where that history matters it is mentioned — clearly labelled
as the Windows result.

> **A note on `~/work` in the commands below.** The `~` symbol is shorthand
> that the Mac's Terminal understands for your *home folder* — the one named
> after your user account. `~/work` is a stand-in for whatever folder you put
> this project in — it is not a special folder, and you do not have to create
> one with that name. If your copy lives in `~/Documents/astronomy`, then read
> every `~/work/...` below as `~/Documents/astronomy/...`. Only that leading
> part changes; the folder names after it are real and must match.

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
10. [The save-file round trip](#10-the-save-file-round-trip)
11. [Speed](#11-speed)
12. [What this proves, and what it does not](#12-what-this-proves-and-what-it-does-not)
13. [Files involved](#13-files-involved)

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
ability to evolve spin vectors as extra differential equations (ODEs —
ordinary differential equations) alongside the orbits. An **integrator** is
the part of the program that advances the simulation through time, one small
step at a time; WHFast takes steps of a fixed size, while IAS15 (Test 2
below) chooses each step's size as it goes.

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
| Spin evolved as extra equations (ODEs) | yes | yes | yes |
| Rotation into the invariable plane | yes | yes | yes |
| Parameter changed mid-run | no | no | **yes** |

One row needs a word: before the run starts, both programs tilt the whole
system so that its total rotation — spins and orbits added together — points
straight up. The reference plane this defines is called the **invariable
plane**, and that setup calculation, too, must agree bit for bit.

## 4. What "bit-for-bit identical" means

Computers store decimal numbers as 64 ones-and-zeros (**bits**), in a format
called IEEE-754 double precision.

**Bit-for-bit identical** means *all 64 bits match*, for every number — not
"agrees to 10 decimal places".

This strictness is necessary because these systems are **chaotic**: a difference
in the very last bit grows with every step until the two runs look completely
different. So bit-for-bit agreement is the only test that really proves the
arithmetic is the same.

Both programs therefore print every number as raw bits in hexadecimal —
base-16 notation, where each character `0`–`9`, `a`–`f` stands for four
bits, so 64 bits fit in 16 characters. The number 1.0 prints as
`3ff0000000000000`. We then compare the files.

## 5. The equipment

| Thing | Version |
|---|---|
| Computer | Apple M5 Max (Apple Silicon, arm64), 128 GB RAM |
| Operating system | macOS Tahoe 26 (26.6.1) |
| C compiler | Apple clang 21.0.0, from the Xcode Command Line Tools |
| Rust | `rustc` 1.94.0, `cargo` 1.94.0 (aarch64-apple-darwin) |
| REBOUND (C) | version 5.1.1, commit `dad5f978` |
| REBOUNDx (C) | version 5.1.0 (the `5.1.0` release tag, commit `e884547a`) |

No Linux, no emulation, no GCC — everything native macOS on Apple Silicon.
All commands below are typed into the **Terminal** application (Applications →
Utilities → Terminal).

If you are starting from a fresh Mac, two installs give you every tool used
here. The first provides `clang` (the C compiler), `ar` (the library bundler)
and `git`; the second provides Rust:

```bash
xcode-select --install
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```

## 6. Step 1 — Build the original C REBOUNDx

### 6a. Get both sources

```bash
cd ~/work
git clone https://github.com/hannorein/rebound.git rebound/rebound
git clone https://github.com/dtamayo/reboundx.git reboundx
git -C rebound/rebound checkout dad5f97806ecbb408dcaff728851c64e67f9f6eb
git -C reboundx checkout 5.1.0
```

The two `git checkout` lines lock the sources to the exact versions this
document measured — REBOUND 5.1.1 at commit `dad5f978…`, and REBOUNDx at its
`5.1.0` release tag (commit `e884547a…`). Without them, the clones would give
you whatever the newest code happened to be on the day you ran them, and the
bit-for-bit results below might silently fail to reproduce. And yes, the
doubled folder name `rebound/rebound` is on purpose: REBOUNDx's build expects
to find REBOUND's source at `../rebound/rebound` relative to itself, so we
match that layout.

### 6b. Build C REBOUND, and make a static library from it

```bash
cd ~/work/rebound/rebound/src
clang -c -DBUILDINGLIBREBOUND -D_GNU_SOURCE -DSERVER \
      -DGITHASH=dad5f97806ecbb408dcaff728851c64e67f9f6eb \
      -O2 -ffp-contract=off *.c
ar rcs librebound_static.a *.o
```

The first command compiles all 31 C source files into object files (`.o`); the
second bundles them into one static library file. The `-DGITHASH=...` part just
stamps the build with the source-control version identifier, which the library
likes to know about.

> **Why a *static* library?** The normal build makes a shared library (a
> `.dylib` on macOS), which only exposes the functions marked as public API.
> Our test harnesses (the small test programs built in step 6d below) call a
> few internal ones, so we bundle all the object files into a static library
> instead.

> ### ⚠ The one flag that really matters: `-ffp-contract=off`
>
> Apple Silicon processors have a **fused multiply-add** instruction: it
> computes a×b+c in a single step, rounding the result once instead of twice.
> That single skipped rounding changes the very last bit of the answer. C
> compilers are normally allowed to "fuse" multiplications and additions this
> way whenever they like; Rust never does it on its own. If the C build fused
> and the Rust build did not, the two could disagree in the last bit while
> both being perfectly reasonable — and in a chaotic system that last bit
> grows. `-ffp-contract=off` forbids the fusing, so both sides round the same
> way. It is the macOS twin of the `/fp:precise` flag the Windows 11 edition
> of this test used with Microsoft's compiler. Every C compile command in this
> document carries it.

### 6c. Build C REBOUNDx

```bash
cd ~/work/reboundx/src
clang -c -I../../rebound/rebound/src -I. -D_GNU_SOURCE -DLIBREBOUNDX \
      -O2 -ffp-contract=off *.c
ar rcs libreboundx.a *.o
```

All 33 source files compile, unmodified, with no warnings.

> ### A piece of Windows history: the two files that would not compile there
>
> On the Windows 11 edition of this test, two of the 33 files — `gr_full.c`
> and `interpolation.c` — would not compile, because they use a C feature
> called **variable-length arrays** (arrays whose size is decided while the
> program runs) that Microsoft's compiler has never supported. Patched copies,
> which take the memory from the heap instead of the stack and change no
> arithmetic, live in `reboundx_rust/porttest/msvc_shim/` and were used for
> the Windows build only. Apple's clang supports variable-length arrays, so on
> macOS **no patching is needed** — the upstream files compile exactly as
> their authors wrote them, and the `msvc_shim` folder is not used at all.

### 6d. Build the three C test harnesses

The shipped examples print rounded numbers, include a header file that Windows
lacks, and delete an old output file by calling out to a shell. Our harnesses
are those examples with exactly three changes, listed at the top of each file:

1. the `<unistd.h>` include dropped (nothing from it is used; it was only ever
   a Windows portability problem),
2. the shell call that deleted an old output file replaced by a direct
   file-delete call (no shell dependency),
3. the text output replaced by a final dump of every state variable as raw
   bits.

**The physics setup is byte-for-byte the stock example.**

```bash
cd ~/work/reboundx_rust/porttest

clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off tides_spin_pseudo_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.c \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a \
      -lm -o tides_spin_pseudo_c

clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off tides_spin_kozai_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.c \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a \
      -lm -o tides_spin_kozai_c

clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off tides_spin_migration_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.c \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a \
      -lm -o tides_spin_migration_c
```

Each command compiles one harness and links it with the REBOUNDx library, the
REBOUND library, and one small extra file explained in the box below. Expect
one harmless warning while compiling `tides_spin_pseudo_c.c` (a pointer-type
cast inside the harness's own bookkeeping code); the libraries themselves
compile warning-free.

> ### Why `rand_r_glibc.c` is on every link line — the macOS `rand_r` story
>
> `rand_r` is a small random-number generator that REBOUND uses when a
> simulation starts from *random* initial conditions. REBOUND wants the same
> seed to produce the same "random" numbers on every platform, so on Windows —
> whose C library has no `rand_r` at all — the upstream source carries its own
> copy of the standard Linux (glibc) algorithm. On macOS that copy is switched
> off, and the C build silently uses Apple's built-in `rand_r` instead — a
> **different** generator that turns the same seed into a different stream of
> numbers. During this port, that mismatch broke a REBOUND-side test that
> builds 1,482 particles from a random seed: the C reference suddenly built
> 1,441 and the comparison died at particle 1. The fix is
> `rebound_rust/porttest/macos_shim/rand_r_glibc.c` — the same glibc algorithm
> as its own little file, linked into every C harness so that the linker
> resolves REBOUND's `rand_r` calls there instead of in Apple's library. The
> upstream source is untouched. The three tides simulations in *this* document
> use no random numbers at all — their starting conditions are written out
> digit by digit — but every harness links the shim anyway, so every harness
> draws from the same generator by construction. The Rust port implements the
> glibc algorithm on every platform, so it needs no shim.

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

Each program takes the simulation end time as its argument.

```bash
cd ~/work/reboundx_rust/porttest

./tides_spin_pseudo_c 62.83185307179586
../target/release/examples/tides_spin_pseudo 62.83185307179586

./tides_spin_kozai_c 1000.0
../target/release/examples/tides_spin_kozai 1000.0

./tides_spin_migration_c 62.83185307179586
../target/release/examples/tides_spin_migration 62.83185307179586
```

Each writes a file of raw bits into the current folder — the C programs write
`state_pseudo_c.txt`, `state_kozai_c.txt` and `state_migration_c.txt`; the
Rust programs write the matching `state_*_rust.txt` files. Compare them with
`cmp`, the Mac's byte-by-byte file comparer — it stays silent when the files
are identical, so we add an `echo` to make success visible:

```bash
cmp state_pseudo_c.txt    state_pseudo_rust.txt    && echo pseudo identical
cmp state_kozai_c.txt     state_kozai_rust.txt     && echo kozai identical
cmp state_migration_c.txt state_migration_rust.txt && echo migration identical
```

**If each line prints its "identical" message, the files match byte for byte**
— that is what you want. If two files differ, `cmp` names the first byte where
they disagree instead.

For the long runs in the results table below, repeat the same pairs of commands
with the long end times: `628.3185307179586` for the first and third tests, and
`100000.0` for the Kozai test.

> **Heads-up: the C prints a lot of warnings.** REBOUNDx warns, on *every single
> timestep*, that you are giving a velocity-dependent force to WHFast. That is
> correct and expected behaviour (our Rust prints it too), but for a long run it
> is hundreds of thousands of lines. Send it to nowhere:
>
> ```bash
> ./tides_spin_pseudo_c 62.83185307179586 2>/dev/null
> ```

## 9. The results

Each simulation was run twice on this Mac: once at a short end time, and once
at a long one. The end times are in the simulation's own time units, in which
a body orbiting at distance 1 takes exactly 2π to go around once — so
t = 62.83 is 10 × 2π, ten of those reference "years", and t = 628.3 is a
hundred. (The close-in planets in these tests orbit much faster than that
reference body, so they complete far more orbits than ten in that span.)

| Test | End time | Result |
|---|---|---|
| pseudo-synchronisation | short run (t = 62.83 = 10 × 2π) | **BIT-IDENTICAL** |
| pseudo-synchronisation | long run (t = 628.3 = 100 × 2π) | **BIT-IDENTICAL** |
| Kozai cycle | t = 1,000 | **BIT-IDENTICAL** |
| Kozai cycle | t = 100,000 (the example's full default) | **BIT-IDENTICAL** |
| migration + obliquity tides | short run (t = 62.83 = 10 × 2π) | **BIT-IDENTICAL** |
| migration + obliquity tides | long run (t = 628.3 = 100 × 2π) | **BIT-IDENTICAL** |

Every position, every velocity, every mass and **every spin vector** matched to
all 64 bits, in all six runs.

The Kozai result deserves a note. Because IAS15 chooses its own step size, and
because those choices depend on the extra forces, matching bit-for-bit at
t = 100,000 means the two programs took **the identical sequence of thousands of
adaptive steps** — not merely that they ended up in the same place.

## 10. The save-file round trip

REBOUNDx can save its entire configuration — which effects are switched on,
every parameter attached to every force, operator and particle — to a **binary
file** (a file of raw bytes rather than readable text), and load it back later.
The Rust port writes and reads the same file format, and this test checks that
claim from both directions.

The test state is deliberately awkward: three particles, two forces
(`gr_potential` and `central_force`), two **operators** — effects that,
instead of applying a continuous force, make a discrete change to the system
just before or just after each timestep (here `modify_mass` and `drift`, one
scheduled before and one after) — and parameters of several
types — decimal numbers, whole numbers, a 3-component vector, and a parameter
that *points at* one of the forces. Two small C harnesses drive the C side; the
Rust side is one example program that writes, reads back, and checks
everything.

Build the two C harnesses (same pattern as section 6d):

```bash
cd ~/work/reboundx_rust/porttest

clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off rebx_binary_roundtrip_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.c \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a \
      -lm -o rebx_binary_roundtrip_c

clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off rebx_binary_read_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.c \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a \
      -lm -o rebx_binary_read_c
```

Now run the round trip. The C writes its file; the Rust writes its own file
*and* reads the C's file back, checking every value by its bit pattern:

```bash
./rebx_binary_roundtrip_c rebx_c_reference.bin
../target/release/examples/rebx_binary_roundtrip rebx_c_reference.bin
```

The Rust program prints one `PASS` line per check and finishes with:

```
25 / 25 checks passed
ROUND TRIP PASSED
```

Those 25 checks cover every list, every ordering, and every parameter value —
and the last of them re-saves the state the Rust just loaded from the C's file
and confirms the re-saved bytes are identical to its own file, byte for byte.

Next, compare the two files directly:

```bash
ls -l  rebx_c_reference.bin rebx_binary_roundtrip.bin
cmp -l rebx_c_reference.bin rebx_binary_roundtrip.bin | wc -l
cmp -l rebx_c_reference.bin rebx_binary_roundtrip.bin
```

Measured on this Mac: **both files are exactly 10,784 bytes**, and `cmp -l`
(which lists every differing byte) reports **26 differences, at byte positions
38 through 63** (`cmp` counts the first byte as 1; counting from zero these
are bytes 37 through 62 of the file) — and nowhere else. Those 26
bytes are the file header's *githash field*, a stamp recording which
source-control version of the writer produced the file: the C writes the text
`notavailable` padded out, and the Rust deliberately writes zeros there. This
is a documented, deliberate deviation — the field is informational, and both
readers ignore it. Every byte that carries physics or parameters is identical.
(An earlier revision of this test program set up a smaller state and wrote a
6,392-byte file; if you see that figure in older notes, it refers to that
earlier revision, not to the current program, which writes 10,784 bytes on
both platforms.)

Finally, the reverse direction: the original **C library reads the Rust's
file**, and as a control also reads its own, dumping everything as raw bits:

```bash
./rebx_binary_read_c rebx_c_reference.bin      > readback_from_c.txt
./rebx_binary_read_c rebx_binary_roundtrip.bin > readback_from_rust.txt
cmp readback_from_c.txt readback_from_rust.txt && echo readbacks identical
```

The two 28-line dumps are **identical**: the C library sees exactly the same
forces, operators, orderings and parameter bit patterns in the Rust-written
file as in its own.

## 11. Speed

The timings below are from the **Windows 11 edition** of this test — measured
on that machine, with Microsoft's compiler. We did not re-time the runs for
this macOS document, and the lesson the table teaches is about printing, not
about platforms or arithmetic, so we quote it as history:

| Test (Windows 11 measurements) | C | Rust |
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
short enough that printing does not dominate. The fair summary is: **the two are
in the same league**, with the Rust carrying extra safety guarantees for free.

## 12. What this proves, and what it does not

**What it proves.** For everything these three simulations compute — the tidal
and spin forces, the spin differential equations, general-relativistic
precession, migration forces, the fixed-step and adaptive integrators, the
rotation into the invariable plane, and mid-run parameter changes — our Rust
REBOUNDx produces results identical to the original C, down to the last bit.
And the save-file format round-trips byte-for-byte in both directions.

**What it does not prove.**

1. **Not every effect is covered.** REBOUNDx has about 20 effects. These three
   tests exercise `tides_spin`, `gr_potential` and `modify_orbits_forces`
   thoroughly, and the machinery all of them share (parameters, forces,
   operators, the extras state). The others are covered by the Rust test suite
   rather than by a bit-exact comparison against C.
2. **On this Mac, no maths-library difference was found — not even `pow`.**
   We measured it directly: 200,000 sampled inputs for each of the nine
   mathematical functions the comparison harness exercises (`sin`, `cos`,
   `tan`, `atan2`, `pow`, `sqrt`, `fmod`, `exp`, `log`), plus a separate
   `exp`/`log` differential check, C against Rust, compared bit for bit —
   **every function was bit-identical, including `pow`** (raise-to-a-power). The
   reason: on macOS, both the clang-built C and the Rust resolve every one of
   those calls to Apple's system maths library, so they cannot disagree. This
   is a genuine platform difference from the Windows 11 edition of this test,
   where `pow` — Microsoft's C-library version against Rust's — disagreed on
   about 0.03% of inputs by at most 2 ULP (the smallest representable
   difference), and was the one known divergence there. That Windows caveat
   simply does not exist on this platform.
3. **One machine, one compiler.** These runs are macOS on Apple Silicon, with
   Apple's clang and Apple's system libraries. A different platform's C
   library would give a different reference — the Windows 11 edition of this
   same test is the demonstration, reaching the same bit-for-bit verdict
   against *its* platform's C. The agreement that matters is always measured
   within a platform, and on this one it is exact.

## 13. Files involved

| File | What it is |
|---|---|
| `reboundx_rust/porttest/tides_spin_pseudo_c.c` | C harness, test 1 |
| `reboundx_rust/porttest/tides_spin_kozai_c.c` | C harness, test 2 |
| `reboundx_rust/porttest/tides_spin_migration_c.c` | C harness, test 3 |
| `reboundx_rust/examples/tides_spin_pseudo.rs` | Rust twin, test 1 |
| `reboundx_rust/examples/tides_spin_kozai.rs` | Rust twin, test 2 |
| `reboundx_rust/examples/tides_spin_migration.rs` | Rust twin, test 3 |
| `reboundx_rust/porttest/rebx_binary_roundtrip_c.c` | C writer for the save-file test (§10) |
| `reboundx_rust/porttest/rebx_binary_read_c.c` | C reader/dumper for the save-file test (§10) |
| `reboundx_rust/examples/rebx_binary_roundtrip.rs` | Rust save-file round-trip checker (§10) |
| `rebound_rust/porttest/macos_shim/rand_r_glibc.c` | the glibc `rand_r` shim linked into every C harness (§6d) |
| `reboundx_rust/porttest/msvc_shim/` | Windows-only patched files (§6c history; not used on macOS) |
| `reboundx_rust/porttest/state_*_c.txt` | the C programs' raw-bit output |
| `reboundx_rust/porttest/state_*_rust.txt` | the Rust programs' raw-bit output |
| `reboundx_rust/porttest/rebx_*.bin`, `readback_from_*.txt` | the save-file test's files and readback dumps |
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
