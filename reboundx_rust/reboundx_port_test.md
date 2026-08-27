# The REBOUNDx Tests: Did the Rust Version Get the Same Answer?

**Short answer: yes — every single bit, in all three test simulations.**

This document explains what REBOUNDx is, which three simulations we used to test
our Rust translation of it, exactly how we ran the test, and what we found. It
is written for someone who has never done this before. Every command is printed
in full, so you never need to open another document.

> **A note on `C:\work` in the commands below.** `C:\work` is a stand-in for
> whatever folder you put this project in — it is not a real folder, and you do
> not have to create one with that name. If your copy lives in
> `C:\Users\Sam\astronomy`, then read every `C:\work\...` below as
> `C:\Users\Sam\astronomy\...`. Only that leading part changes; the folder
> names after it are real and must match.

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
| Computer | Windows 11 Pro for Workstations, 64-bit (x86-64) |
| C compiler | Microsoft Visual C++ (`cl`) 19.51.36256 |
| Rust | `rustc` 1.91.1, `cargo` 1.91.1 |
| REBOUND (C) | version 5.1.1, commit `dad5f978` |
| REBOUNDx (C) | version 5.1.0 |

No Linux, no WSL2, no GCC, no Clang — everything native Windows.

## 6. Step 1 — Build the original C REBOUNDx

### 6a. Get both sources

```bash
cd C:\work
git clone https://github.com/hannorein/rebound.git rebound\rebound
git clone https://github.com/dtamayo/reboundx.git reboundx
```

### 6b. Build C REBOUND, and make a static library from it

```bash
cd C:\work\rebound\rebound\examples\shearing_sheet
cmd /c 'set PATH=C:\Program Files (x86)\GnuWin32\bin;%PATH% && "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" && make'

cd C:\work\rebound\rebound\src
cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && lib /nologo /OUT:librebound_static.lib *.obj'
```

> **Why a *static* library?** The normal build makes a DLL, which only exposes
> the functions marked as public API. Our test harnesses call a few internal
> ones, so we bundle all the object files into a static library instead.

### 6c. Build C REBOUNDx

```bash
cd C:\work\reboundx\src
cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /c /I"..\..\rebound\rebound\src" /I"." /D_GNU_SOURCE /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_WARNINGS /DLIBREBOUNDX /Ox /fp:precise *.c'
```

> ### ⚠ Two of the 33 files will not compile, and that is expected
>
> `gr_full.c` and `interpolation.c` use a C99 feature called **variable-length
> arrays** — arrays whose size is decided while the program runs:
>
> ```c
> double a_const[N][3];   /* gr_full.c   */
> double u[n];            /* interpolation.c */
> ```
>
> Microsoft's C compiler has never supported these (GCC and Clang do). This is
> the only genuine portability problem in either project.
>
> The fix changes **only where the memory comes from** — the heap instead of
> the stack — and touches no arithmetic:
>
> ```c
> double (*a_const)[3] = malloc((size_t)N*sizeof(*a_const));
> double* u            = malloc((size_t)n*sizeof(double));
> ```
>
> with matching `free()` calls so the lifetime matches. Same element type, same
> indices, same values, same order.
>
> To keep the upstream source untouched, the patched copies live separately, in
> `reboundx_rust/porttest/msvc_shim/`. Build them and add their object files:
>
> ```bash
> cd C:\work\reboundx_rust\porttest\msvc_shim
> cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /c /I"..\..\..\rebound\rebound\src" /I"..\..\..\reboundx\src" /D_GNU_SOURCE /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_WARNINGS /DLIBREBOUNDX /Ox /fp:precise gr_full.c interpolation.c'
> copy gr_full.obj ..\..\..\reboundx\src\
> copy interpolation.obj ..\..\..\reboundx\src\
> ```
>
> Note this affects the **C reference only**. The Rust port has no such problem:
> Rust's `Vec` is a growable array and needs no workaround.

Now bundle all 33 object files into the library:

```bash
cd C:\work\reboundx\src
cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && lib /nologo /OUT:libreboundx.lib *.obj'
```

### 6d. Build the three C test harnesses

The shipped examples print rounded numbers, use `<unistd.h>` (which does not
exist on Windows), and call `system("rm ...")`. Our harnesses are those examples
with exactly three changes, listed at the top of each file:

1. `<unistd.h>` dropped (nothing from it is used),
2. `system("rm ...")` removed,
3. the text output replaced by a final dump of every state variable as raw bits.

**The physics setup is byte-for-byte the stock example.**

```bash
cd C:\work\reboundx_rust\porttest
copy ..\..\reboundx\src\libreboundx.lib .
copy ..\..\rebound\rebound\src\librebound_static.lib .

cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /I"..\..\rebound\rebound\src" /I"..\..\reboundx\src" /D_GNU_SOURCE /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_WARNINGS /Ox /fp:precise tides_spin_pseudo_c.c libreboundx.lib librebound_static.lib /Fe:tides_spin_pseudo_c.exe'

cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /I"..\..\rebound\rebound\src" /I"..\..\reboundx\src" /D_GNU_SOURCE /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_WARNINGS /Ox /fp:precise tides_spin_kozai_c.c libreboundx.lib librebound_static.lib /Fe:tides_spin_kozai_c.exe'

cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /I"..\..\rebound\rebound\src" /I"..\..\reboundx\src" /D_GNU_SOURCE /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_WARNINGS /Ox /fp:precise tides_spin_migration_c.c libreboundx.lib librebound_static.lib /Fe:tides_spin_migration_c.exe'
```

> **Use PowerShell, not Git-Bash, for these.** Git-Bash's `cmd //c` can fail
> *silently* and leave an old program file in place, so you end up testing
> yesterday's build. This actually happened to us and cost real time.

## 7. Step 2 — Build the Rust version

```bash
cd C:\work\reboundx_rust
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
cd C:\work\reboundx_rust\porttest

.\tides_spin_pseudo_c.exe 62.83185307179586
..\target\release\examples\tides_spin_pseudo.exe 62.83185307179586

.\tides_spin_kozai_c.exe 1000.0
..\target\release\examples\tides_spin_kozai.exe 1000.0

.\tides_spin_migration_c.exe 62.83185307179586
..\target\release\examples\tides_spin_migration.exe 62.83185307179586
```

Each writes a file of raw bits. Compare them in PowerShell:

```powershell
Compare-Object (Get-Content state_pseudo_c.txt)    (Get-Content state_pseudo_rust.txt)
Compare-Object (Get-Content state_kozai_c.txt)     (Get-Content state_kozai_rust.txt)
Compare-Object (Get-Content state_migration_c.txt) (Get-Content state_migration_rust.txt)
```

**If `Compare-Object` prints nothing, the files are identical** — that is what
you want.

> **Heads-up: the C prints a lot of warnings.** REBOUNDx warns, on *every single
> timestep*, that you are giving a velocity-dependent force to WHFast. That is
> correct and expected behaviour (our Rust prints it too), but for a long run it
> is hundreds of thousands of lines. Send it to nowhere:
>
> ```powershell
> .\tides_spin_pseudo_c.exe 62.83185307179586 2>$null
> ```

## 9. The results

Each simulation was run twice: once at a short end time, and once at a long one.

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
t = 100,000 means the two programs took **the identical sequence of thousands of
adaptive steps** — not merely that they ended up in the same place.

## 10. Speed

Measured on the same machine, same runs:

| Test | C | Rust |
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
2. **`pow` can differ.** Rust's own `pow` (raise-to-a-power) and Microsoft's C
   library disagree on about 0.03% of inputs, by at most 2 ULP (the smallest
   representable difference). Every other maths function tested — `sin`, `cos`,
   `tan`, `atan2`, `sqrt`, `fmod`, `exp`, `log`, `cbrt` — is bit-identical. If
   an effect you enable calls `pow` at run time, expect eventual divergence in a
   chaotic system. The measurement behind that figure: 200,000 sampled inputs
   per function, compared bit for bit against the Microsoft C library on this
   machine. `pow` was the only function to disagree, and every disagreement
   found in the shaped (physically realistic) sample was exactly 1 ULP.
3. **One machine, one compiler.** These runs are Windows + MSVC. A different
   platform's C library would give a different reference.

## 12. Files involved

| File | What it is |
|---|---|
| `reboundx_rust/porttest/tides_spin_pseudo_c.c` | C harness, test 1 |
| `reboundx_rust/porttest/tides_spin_kozai_c.c` | C harness, test 2 |
| `reboundx_rust/porttest/tides_spin_migration_c.c` | C harness, test 3 |
| `reboundx_rust/examples/tides_spin_pseudo.rs` | Rust twin, test 1 |
| `reboundx_rust/examples/tides_spin_kozai.rs` | Rust twin, test 2 |
| `reboundx_rust/examples/tides_spin_migration.rs` | Rust twin, test 3 |
| `reboundx_rust/porttest/msvc_shim/` | the two MSVC-portable C files (§6c) |
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
