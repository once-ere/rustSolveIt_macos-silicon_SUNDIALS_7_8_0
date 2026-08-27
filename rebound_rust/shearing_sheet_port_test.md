# The Shearing-Sheet Test: Did the Rust Version Get the Same Answer as the Original?

**Short answer: yes — every single bit, for all 1,482 particles, after 400
time steps and 102,533 collisions.**

This document explains what we tested, how we tested it, what went wrong the
first time, how we tracked down the cause, and how you can run the whole test
yourself. It is written for someone who has never done this before. Every
command you need is printed in full. You will not have to look anything up in
another document.

> **A note on `C:\work` in the commands below.** `C:\work` is a stand-in for
> whatever folder you put this project in — it is not a real folder, and you do
> not have to create one with that name. If your copy lives in
> `C:\Users\Sam\astronomy`, then read every `C:\work\...` below as
> `C:\Users\Sam\astronomy\...`. Only that leading part changes; the folder
> names after it are real and must match.

---

## Table of contents

1. [What is REBOUND, and what is a "port"?](#1-what-is-rebound-and-what-is-a-port)
2. [What is the shearing sheet?](#2-what-is-the-shearing-sheet)
3. [Why this particular test?](#3-why-this-particular-test)
4. [What "bit-for-bit identical" means](#4-what-bit-for-bit-identical-means)
5. [The equipment: what is installed on this computer](#5-the-equipment-what-is-installed-on-this-computer)
6. [Step 1 — Build the original C program](#6-step-1--build-the-original-c-program)
7. [Step 2 — Build the Rust program](#7-step-2--build-the-rust-program)
8. [Step 3 — Run both and compare](#8-step-3--run-both-and-compare)
9. [The result](#9-the-result)
10. [The detective story: the first attempt FAILED](#10-the-detective-story-the-first-attempt-failed)
11. [What this proves, and what it does not](#11-what-this-proves-and-what-it-does-not)
12. [Files involved](#12-files-involved)

---

## 1. What is REBOUND, and what is a "port"?

**REBOUND** is a professional astronomy program. Astronomers use it to
calculate how objects move under gravity — planets around a star, moons around
a planet, or billions of ice chunks in Saturn's rings. It is written in the C
programming language by Hanno Rein and his collaborators, and it is used in
hundreds of published scientific papers.

A **port** means rewriting a program in a different programming language while
keeping exactly what it does. We rewrote REBOUND from **C** into **Rust**.

**Why bother?** C is fast but lets programmers make dangerous mistakes — for
example, reading past the end of a list, which can silently corrupt results or
crash the program. Rust refuses to compile code that could do those things. So
a Rust version gives you the same science with a large class of bugs made
impossible by the compiler.

**The catch:** if the rewrite changes even one number, it is not the same
program any more, and you could not trust its results. So we had to prove that
our Rust version computes *exactly* the same numbers. That is what this test
does.

---

## 2. What is the shearing sheet?

Imagine you want to simulate Saturn's rings. The rings contain something like
10^14 chunks of ice. No computer can track that many.

The **shearing sheet** is a clever shortcut. Instead of the whole ring, you
simulate a small rectangular patch of it — a "box" — and you assume the rest of
the ring looks the same. When a particle drifts out the left side of the box, a
copy comes back in on the right. Because material closer to Saturn orbits
faster than material further out, the box is *sheared*: the neighbouring boxes
slide past yours. Hence "shearing sheet".

This is the standard test problem that ships with REBOUND. In our run it
contains **1,482 ice particles** that:

- pull on each other by gravity,
- **bounce off each other** when they touch (102,533 bounces in our run),
- wrap around the box edges with the shear offset applied.

---

## 3. Why this particular test?

Because it uses almost every part of the program at once. Look at what has to
be exactly right for the final answer to match:

| Part of the program | What it does here |
|---|---|
| SEI integrator | moves particles forward in time in the rotating frame |
| Octree gravity | groups distant particles to compute gravity quickly |
| Tree collision search | finds which particles are touching |
| Hard-sphere collisions | makes them bounce, using a bounce law that depends on speed |
| Shear-periodic boundary | wraps particles around the box with the shear offset |
| Random number generator | places the particles at the start |

If **any single one** of those had a mistake — one wrong sign, one wrong loop
limit, one number computed in a different order — the final positions would
differ. With 102,533 collisions, an error introduced anywhere would spread to
every particle long before the end.

This is why we chose it as the acceptance test. It is the hardest single thing
in the program to get exactly right.

---

## 4. What "bit-for-bit identical" means

Computers store decimal numbers (like 3.14159...) in a format called
**IEEE-754 double precision** — 64 ones-and-zeros, called *bits*.

When we say two results are **bit-for-bit identical**, we mean all 64 bits are
the same for every number. Not "the same to 10 decimal places" — *the same*.

This matters more than it might sound. Simulations like this one are
**chaotic**: a difference in the very last bit (about 1 part in 10 million
billion) grows larger every step until the two runs look completely different.
So bit-for-bit agreement after 400 steps and 102,533 collisions is very strong
evidence that the two programs are doing identical arithmetic.

To check this, both programs dump every number as its raw bits, written in
hexadecimal (base 16). For example the number 1.0 is written
`3ff0000000000000`. Then we compare the two files.

> **A note on "one ULP".** ULP stands for *Unit in the Last Place* — the
> smallest possible difference between two neighbouring numbers the computer
> can represent. A "one ULP difference" is the smallest disagreement that can
> exist. It comes up in section 10.

---

## 5. The equipment: what is installed on this computer

This test was run on:

| Thing | Version |
|---|---|
| Computer | Windows 11 Pro for Workstations, 64-bit Intel/AMD (x86-64) |
| C compiler | Microsoft Visual C++ (`cl`) version 19.51.36256 |
| Where the C compiler lives | `C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\` |
| `make` (a build helper) | GnuWin32 Make 3.81 |
| Rust compiler | `rustc` 1.91.1 |
| Rust build tool | `cargo` 1.91.1 |
| Original REBOUND source | github.com/hannorein/rebound, version 5.1.1, commit `dad5f978` |

**No Linux, no WSL2, no GCC and no Clang were used anywhere.** Everything is
native Windows.

If you want to reproduce this you need:

1. **Visual Studio Build Tools** (free) — provides the C compiler `cl`.
   Download from <https://visualstudio.microsoft.com/downloads/>, pick
   "Build Tools for Visual Studio", and tick "Desktop development with C++".
2. **Rust** — download `rustup-init.exe` from <https://rustup.rs> and run it.
   Choose the default (MSVC) toolchain.
3. **GnuWin32 Make** — in a terminal, run:

   ```bash
   winget install GnuWin32.Make
   ```

---

## 6. Step 1 — Build the original C program

We need the original C program first, because it produces the "right answers"
that we compare against.

### 6a. Get the source code

```bash
cd C:\work
git clone https://github.com/hannorein/rebound.git rebound\rebound
```

### 6b. Compile it

```bash
cd C:\work\rebound\rebound\examples\shearing_sheet
cmd /c 'set PATH=C:\Program Files (x86)\GnuWin32\bin;%PATH% && "C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" && make'
```

That one line does three things: it puts `make` on the search path, it runs
`vcvars64.bat` (which tells the terminal where the C compiler is), and then it
runs `make`.

> ### ⚠ Two traps that will waste your afternoon
>
> **Trap 1 — the order matters.** You must add GnuWin32 to `PATH` *before*
> running `vcvars64.bat`, exactly as written above. If you do it the other way
> round it will fail with "cl is not recognized". The reason: Windows expands
> `%PATH%` when it *reads* the line, not when it runs it, so putting the
> GnuWin32 part second silently erases everything `vcvars64.bat` just added.
>
> **Trap 2 — use PowerShell, not Git-Bash.** Compile C only from PowerShell
> using `cmd /c '...'`. Git-Bash's `cmd //c` can fail *silently* and leave an
> old program file behind, so you end up testing yesterday's build and getting
> nonsense differences. This actually happened to us and cost real time.

When it finishes you will have `rebound.exe` plus the library files
`librebound.lib` and `librebound.dll` in that folder.

Each source file is compiled with these options, which matter for the test:

```bash
cl -c /DBUILDINGLIBREBOUND /D_GNU_SOURCE /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_WARNINGS /Ox /fp:precise -DSERVER -DGITHASH=dad5f97806ecbb408dcaff728851c64e67f9f6eb /Fo:rebound.obj rebound.c
```

- `/Ox` = optimise for speed.
- `/fp:precise` = **do not** rearrange floating-point arithmetic. This is
  essential: it tells the compiler to compute things in the order written.
- `-DSERVER` = include the built-in web viewer.
- OpenGL is switched off automatically on Windows by REBOUND's own build files.

### 6c. Build the C test harness

The stock example prints rounded numbers and uses a random starting seed. For a
bit-exact comparison we need a fixed seed and raw bits. `porttest/problem_test.c`
is the stock example with exactly three changes, listed at the top of that file.
Build it:

```bash
cd C:\work\rebound_rust\porttest
cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /I"..\..\rebound\rebound\src" /D_GNU_SOURCE /D_CRT_SECURE_NO_WARNINGS /D_CRT_NONSTDC_NO_WARNINGS /Ox /fp:precise problem_test.c librebound.lib /Fe:rebound_test.exe'
```

---

## 7. Step 2 — Build the Rust program

```bash
cd C:\work\rebound_rust
cargo build --release --example shearing_sheet_test
```

`cargo` is Rust's build tool. `--release` means "optimise for speed" (like
`/Ox`). `--example shearing_sheet_test` builds that one test program.

This finishes with **zero warnings**. The Rust code is compiled with
`#![forbid(unsafe_code)]`, which means the compiler *rejects* any use of Rust's
escape hatch for unchecked operations, and `#![deny(warnings)]`, which turns
every warning into a hard error. It uses **no third-party libraries at all**.

---

## 8. Step 3 — Run both and compare

Run each program for 400 time steps:

```bash
cd C:\work\rebound_rust\porttest
.\rebound_test.exe 400
..\target\release\examples\shearing_sheet_test.exe 400
```

Each writes two files: the starting state and the final state, as raw bits.

Now compare them. In PowerShell:

```powershell
Compare-Object (Get-Content state_c_init.txt)  (Get-Content state_rust_init.txt)
Compare-Object (Get-Content state_c_final.txt) (Get-Content state_rust_final.txt)
```

`Compare-Object` prints the lines that differ. **If it prints nothing, the
files are identical** — that is what you want.

As a second, independent check, take a *fingerprint* of each file. A SHA-256
hash is a 64-character code; if even one bit of the file changed, the code
changes completely:

```powershell
Get-FileHash state_c_final.txt    -Algorithm SHA256
Get-FileHash state_rust_final.txt -Algorithm SHA256
```

---

## 9. The result

```
shearing_sheet 400 steps: init=IDENTICAL final=IDENTICAL
C    SHA256: 75BDAAB7109F125192F56AEB0CCDCAC554AF88A165FE11075EC5A871178521F0
Rust SHA256: 75BDAAB7109F125192F56AEB0CCDCAC554AF88A165FE11075EC5A871178521F0
```

The two fingerprints are the same. **Every one of the 1,482 particles has the
identical position and velocity — all 64 bits of each of the 6 numbers — after
400 steps and 102,533 collisions.**

The run itself:

```
Toomre wavelength: 61.009899
N after init: 1482
final: t=1.91217633050229269e+04 steps=400 collisions=102533
```

**Speed.** On the same computer, the same 400 steps took:

| | Time |
|---|---|
| Original C | 2.2 seconds |
| Our Rust | 1.7 seconds |

The Rust version is slightly *faster* here, while doing the same arithmetic and
carrying the extra safety guarantees. (Do not read too much into one
measurement — across a range of integrators the two are within about ±30% of
each other, sometimes one ahead, sometimes the other.)

---

## 10. The detective story: the first attempt FAILED

It is worth telling you what happened the first time, because it shows how the
test is supposed to be used.

**Symptom.** The starting positions matched perfectly. But after 400 steps,
**330 of the 1,482 particles differed** — always in the last one or two bits.

Here is how we found the cause.

### Clue 1 — when does it first go wrong?

We ran both programs for 1 step, 2 steps, 4, 8, ... and compared each time.
Everything matched up to and including **step 77**. Step 78 was the first
mismatch. So one single event during step 78 started it.

### Clue 2 — what happened at step 78?

We printed every collision resolved during step 78 and compared. Exactly **one**
pair of particles — numbers **390 and 1456** — came out with different
velocities. Every input to that collision was bit-identical. Only the output
differed.

### Clue 3 — what is different about that collision?

When two ice particles bounce, REBOUND uses the **Bridges bounce law**, which
says how bouncy the collision is depending on the impact speed:

```
bounciness = 0.32 × (100 × speed) ^ (−0.234)
```

That `^` is the C function `pow` (raise to a power). It was the only unusual
function in the collision path.

### Clue 4 — testing every maths function

We wrote a test that calls each maths function 200,000 times in both C and Rust
and compares the raw bits:

```bash
cd C:\work\rebound_rust
cargo build --release --example libm_diff
cd porttest
cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /Ox /fp:precise libm_diff.c /Fe:libm_diff.exe'
.\libm_diff.exe
..\target\release\examples\libm_diff.exe
```

The result:

| Function | C vs Rust |
|---|---|
| `sin`, `cos`, `tan`, `atan2` | identical, all 200,000 |
| `sqrt`, `fmod`, `exp`, `log` | identical, all 200,000 |
| `cbrt` (cube root) | identical |
| **`pow`** | **60 of 200,000 differed** (0.03%), never by more than 2 ULP |

**`pow` is the one and only maths function where Rust and Microsoft's C library
disagree.** Rust ships its own `pow`; Microsoft's C library has a different
one. Both are correct to within the accuracy the standard requires — they just
round a handful of cases differently.

So the bug was not in our port at all. It was a difference between two vendors'
maths libraries.

### Clue 5 — the proof (the control experiment)

To be certain, we rewrote the bounce law using functions we had just proven
identical. These two lines compute the same thing:

```c
eps = 0.32 * pow(fabs(v)*100., -0.234);           /* uses pow  */
eps = 0.32 * exp(-0.234 * log(fabs(v)*100.));     /* uses exp and log */
```

We made that change **in both programs identically**, and re-ran.

**Result: all 400 steps bit-for-bit identical, matching SHA-256.**

That is the run reported in section 9, and it is why `porttest/problem_test.c`
and `examples/shearing_sheet_test.rs` both use the `exp`/`log` form. Swapping
the last remaining vendor difference out of the test made everything match —
which proves the port itself is exact, and the only difference in the whole
system was `pow`.

### The same story, seen a second time

Much later we found the identical fingerprint in a completely different place,
which is nice independent confirmation. The **BS integrator** (a different
method, not used in the shearing sheet) chooses its own step size using `pow`.
Running it for thousands of steps:

- steps 1 to 2,559: everything bit-identical;
- step 2,560: **all particle positions and velocities still bit-identical**,
  but the *proposed next step size* differed by exactly **one ULP**.

We then tested `pow` with exactly the arguments that step-size chooser uses
(200,000 samples):

```bash
cd C:\work\rebound_rust
cargo build --release --example bs_pow_diff
cd porttest
cmd /c '"C:\Program Files (x86)\Microsoft Visual Studio\18\BuildTools\VC\Auxiliary\Build\vcvars64.bat" >nul && cl /nologo /Ox /fp:precise bs_pow_diff.c /Fe:bs_pow_diff.exe'
.\bs_pow_diff.exe
..\target\release\examples\bs_pow_diff.exe
```

Result: **56 disagreements out of 200,000 (0.0280%), every single one exactly
one ULP.** Same function, same rate, same size. Case closed.

---

## 11. What this proves, and what it does not

**What it proves.** For everything REBOUND itself calculates in this
simulation — gravity, the tree, collisions, the boundary, the integrator, the
random numbers — the Rust port and the original C produce identical results
down to the last bit.

**What it does not prove.** It does not prove the two are identical when the
program calls `pow` at run time. They agree on about 99.97% of inputs and
differ by at most 2 ULP on the rest. Where a simulation is chaotic, that tiny
difference will eventually grow, and the two runs will drift apart — while both
remain equally valid answers to the same problem.

Places where `pow` gets called:

- user-written physics such as the stock Bridges bounce law,
- `reb_random_powerlaw` (drawing random sizes from a power law),
- the BS integrator's step-size chooser.

Everywhere else in REBOUND, agreement is exact.

**One honest caveat about this document's scope.** This is one test problem on
one computer. It is a very demanding test, but it is still a single
configuration. The broader evidence — 63 different integrator settings, all
bit-identical — is in the main provenance document `rebound_rust.md`, which
also repeats all of these commands so you never need two documents open.

---

## 12. Files involved

| File | What it is |
|---|---|
| `porttest/problem_test.c` | the C test harness (stock example + fixed seed + bit dump) |
| `examples/shearing_sheet_test.rs` | the Rust twin of that harness |
| `examples/shearing_sheet.rs` | a straight port of the *stock* example (uses `pow`, has the web viewer) |
| `porttest/libm_diff.c`, `examples/libm_diff.rs` | the maths-library comparison from section 10 |
| `porttest/bs_pow_diff.c`, `examples/bs_pow_diff.rs` | the BS step-size `pow` comparison |
| `porttest/state_c_init.txt`, `state_c_final.txt` | the C program's raw-bit output |
| `porttest/state_rust_init.txt`, `state_rust_final.txt` | the Rust program's raw-bit output |
| `notebooks/shearing_sheet_test.ipynb` | a Jupyter notebook that runs all of this and shows the result |

---

## Credit where it is due

REBOUND was written by **Hanno Rein** and collaborators and is published under
the GNU General Public License v3. The shearing-sheet method it implements is
described in **Rein & Tremaine 2011** (*Monthly Notices of the Royal
Astronomical Society*, volume 415, pages 3168–3176), and REBOUND itself in
**Rein & Liu 2012** (*Astronomy & Astrophysics*, volume 537, A128). If you use
this for published work, please cite them.

Everything here is a translation of their work. All the science is theirs.
