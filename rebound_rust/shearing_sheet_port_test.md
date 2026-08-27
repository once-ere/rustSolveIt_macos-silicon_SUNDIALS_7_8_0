# The Shearing-Sheet Test: Did the Rust Version Get the Same Answer as the Original?

**Short answer: yes — every single bit, for all 1,482 particles, after 400
time steps and every one of the enormous number of collisions along the way.**

This document explains what we tested, how we tested it, what went wrong the
first time — twice, once on each platform this port has lived on — how we
tracked down each cause, and how you can run the whole test yourself on a Mac.
It is written for someone who has never done this before. Every command you
need is printed in full. You will not have to look anything up in another
document.

> **A note on `~/work` in the commands below.** `~` is macOS shorthand for
> your home folder (for example `/Users/sam`), and `~/work` is a stand-in for
> whatever folder you put this project in — you do not have to create one with
> that exact name. If your copy lives in `~/Developer/astronomy`, then read
> every `~/work/...` below as `~/Developer/astronomy/...`. Only that leading
> part changes; the folder names after it are real and must match.

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
10. [The detective story, part one: Windows and the `pow` function](#10-the-detective-story-part-one-windows-and-the-pow-function)
11. [The detective story, part two: macOS and the random numbers](#11-the-detective-story-part-two-macos-and-the-random-numbers)
12. [What this proves, and what it does not](#12-what-this-proves-and-what-it-does-not)
13. [Files involved](#13-files-involved)

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

This port was first built and verified on Windows 11. The machine in front of
you now is a Mac with an Apple Silicon processor, so everything was re-built
and re-measured here from scratch. Every headline number in this document is
the macOS measurement; where a Windows number appears, it is clearly labelled
as history.

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
- **bounce off each other** when they touch — the program counts the bounces
  and prints the total at the end; on this Mac the run resolved **102,478**
  of them (the number is recorded in the `collisions_log_n` line both
  programs write into their final-state files; the Windows edition of this
  same test counted 102,533),
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
differ. With a hundred-thousand-odd collisions in the run, an error introduced
anywhere would spread to every particle long before the end.

This is why we chose it as the acceptance test. It is the hardest single thing
in the program to get exactly right. And as you will see in sections 10 and
11, two of the rows in that table have each starred in their own failure: the
hard-sphere bounce law on Windows, and the random number generator on macOS.

---

## 4. What "bit-for-bit identical" means

Computers store decimal numbers (like 3.14159...) in a format called
**IEEE-754 double precision** — 64 ones-and-zeros, called *bits*.

When we say two results are **bit-for-bit identical**, we mean all 64 bits are
the same for every number. Not "the same to 10 decimal places" — *the same*.

This matters more than it might sound. Simulations like this one are
**chaotic**: a difference in the very last bit (about 1 part in 10 million
billion) grows larger every step until the two runs look completely different.
So bit-for-bit agreement after 400 steps and 102,478 collisions is
very strong evidence that the two programs are doing identical arithmetic.

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
| Computer | Apple M5 Max (Apple Silicon, arm64), 128 GB RAM |
| Operating system | macOS Tahoe 26 (26.6.1, build 25G76) |
| C compiler | Apple clang 21.0.0 (clang-2100.1.1.101), from the Xcode Command Line Tools |
| Rust compiler | `rustc` 1.94.0 (aarch64-apple-darwin) |
| Rust build tool | `cargo` 1.94.0 |
| Original REBOUND source | github.com/hannorein/rebound, version 5.1.1, commit `dad5f978` |

**Everything is native macOS on Apple Silicon.** No virtual machines, no
Rosetta translation, no Linux — the C compiler and the Rust compiler both
produce arm64 code that runs directly on the M5 Max.

If you want to reproduce this you need two free things:

1. **Xcode Command Line Tools** — provides the C compiler `clang`, the
   library tool `ar`, and `git`. In Terminal, run:

   ```bash
   xcode-select --install
   ```

   and click Install in the window that appears.

2. **Rust** — install with the one-line command from <https://rustup.rs>:

   ```bash
   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
   ```

   Accept the default (aarch64-apple-darwin) toolchain.

---

## 6. Step 1 — Build the original C program

We need the original C program first, because it produces the "right answers"
that we compare against.

### 6a. Get the source code

```bash
cd ~/work
git clone https://github.com/hannorein/rebound.git rebound/rebound
cd rebound/rebound
git checkout dad5f97806ecbb408dcaff728851c64e67f9f6eb
```

The `git checkout` line pins the source to the exact version this port was
made from (5.1.1), so you are comparing against the same code we did.

### 6b. Compile the REBOUND library

```bash
cd ~/work/rebound/rebound/src
clang -c -DBUILDINGLIBREBOUND -D_GNU_SOURCE -DSERVER \
      -DGITHASH=dad5f97806ecbb408dcaff728851c64e67f9f6eb \
      -O2 -ffp-contract=off *.c
ar rcs librebound_static.a *.o
```

The first command compiles all 31 C source files (it finishes with no
errors); the second bundles the results into one library file,
`librebound_static.a` — think of it as a zip file of compiled code that other
programs can link against. What the options mean:

- `-c` = compile each file, don't try to make a runnable program yet.
- `-O2` = optimise for speed.
- `-DSERVER` = include the built-in web viewer.
- `-DGITHASH=...` = stamp the version number into the library.
- `-ffp-contract=off` = **do not fuse floating-point operations.** This is
  the essential one. Apple Silicon chips have a native instruction called FMA
  (*fused multiply-add*) that computes `a*b+c` in one step with one rounding
  instead of two — slightly *more* accurate, but with different last bits
  than doing the multiply and add separately. Compilers love to use it. Rust
  never fuses, so the C build must not either, or the two would disagree in
  the last bit everywhere. (On Windows the same job was done by the
  `/fp:precise` flag of MSVC, Microsoft's C compiler.)

### 6c. Build the C test harness — and the shim

The stock example prints rounded numbers and uses a random starting seed — a
*seed* is the starting number that determines the entire "random" sequence,
explained fully in section 11. For a bit-exact comparison we need a fixed
seed (we use **42**) and raw bits. `porttest/problem_test.c` is the stock
example with exactly three changes, listed at the top of that file.

On macOS there is one extra ingredient: a small helper file called
`macos_shim/rand_r_glibc.c` — its entire algorithm is about twenty lines of
arithmetic. Compile it first:

```bash
cd ~/work/rebound_rust/porttest
clang -c -O2 -ffp-contract=off macos_shim/rand_r_glibc.c -o macos_shim/rand_r_glibc.o
```

> ### ⚠ The trap that will waste your afternoon
>
> **Never build a C harness without linking `macos_shim/rand_r_glibc.o`.**
> If you leave it out, the program still builds and still runs — but it
> silently uses Apple's random number generator instead of the one the test
> was defined with, creates **1,441** particles instead of 1,482, and the
> comparison fails at the very first particle. This actually happened to us
> and is the whole subject of section 11. The command below includes the shim;
> keep it that way.

Now build the harness, linking the shim and the library from step 6b:

```bash
cd ~/work/rebound_rust/porttest
clang -I../../rebound/rebound/src -D_GNU_SOURCE -O2 -ffp-contract=off \
      problem_test.c macos_shim/rand_r_glibc.o \
      ../../rebound/rebound/src/librebound_static.a -lm -o problem_test
```

(`-I...` tells the compiler where REBOUND's header files are; `-lm` links the
system maths library; `-o problem_test` names the finished program.) You now
have a runnable program called `problem_test` in the `porttest` folder.

---

## 7. Step 2 — Build the Rust program

```bash
cd ~/work/rebound_rust
cargo build --release --example shearing_sheet_test
```

`cargo` is Rust's build tool. `--release` means "optimise for speed" (like
`-O2`). `--example shearing_sheet_test` builds that one test program.

This finishes with **zero warnings**. The Rust code is compiled with
`#![forbid(unsafe_code)]`, which means the compiler *rejects* any use of
Rust's escape hatch for unchecked operations, and `#![deny(warnings)]`, which
turns every warning into a hard error. It uses **no third-party libraries at
all**. And note: the Rust program needs no shim — you will see why in
section 11.

---

## 8. Step 3 — Run both and compare

Run each program for 400 time steps:

```bash
cd ~/work/rebound_rust/porttest
./problem_test 400
../target/release/examples/shearing_sheet_test 400
```

Each writes two files: the starting state and the final state, as raw bits.
Each also prints its particle count (`N after init: 1482`) and, at the end,
its step and collision totals.

Now compare the files. `cmp -s` compares two files byte by byte and says
nothing — it only sets a success/failure flag, which the `&& echo` turns into
a visible word:

```bash
cmp -s state_c_init.txt  state_rust_init.txt  && echo init IDENTICAL
cmp -s state_c_final.txt state_rust_final.txt && echo final IDENTICAL
```

**If a line prints `IDENTICAL`, that pair of files matches byte for byte.**
(If you ever want to *see* a difference, run `diff file1 file2` instead — it
prints the lines that disagree.)

As a second, independent check, take a *fingerprint* of each file. A SHA-256
hash is a 64-character code; if even one bit of the file changed, the code
changes completely:

```bash
shasum -a 256 state_c_final.txt state_rust_final.txt
```

---

## 9. The result

Measured on this machine:

```
init  IDENTICAL
final IDENTICAL

418c864dd1a610cbe8ea6d81ecafa1e4ce6d36837494177d9875ee820ef0766f  state_c_final.txt
418c864dd1a610cbe8ea6d81ecafa1e4ce6d36837494177d9875ee820ef0766f  state_rust_final.txt
```

The two fingerprints are the same. **Every one of the 1,482 particles has the
identical position and velocity — all 64 bits of each of the 6 numbers —
after 400 steps and all 102,478 collisions the run resolved.**

> **Why is this fingerprint different from the Windows one?** The Windows
> edition of this test reported SHA-256
> `75bdaab7109f125192f56aeb0ccdcac554af88a165fe11075ec5a871178521f0`. Both
> are correct. The bounce law calls the system maths library's `exp` and
> `log` functions, and Microsoft's and Apple's maths libraries round a few
> inputs differently — so the *trajectory* differs between platforms after
> the first bounce, while on each platform the C program and the Rust
> program (which share that platform's maths library) agree exactly. The
> acceptance criterion is C-versus-Rust agreement *on the same machine*, and
> on this Mac that agreement is bit-perfect.

**Speed.** We did not re-time the run for this document. For the record, on
the Windows 11 machine where the port was first developed, the same 400 steps
took 2.2 seconds in C and 1.7 seconds in Rust — the Rust version slightly
faster while doing the same arithmetic and carrying the extra safety
guarantees. (Do not read too much into one measurement; on that Windows
machine, across a range of integrators, the two were within about ±30% of
each other, sometimes one ahead, sometimes the other.)

---

## 10. The detective story, part one: Windows and the `pow` function

This test has failed twice in its life, once per platform, and both failures
are worth telling in full because they show how the test is supposed to be
used. Part one happened on Windows, where the port was first built. It is
history now — nothing in it needs fixing on your Mac — but it explains a line
of code you will find in both test programs to this day.

**Symptom (Windows, first ever run).** The starting positions matched
perfectly. But after 400 steps, **330 of the 1,482 particles differed** —
always in the last one or two bits.

### Clue 1 — when does it first go wrong?

We ran both programs for 1 step, 2 steps, 4, 8, ... and compared each time.
Everything matched up to and including **step 77**. Step 78 was the first
mismatch. So one single event during step 78 started it.

### Clue 2 — what happened at step 78?

We printed every collision resolved during step 78 and compared. Exactly
**one** pair of particles — numbers **390 and 1456** — came out with different
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

We wrote a test that calls each of nine maths functions 200,000 times in
both C and Rust and compares the raw bits. On Windows, the result was:

| Function | C vs Rust (Windows / MSVC) |
|---|---|
| `sin`, `cos`, `tan`, `atan2` | identical, all 200,000 |
| `sqrt`, `fmod`, `exp`, `log` | identical, all 200,000 |
| **`pow`** | **60 of 200,000 differed** (0.03%), never by more than 2 ULP |

**On Windows, `pow` was the one and only maths function where Rust and
Microsoft's C library disagreed.** Rust ships its own `pow` there;
Microsoft's C library has a different one. Both are correct to within the
accuracy the standard requires — they just round a handful of cases
differently.

So the bug was not in the port at all. It was a difference between two
vendors' maths libraries.

### Clue 5 — the proof (the control experiment)

To be certain, we rewrote the bounce law using functions we had just proven
identical. These two lines compute the same thing:

```c
eps = 0.32 * pow(fabs(v)*100., -0.234);           /* uses pow  */
eps = 0.32 * exp(-0.234 * log(fabs(v)*100.));     /* uses exp and log */
```

We made that change **in both programs identically**, and re-ran.

**Result: all 400 steps bit-for-bit identical, matching SHA-256.** That is
why `porttest/problem_test.c` and `examples/shearing_sheet_test.rs` both use
the `exp`/`log` form to this day.

The same fingerprint later showed up independently in the **BS integrator**
(a different method, not used in the shearing sheet), which chooses its own
step size using `pow`: on Windows, steps 1 to 2,559 were bit-identical, and
at step 2,560 all particle positions and velocities were *still*
bit-identical but the proposed next step size differed by exactly one ULP. A
200,000-sample test of `pow` with exactly the arguments that step-size
chooser uses found 56 disagreements (0.0280%), every one exactly one ULP.
Same function, same rate, same size. Case closed.

### And on macOS? The story evaporates.

On this Mac we re-ran the same 200,000-sample maths-function comparison. The
test dumps **nine functions** — `sin`, `cos`, `tan`, `atan2`, `pow`, `sqrt`,
`fmod`, `exp` and `log` — plus an extra `exp`/`log` cross-check appended to
the same dump. You can too — no shim needed, these programs call only the
maths library. Each program writes its dump to a file (`libm_c.txt` and
`libm_rust.txt`), and the final `cmp` line compares them:

```bash
cd ~/work/rebound_rust
cargo build --release --example libm_diff
cd porttest
clang -O2 -ffp-contract=off libm_diff.c -lm -o libm_diff
./libm_diff
../target/release/examples/libm_diff
cmp -s libm_c.txt libm_rust.txt && echo libm IDENTICAL
```

The measured macOS result: `libm IDENTICAL` — **C and Rust are bit-identical
for all nine functions, including `pow`.** On this platform both the
clang-compiled C and
the Rust program resolve every maths call, `pow` included, to Apple's maths
library, so there is nothing left to disagree about. The Windows port's "one
known difference" simply does not exist here.

The BS step-size check confirms it the same way:

```bash
cd ~/work/rebound_rust
cargo build --release --example bs_pow_diff
cd porttest
clang -O2 -ffp-contract=off bs_pow_diff.c -lm -o bs_pow_diff
./bs_pow_diff
../target/release/examples/bs_pow_diff
cmp -s bs_pow_c.txt bs_pow_rust.txt && echo bs_pow IDENTICAL
```

Measured macOS result: `bs_pow IDENTICAL` — **bit-identical, all 200,000
samples.**

We kept the `exp`/`log` form of the bounce law anyway — it is proven
identical on both platforms, and it keeps the two editions of the test
directly comparable.

---

## 11. The detective story, part two: macOS and the random numbers

So with the `pow` lesson already learned, bringing the test to the Mac should
have been routine: build the C library with clang, build the harness, run,
compare. Here is what the very first macOS run printed instead.

**Symptom.** The Rust program said `N after init: 1482`, as always. The C
program said `N after init: 1441`. Forty-one particles were simply missing —
and the comparison failed at **particle 1**, before a single time step had
been taken.

### Clue 1 — the failure is at the starting line

On Windows the starting states had always matched and the trouble began
mid-run. This was the opposite: the two programs disagreed about the *initial
conditions*. Only one part of the program builds the initial conditions — the
**random number generator** that places the particles.

"Random" here does not mean unpredictable. Simulations use a
**pseudo-random number generator**: a simple formula that, starting from a
chosen number called the **seed**, produces a long, fixed, repeatable stream
of scrambled-looking numbers. Same seed, same generator ⇒ same stream, every
time, on every machine. Both programs use seed **42** precisely so that they
will lay down identical particles. The shearing sheet keeps drawing random
positions until the patch reaches a target surface density, so the particle
*count* itself depends on the stream — 1,482 with the expected stream. A
count of 1,441 means the C program was drawing from a *different stream*.
Same seed, different generator.

### Clue 2 — read the original source

The generator REBOUND uses is a standard C function called `rand_r`. In the
original `rebound.c` there is a curious passage: a complete, hand-vendored
copy of `rand_r` from **glibc** (the GNU C library used on Linux), with a
comment citing its source
(`codebrowser.dev/glibc/glibc/stdlib/rand_r.c.html`) — wrapped in the
preprocessor guard `#ifdef _WIN32`, meaning "only compile this on Windows".

That explains everything:

- **On Linux**, `rand_r` comes from glibc itself.
- **On Windows**, Microsoft's C library has no `rand_r` at all, so REBOUND
  ships the glibc copy — the `#ifdef _WIN32` turns it on, and our Windows
  build had used it without us ever noticing.
- **On macOS**, the guard is false, so the vendored copy is *not* compiled —
  and Apple's C library *does* have a function named `rand_r`, so the build
  succeeds silently. But Apple's `rand_r` is a **different formula with the
  same name**. Seed 42 in, a completely different stream out.

Neither generator is "wrong" — both are legitimate `rand_r` implementations.
But our Rust port implements the *glibc* formula on every platform (that is
what all the reference results were made with, and it is why the Rust side
needs no shim and printed 1,482 everywhere). For the comparison to mean
anything, the C reference on macOS has to use the glibc formula too.

### The fix — a shim, not an edit

The original REBOUND source tree is a read-only reference: we do not edit it,
ever, because the whole point is to compare against *unmodified* upstream
code. So instead of touching `rebound.c`, we wrote
`porttest/macos_shim/rand_r_glibc.c`: the same vendored glibc algorithm,
compiled as its own little object file and placed on the link line of every C
harness (you did this in section 6c).

Why does that work? When the linker assembles the final program, it satisfies
each function name with the first definition it encounters. Our
`rand_r_glibc.o` appears on the command line *before* the system C library,
so REBOUND's call to `rand_r` is wired to the glibc formula and Apple's
version is never used. The upstream source stays untouched.

(If you are keeping score: Windows needed a small compatibility shim of its
own for the REBOUNDx half of this project — REBOUNDx is an add-on effects
package ported alongside REBOUND, covered in its own document — because MSVC
lacks a C99 feature called variable-length arrays that clang supports fine.
Each platform needs exactly one tiny shim, for opposite reasons — Windows
because a thing was missing, macOS because a thing was present but
different.)

### The result

With the shim linked: `N after init: 1482`, starting states identical, and —
as section 9 records — all 400 steps bit-for-bit identical with matching
SHA-256 fingerprints. The same fix also cured the only other test in the
suite that draws random numbers — the MEGNO/variational test
`movetocom_var` (MEGNO is a standard measure of how chaotic an orbit is) —
which is now bit-identical too.

Two platforms, two detective stories, one moral: **every mismatch this test
has ever shown traced back to the platform's system libraries, never to the
port.** Which is exactly what the test is for.

---

## 12. What this proves, and what it does not

**What it proves.** For everything REBOUND itself calculates in this
simulation — gravity, the tree, collisions, the boundary, the integrator, the
random numbers — the Rust port and the original C produce identical results
down to the last bit on this machine. And on macOS that agreement extends to
*every* maths-library function we tested, `pow` included: nine functions,
200,000 samples each, zero differing bits.

**What it does not prove.** It does not prove that a run on this Mac will
match a run on a Windows or Linux machine. Each platform's maths library
(`exp`, `log`, `pow`, ...) rounds a few inputs differently from the others,
and in a chaotic simulation those last-bit differences grow until the
trajectories visibly part ways — while both remain equally valid answers to
the same problem. That is why the Windows edition of this test has a
different SHA-256 fingerprint (section 9): the guarantee is C-equals-Rust
*within* a platform, not platform-equals-platform. On the Windows edition
there was the further caveat that C and Rust disagreed about `pow` at run
time; on macOS even that caveat is gone.

**One honest caveat about this document's scope.** This is one test problem
on one computer. It is a very demanding test, but it is still a single
configuration. The broader evidence — 63 different integrator configurations,
all 63 measured bit-identical on this Mac — is in the main provenance
document `rebound_rust.md`, which also repeats all of these commands so you
never need two documents open.

---

## 13. Files involved

| File | What it is |
|---|---|
| `porttest/problem_test.c` | the C test harness (stock example + fixed seed + bit dump) |
| `porttest/macos_shim/rand_r_glibc.c` | the glibc `rand_r` shim from section 11 — linked into every C harness on macOS |
| `examples/shearing_sheet_test.rs` | the Rust twin of that harness |
| `examples/shearing_sheet.rs` | a straight port of the *stock* example (uses `pow`, has the web viewer) |
| `porttest/libm_diff.c`, `examples/libm_diff.rs` | the maths-library comparison from section 10 |
| `porttest/bs_pow_diff.c`, `examples/bs_pow_diff.rs` | the BS step-size `pow` comparison |
| `porttest/state_c_init.txt`, `state_c_final.txt` | the C program's raw-bit output |
| `porttest/state_rust_init.txt`, `state_rust_final.txt` | the Rust program's raw-bit output |
| `notebooks/shearing_sheet_test.ipynb` | a Jupyter notebook that runs all of this and shows the result |

The four `state_*.txt` files are *outputs*, not sources: every run of the
section 8 commands recreates them. If one of them is missing from your copy
(or looks stale), that just means the programs have not been run yet — run
`./problem_test 400` and `../target/release/examples/shearing_sheet_test 400`
from the `porttest` folder and all four appear fresh.

---

## Credit where it is due

REBOUND was written by **Hanno Rein** and collaborators and is published under
the GNU General Public License v3. The shearing-sheet method it implements is
described in **Rein & Tremaine 2011** (*Monthly Notices of the Royal
Astronomical Society*, volume 415, pages 3168–3176), and REBOUND itself in
**Rein & Liu 2012** (*Astronomy & Astrophysics*, volume 537, A128). If you use
this for published work, please cite them.

Everything here is a translation of their work. All the science is theirs.
