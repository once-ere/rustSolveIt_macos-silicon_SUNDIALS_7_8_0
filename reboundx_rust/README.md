# reboundx_rust (`reboundx_rs`)

A pure-Rust translation of **[REBOUNDx](https://github.com/dtamayo/reboundx)
5.1.0** — "REBOUND eXtras", the library by **Dan Tamayo**, Hanno Rein and
collaborators that adds extra physics (general relativity, tides, migration,
radiation pressure, and more) to an N-body simulation run by
**[REBOUND](https://github.com/hannorein/rebound)**.

It is the companion of
**[`rebound_rs`](https://github.com/once-ere/rebound_rust)**, the pure-Rust
translation of REBOUND itself, and it depends on it exactly as the C
`libreboundx` links against `librebound`.

- Zero `unsafe`, zero external dependencies (Rust's standard library only),
  zero build warnings, zero `clippy` warnings.
- C function and struct names preserved (`rebx_attach`, `rebx_add_force`,
  `rebx_set_param_double`, `rebx_tides_spin`, …), so the upstream REBOUNDx
  documentation still tells you what everything does.
- Verified against the C REBOUNDx compiled with `clang` — Apple's C
  compiler — on the same Apple-Silicon Mac: every tested result matches
  **bit for bit** — every single bit of every number, not "agrees to ten
  decimal places". (The same comparison was first done on Windows 11
  against Microsoft's compiler; the verification section below tells both
  stories.)

**All of the science, and all of the credit, belongs to the original
authors.** This is a translation, not new research.

---

## Before you start: this needs TWO folders

**If you cloned the `rustSolveIt_macos-silicon_SUNDIALS_7_8_0` repository,
you already have the right arrangement.** `rebound_rust/` and
`reboundx_rust/` are both top-level folders of that repository, sitting side
by side, and `cargo build` inside `reboundx_rust` works immediately. The
rest of this section is for anyone who copies `reboundx_rust` out on its
own, or clones the two crates from their stand-alone repositories.

**`cargo build` will fail if you take only this folder.** That is not a
bug — it is the same arrangement the original C library uses, and it is worth
one minute of your time to understand.

REBOUNDx is not a program on its own. It is an *add-on* to REBOUND, and it
needs REBOUND's code sitting right next to it. The C version enforces exactly
the same thing: if you try to build C REBOUNDx without REBOUND beside it, its
Makefile stops with

> `REBOUNDx not in the same directory as REBOUND.`

This Rust port keeps that convention. Its `Cargo.toml` says:

```toml
[dependencies]
rebound_rs = { path = "../rebound_rust" }
```

`../rebound_rust` means "the folder named `rebound_rust`, one level up from
here". So the two folders must be **siblings** — side by side in the same
parent folder, like this:

```
some-folder/
├── rebound_rust/     <- github.com/once-ere/rebound_rust
└── reboundx_rust/    <- this crate
```

### Do this, and it will work

```bash
mkdir -p ~/work/rebound-rust-ports
cd ~/work/rebound-rust-ports
git clone https://github.com/once-ere/rebound_rust.git
git clone https://github.com/once-ere/reboundx_rust.git
cd reboundx_rust
cargo build --release
```

If you skip the first `git clone`, you get this error, which names a folder
that does not exist:

```
error: failed to load manifest for dependency `rebound_rs`

Caused by:
  failed to read `/Users/you/work/rebound-rust-ports/rebound_rust/Cargo.toml`

Caused by:
  No such file or directory (os error 2)
```

The fix is always the same: put `rebound_rust` next to `reboundx_rust`.

---

## What you need installed

| You need | Why | Get it |
|---|---|---|
| **Rust** (1.91 or newer; this port was built and verified with 1.94.0) | compiles the code | <https://rustup.rs> — the installer does everything |
| **Git** | downloads the repositories | comes with Apple's Xcode Command Line Tools: run `xcode-select --install` in Terminal |

Nothing else. No C compiler, no Python, no libraries to hunt down — this crate
depends on nothing outside Rust's own standard library and its sibling
`rebound_rs`. (One exception: if you want to *re-run the C-versus-Rust
verification* described below, you need `clang`, Apple's C compiler — which
the same `xcode-select --install` already gave you.)

Check Rust is installed:

```bash
cargo --version
```

---

## Five-minute start: watch Mercury's orbit precess

This is the classic test of general relativity. Mercury's orbit slowly rotates
in a way Newton's gravity alone cannot explain — 43 arcseconds per century.
Switch on REBOUNDx's `gr_potential` force and you can watch it happen.

Make a new project **beside the two folders**:

```bash
cargo new mercury
cd mercury
```

Put this in `mercury/Cargo.toml`:

```toml
[workspace]

[package]
name = "mercury"
version = "0.1.0"
edition = "2021"

[dependencies]
rebound_rs  = { path = "../rebound_rust" }
reboundx_rs = { path = "../reboundx_rust" }
```

(The bare `[workspace]` line on top is worth keeping. It tells Cargo "this
project stands alone", which stops it from getting confused by unrelated
project files that may sit higher up in your folders.)

Put this in `mercury/src/main.rs`:

```rust
//! Minimal REBOUNDx demo: Mercury, plus the general-relativistic
//! correction to the Sun's potential, which makes its perihelion precess.

use rebound_rs::*;
use reboundx_rs::*;

fn main() {
    // 1. A simulation. G = 1, so the units are AU, solar masses, yr/2pi.
    let mut sim = reb_simulation_create();

    // 2. Two particles: a 1 solar-mass star, and Mercury on its real orbit.
    reb_simulation_add_fmt(&mut sim, "m", &[reb_fmt_arg::d(1.0)]);
    reb_simulation_add_fmt(
        &mut sim,
        "m a e",
        &[
            reb_fmt_arg::d(1.66e-7),   // mass, in solar masses
            reb_fmt_arg::d(0.387098),  // semi-major axis, in AU
            reb_fmt_arg::d(0.205630),  // eccentricity (0 = circle)
        ],
    );
    reb_simulation_move_to_com(&mut sim);

    // 3. Attach REBOUNDx and switch on one extra force.
    rebx_attach(&mut sim);
    let gr = rebx_load_force(&mut sim, "gr_potential").expect("gr_potential is a known force");
    rebx_add_force(&mut sim, gr);

    // 4. gr_potential needs exactly one parameter: the speed of light,
    //    which is 10065.32 in these units.
    if let Some(rebx) = rebx_extras_mut(&mut sim) {
        rebx_set_param_double(rebx, rebx_ap::force(gr), "c", 10065.32);
    }

    // 5. Integrate for 200 orbits of Mercury.
    let orbits = 200.0;
    let before = reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0]);
    reb_simulation_integrate(&mut sim, orbits * before.P);

    // 6. Print what general relativity did to the orbit.
    let after = reb_orbit_from_particle(sim.G, sim.particles[1], sim.particles[0]);
    println!("after t = {:.4} ({} orbits)", sim.t, orbits);
    println!("  a      = {:.9} AU   (was {:.9})", after.a, before.a);
    println!("  e      = {:.9}      (was {:.9})", after.e, before.e);
    println!(
        "  pomega = {:.6e} rad  (was {:.6e})  <- GR perihelion precession",
        after.pomega, before.pomega
    );
}
```

Run it:

```bash
cargo run --release
```

You will see:

```
after t = 302.6504 (200 orbits)
  a      = 0.387098000 AU   (was 0.387098000)
  e      = 0.205630000      (was 0.205630000)
  pomega = 1.003740e-4 rad  (was 0.000000e0)  <- GR perihelion precession
```

Read that last line. `pomega` is the direction the orbit points. It started at
zero and has rotated. The size and shape of the orbit (`a` and `e`) did not
change at all — general relativity does not shrink Mercury's orbit, it turns
it. Divide the rotation by 200 orbits and convert to arcseconds per century
and you get the textbook 43 arcsec/century. **The force really is doing
physics.**

---

## What you can switch on

REBOUNDx effects come in two kinds. A **force** adds an acceleration to the
particles at every step. An **operator** does something else to the simulation
between steps — changing a mass, nudging an orbit, recording a distance.

You choose one by passing its name as text.

### The 17 forces

Load with `rebx_load_force(&mut sim, "name")`, then `rebx_add_force`.

| Area | Names you can use |
|---|---|
| General relativity | `gr`, `gr_full`, `gr_potential`, `lense_thirring` |
| Tides and spin | `tides_spin`, `tides_constant_time_lag`, `tides_dynamical` |
| Migration and disks | `modify_orbits_forces`, `type_I_migration`, `exponential_migration`, `gas_damping_timescale`, `gas_dynamical_friction` |
| Radiation | `radiation_forces`, `yarkovsky_effect` |
| Shape of the central body | `gravitational_harmonics` |
| Other | `central_force`, `stochastic_forces` |

### The 10 operators

Load with `rebx_load_operator(&mut sim, "name")`, then `rebx_add_operator_step`.

| Area | Names you can use |
|---|---|
| Mass change | `modify_mass` |
| Orbit change | `modify_orbits_direct` |
| Measurement | `track_min_distance` |
| Applying a force as an operator | `integrate_force` |
| Building your own integrator | `drift`, `kick`, `kepler`, `jump`, `interaction`, `ias15` |

These names are matched **exactly and case-sensitively**. If you misspell one,
the library prints

```
REBOUNDx error: Force 'gr_potentail' not found in REBOUNDx library.
```

and returns nothing, rather than silently doing the wrong thing. That is why
the starter program above calls `.expect(...)` — it turns that "nothing" into
a clear stop.

There is also `inner_disk_edge`, which is not loaded by name: it is a helper
that `type_I_migration` uses.

---

## How to use it

### The shape of a program

Every REBOUNDx program follows the same five steps:

```rust
let mut sim = reb_simulation_create();   // 1. make a simulation
// ... add particles ...
rebx_attach(&mut sim);                   // 2. attach REBOUNDx to it
let f = rebx_load_force(&mut sim, "gr_potential").unwrap();  // 3. pick an effect
rebx_add_force(&mut sim, f);             // 4. switch it on
// ... set its parameters ...
reb_simulation_integrate(&mut sim, 100.0);  // 5. run
```

### Forces and operators are numbers, not objects

In C, `rebx_load_force` hands you a pointer. Safe Rust cannot pass raw
pointers around, so this port hands you an **index** — an ordinary `usize`
number — and you pass that number back in:

```rust
let gr = rebx_load_force(&mut sim, "gr_potential").expect("known force");
rebx_add_force(&mut sim, gr);
```

Keep that number in a variable. You need it again to set the effect's
parameters.

### Setting parameters

Most effects need you to tell them something — the speed of light, a
timescale, a radius. Each parameter has a name (text) and a type.

You reach the parameters through `rebx_extras_mut`, and you say *what* you are
attaching the parameter to with `rebx_ap`:

- `rebx_ap::force(index)` — a parameter of a force
- `rebx_ap::operator_(index)` — a parameter of an operator
- `rebx_ap::particle(index)` — a parameter of a particle

```rust
if let Some(rebx) = rebx_extras_mut(&mut sim) {
    rebx_set_param_double(rebx, rebx_ap::force(gr), "c", 10065.32);
    rebx_set_param_double(rebx, rebx_ap::particle(1), "tau_a", -1e4);
}
```

The `if let Some(...)` block matters: it must **end** before you touch `sim`
again. That is Rust making sure two parts of your program never change the
same thing at the same time.

Setters and getters, one per type:

| Type | Set | Get |
|---|---|---|
| decimal number | `rebx_set_param_double` | `rebx_get_param_double` |
| whole number | `rebx_set_param_int` | `rebx_get_param_int` |
| unsigned whole number | `rebx_set_param_uint32` | `rebx_get_param_uint32` |
| 3-D vector | `rebx_set_param_vec3d` | `rebx_get_param_vec3d` |
| text | `rebx_set_param_string` | `rebx_get_param_string` |
| a force | `rebx_set_param_force` | `rebx_get_param_force` |

Every getter returns `Option`: `Some(value)` if the parameter was set, `None`
if it was not. Reading a parameter of the wrong type gives `None` rather than
nonsense — in C that same mistake is undefined behaviour.

### `rebx_with`: when you need both at once

A few functions need the simulation **and** the REBOUNDx state changed
together. In C that is easy, because `sim->extras` is a raw pointer. Safe Rust
will not let you hold two changeable references to overlapping data, so this
port gives you `rebx_with`, which lends you both inside a block:

```rust
let l_spin = rebx_with(&mut sim, |sim, rebx| {
    rebx_tools_spin_angular_momentum(sim, rebx)
}).expect("REBOUNDx is attached");
```

You need `rebx_with` for exactly these:

- `rebx_tools_spin_angular_momentum`
- `rebx_tools_spin_energy`
- `rebx_simulation_irotate`
- `rebx_spin_initialize_ode`

Everything else in the ordinary path — attaching, loading, adding, setting
parameters, integrating, saving — does not need it.

### Saving and loading

A whole REBOUNDx configuration (every force, operator, step and parameter) can
be written to a file and read back:

```rust
rebx_output_binary(&mut sim, "state.bin");
```

These files are **interchangeable with the C REBOUNDx** in both directions —
see the verification section below.

### Complete, runnable examples

Four full programs live in `examples/`. Run any of them with:

```bash
cargo run --release --example tides_spin_pseudo
```

| Example | What it shows |
|---|---|
| `tides_spin_pseudo` | A hot Jupiter's tides circularise its orbit and drive its spin to the pseudo-synchronous value (Hut 1981) |
| `tides_spin_kozai` | A Lidov–Kozai cycle with tides, spin and GR, under the adaptive IAS15 integrator |
| `tides_spin_migration` | Two Earth-like planets migrating in a disk with tidal spin evolution, changing a parameter mid-run |
| `rebx_binary_roundtrip` | Writes a configuration to a file, reads it back, and checks all 25 items bit for bit |

---

## Building and testing

```bash
cargo build --release
```

```bash
cargo test --release
```

Expect **137 tests, 0 failures**. (The run also lists two items as
"ignored": those are illustrative code fragments inside documentation
comments, marked so they are shown but not executed. They are not skipped
tests.)

| Test file | Tests | Covers |
|---|---|---|
| `tests/core_params.rs` | 34 | Attach/detach, the parameter store, all six parameter types |
| `tests/forces_conservative.rs` | 32 | `gr`, `gr_full`, `gr_potential`, harmonics — against analytic results |
| `tests/forces_dissipative.rs` | 28 | Migration and damping effects |
| `tests/interpolation_and_spin.rs` | 17 | Cubic-spline parameter interpolation, spin machinery |
| `tests/operators_and_integrators.rs` | 26 | Operators and REBOUNDx's own sub-integrators |

The strict style checker is also clean:

```bash
cargo clippy --release --all-targets
```

It prints nothing. Where a suggestion is deliberately not applied — because
following it would change the order of floating-point arithmetic, and floating
point is **not** associative, so the answer would change — the waiver is
written on its own line at the top of the file with the reason beside it,
rather than hidden in a configuration file.

---

## How this was verified

The method: run the identical experiment in the C REBOUNDx — compiled with
`clang`, Apple's C compiler — and in this Rust code, dump every number as its
raw 64-bit pattern, and compare byte for byte. A test passes only if **every
bit of every number** matches.

Every figure below was measured on one machine: an Apple M5 Max running
macOS Tahoe 26 (26.6.1), with Apple clang 21.0.0 and Rust 1.94.0. (This port
was first made and verified on Windows 11 against Microsoft's C compiler;
where the two platforms' stories differ, both are told below.)

| What was checked | Result on this Mac |
|---|---|
| `tides_spin_pseudo_synchronization`, 10 and 100 orbits (t = 62.83…, t = 628.3…) | **bit-identical** |
| `tides_spin_kozai`, t = 1,000 and t = 100,000 (adaptive IAS15) | **bit-identical** |
| `tides_spin_migration_driven_obliquity_tides`, 10 and 100 orbits | **bit-identical** |
| REBOUNDx binary file written by C, read by Rust | 25/25 items recovered exactly; re-serialising reproduces every byte |
| REBOUNDx binary file written by Rust, read by C | C's 28-line dump identical to C reading its own file |
| Byte comparison of the two written files (10,784 bytes each) | 10,758 of 10,784 bytes identical |
| Automated tests | 137 pass, 0 fail |

Six dynamical comparison runs, six exact matches. Two rows deserve a
sentence each.

**The Kozai result is the strongest.** IAS15 chooses its own step sizes, and
those choices depend on the REBOUNDx forces. Matching bit for bit at
t = 100,000 means both programs took the *identical sequence of thousands of
adaptive steps* — not merely that they arrived at the same place.

**The 26 differing bytes** in the binary file sit at offsets 37–62, inside
one field of the header: the git hash, a stamp recording which source
revision the library was compiled from. The C library writes the text
`notavailable…` there; this port deliberately writes zeros (a documented
deviation — see the differences section below). The field contains no
simulation data and neither library reads it back. Every byte of actual
physics — masses, positions, velocities, parameters, names, types and
orderings — matches. (The Windows-era edition of this document said the
files were 6,392 bytes; that number came from an earlier revision of the
round-trip example. The current source writes 10,784-byte files on both
platforms.)

### Two platforms, two shims — both on the C side

The Rust code is identical on every platform. It is the *C reference* that
needs a small piece of platform help, and each platform needs a different
piece.

**Windows needed a VLA shim; macOS does not.** Two REBOUNDx source files,
`gr_full.c` and `interpolation.c`, use a C feature called variable-length
arrays (a "VLA" is an array whose size is decided while the program runs,
rather than fixed in advance). Microsoft's compiler does not support VLAs, so
the Windows verification carried modified copies of those two files — still
kept in `porttest/msvc_shim/` as the record of that port. Apple's `clang`
supports VLAs, so on macOS both files compile completely unmodified.

**macOS needs a rand_r shim; Windows does not.** Some REBOUND setups place
particles at random using a function called `rand_r` — a recipe that turns a
starting number (the "seed") into a stream of pseudo-random numbers. Upstream
REBOUND carries its own copy of the GNU C library's `rand_r` recipe, but only
switches it on when compiling for Windows. On macOS that switch is off, and
the C build silently gets Apple's `rand_r` — a *different* recipe. Same seed,
different random stream, different initial particles: a bit-for-bit
comparison dies at the very first particle, before any physics has happened.
The fix lives in the sibling crate, at
`rebound_rust/porttest/macos_shim/rand_r_glibc.c`: the same GNU recipe
compiled as its own object file and linked into every C comparison harness,
so the C reference on macOS produces the same stream the identical C code
produces on Linux and Windows. The Rust port implements the GNU recipe on
every platform and needs no shim anywhere.

**And one difference that vanished entirely: the math library.** On Windows,
exactly one math function — `pow`, which raises a number to a power —
disagreed between Microsoft's math library and Rust's, for about 0.03% of
inputs, by at most 2 units in the last place. On this Mac that difference
does not exist: a 200,000-sample sweep of 21 math-library functions (run as
part of the sibling `rebound_rust` verification) found C and Rust
**bit-identical on every function, `pow` included** — both sides resolve
every math call to Apple's math library.

### Re-running the C comparison yourself

You do not need any of this to *use* the crate — it is how the claims above
were established, complete enough to repeat. `~/work` stands in for wherever
you keep the repository.

First, the C reference sources. They are reference material, not part of the
port, so they are cloned at the repository root (the repository's
`.gitignore` already expects them there):

```bash
cd ~/work/rustSolveIt_macos-silicon_SUNDIALS_7_8_0
git clone https://github.com/hannorein/rebound.git rebound/rebound
git -C rebound/rebound checkout dad5f97806ecbb408dcaff728851c64e67f9f6eb  # REBOUND 5.1.1
git clone https://github.com/dtamayo/reboundx.git reboundx
# this port translates REBOUNDx 5.1.0 — check out that release in reboundx/
```

Compile the two C libraries and the rand_r shim:

```bash
cd ~/work/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/rebound/rebound/src
clang -c -DBUILDINGLIBREBOUND -D_GNU_SOURCE -DSERVER \
      -DGITHASH=dad5f97806ecbb408dcaff728851c64e67f9f6eb \
      -O2 -ffp-contract=off *.c
ar rcs librebound_static.a *.o

cd ~/work/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/reboundx/src
clang -c -I../../rebound/rebound/src -I. -D_GNU_SOURCE -DLIBREBOUNDX \
      -O2 -ffp-contract=off *.c
ar rcs libreboundx.a *.o

cd ~/work/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/rebound_rust/porttest/macos_shim
clang -c -O2 rand_r_glibc.c
```

One flag deserves a sentence: `-ffp-contract=off`. Apple-Silicon processors
have a "fused multiply-add" instruction that computes `a*b + c` in one step
with one rounding instead of two. Faster, and slightly *more* accurate — but
"slightly more accurate" means "different in the last bit", and the Rust
compiler never fuses, so the comparison must forbid `clang` from fusing too.
It plays the same role `/fp:precise` played for Microsoft's compiler on
Windows.

Now build one comparison harness and race it against the Rust example:

```bash
cd ~/work/rustSolveIt_macos-silicon_SUNDIALS_7_8_0/reboundx_rust
cargo build --release --examples
cd porttest
clang -I../../rebound/rebound/src -I../../reboundx/src -D_GNU_SOURCE \
      -O2 -ffp-contract=off tides_spin_kozai_c.c \
      ../../rebound_rust/porttest/macos_shim/rand_r_glibc.o \
      ../../reboundx/src/libreboundx.a \
      ../../rebound/rebound/src/librebound_static.a -lm -o tides_spin_kozai_c

./tides_spin_kozai_c 1000                          # writes state_kozai_c.txt
../target/release/examples/tides_spin_kozai 1000   # writes state_kozai_rust.txt

cmp -s state_kozai_c.txt state_kozai_rust.txt && echo bit-identical
shasum -a 256 state_kozai_c.txt state_kozai_rust.txt   # same hash = same bytes
```

`cmp -s` compares two files byte by byte and says nothing unless they differ;
`shasum -a 256` prints each file's SHA-256 fingerprint, so equal fingerprints
mean equal files. Run without the `1000` for the full t = 100,000 default.
The other harnesses in `porttest/` (`tides_spin_pseudo_c.c`,
`tides_spin_migration_c.c`, `rebx_binary_roundtrip_c.c`,
`rebx_binary_read_c.c`) build with the same command shape. Expect a handful
of harmless `clang` warnings about `void**` pointer casts in the harness
files themselves — the libraries compile warning-free.

The full narrative of how the port was made and checked, with every
comparison written out, is in `reboundx_port_test.md` in this folder — a
companion document, not a prerequisite: everything you need to repeat the
comparison is on this page.

---

## What is inside

34 Rust source files, 10,692 lines, translated from 8,452 lines of C.

| File | Lines | What it does |
|---|---|---|
| `core.rs` | 1,507 | Lifecycle, the parameter store, the force/operator dispatch |
| `input.rs` | 1,238 | Reads a REBOUNDx binary file back in |
| `tides_spin.rs` | 643 | Self-consistent spin, tidal and dynamical evolution |
| `output.rs` | 636 | Writes the state to a REBOUNDx binary file |
| `rebxtools.rs` | 475 | Shared helpers: centre of mass, rotations, spin energy |
| `gr_full.rs` | 459 | Post-Newtonian GR for all bodies |
| `gravitational_harmonics.rs` | 454 | J2/J4 — the effect of a squashed central body |
| `yarkovsky_effect.rs` | 418 | Thermal recoil on small asteroids |
| `tides_dynamical.rs` | 409 | Orbital and modal evolution from dynamical tides |
| `stochastic_forces.rs` | 385 | Turbulent forcing from a gas disk |
| `types.rs` | 270 | The REBOUNDx data structures |
| *(23 more)* | | one file per effect, matching the C file names |

Each file keeps the copyright notice of the person who wrote the C original —
Mohamad Ali-Dib, Arya Akmal, Tiger Lu, Stanley A. Baronett, Noah Ferich,
Aleksey Generozov, Kaltrina Kajtazi, Gabriele Pichierri, Donald J. Liveoak,
Phoebe Sandhaus, Pengshuai (Sam) Shi and others.

---

## How this differs from the C, and why each difference is safe

Every difference is mechanical — a consequence of Rust's ownership rules — and
none changes a computed number.

1. **Parameters carry their type.** C stores a `void*` value plus a type tag
   and casts on read; a wrong-type read is undefined behaviour. Rust fuses the
   two into one value, so a wrong-type read returns `None` instead. For correct
   programs the behaviour is identical; for incorrect ones Rust reports rather
   than corrupts.

2. **Lists are vectors, in the same order.** C's `rebx_add_node` *prepends* to
   its linked lists, so a list iterates in reverse insertion order — and that
   order decides the order accelerations are summed, which changes
   floating-point results. This port's lists are vectors whose **index 0 is the
   head**, and the add helper inserts at position 0. Iteration order is
   identical, element for element.

3. **Forces and operators are indices, not pointers** (explained above).

4. **The simulation is passed explicitly.** C reaches the simulation through a
   `rebx->sim` back-pointer. Safe Rust cannot hold that back-pointer, so the
   affected functions take the simulation as their first argument — which is
   why `rebx_with` exists.

5. **The git-hash stamp in binary files** is written as zeros rather than
   filled in from `git` at compile time (explained above).

---

## Credit and citation

**This is a translation. The science, and the credit, belong to the people
below.** If this port is useful to you, cite *them* — a translation earns no
citation of its own.

### REBOUNDx

REBOUNDx is © **Dan Tamayo**, Hanno Rein and the REBOUNDx contributors.
GPL-3.0-or-later.

> REBOUNDx allows you to easily incorporate additional physics into your
> REBOUND simulations. For an overview of the technical details and some
> practical recommendations, see **Tamayo, Rein, Shi and Hernandez 2019**. The
> paper publication lines up with REBOUNDx version 3.0.0.

**Cite for any use of REBOUNDx:** Tamayo, Rein, Shi & Hernandez 2019, MNRAS
491, 2885 — [ADS](https://ui.adsabs.harvard.edu/abs/2020MNRAS.491.2885T) ·
[arXiv:1908.05634](https://arxiv.org/abs/1908.05634).

### REBOUND

REBOUND is © **Hanno Rein**, Shangfei Liu, Dan Tamayo, David S. Spiegel,
Daniel Tamayo, Tiger Lu, Pejvak Javaheri, Rishit Dagli, Dave O'Hallaron,
Ernst Hairer and the REBOUND contributors. GPL-3.0-or-later.

From the REBOUND project:

> If you use this code or parts of this code for results presented in a
> scientific publication, we would greatly appreciate a citation. The simplest
> way to find the citations relevant to the specific setup of your REBOUND
> simulation is:
>
> ```python
> sim = rebound.Simulation()
> # -your setup-
> sim.cite()
> ```

**Cite for any use of REBOUND:** Rein & Liu 2012, A&A 537, A128.

### `sim.cite()` in this port

`sim.cite()` lives in the upstream **Python** package (`pip install rebound`),
which this Rust port does not reimplement. The table below stands in for it:
find the rows matching the effects your simulation uses, and cite those papers.

**At minimum, cite Rein & Liu 2012 for REBOUND, and additionally
Tamayo et al. 2019 if you use REBOUNDx.**

| If you use the effect… | Cite (in addition to Tamayo et al. 2019) |
|---|---|
| `tides_spin` | Lu, Hernandez & Rein 2023, MNRAS 526, 66; Eggleton, Kiseleva & Hut 1998, ApJ 499, 853; Hut 1981, A&A 99, 126 |
| `gr`, `gr_full`, `gr_potential` | Tamayo et al. 2019 §2; Anderson et al. 1975; Nobili & Roxburgh 1986; Newhall, Standish & Williams 1983 |
| `tides_constant_time_lag` | Hut 1981, A&A 99, 126; Bolmont et al. 2015, A&A 583, A116 |
| `radiation_forces` | Burns, Lamy & Soter 1979, Icarus 40, 1 |
| `modify_orbits_forces`, `modify_orbits_direct` | Papaloizou & Larwood 2000, MNRAS 315, 823; Kominami & Ida 2002, Icarus 157, 43 |
| `type_I_migration` | Cresswell & Nelson 2008, A&A 482, 677; Pichierri et al. 2018, CeMDA 130, 54 |
| `exponential_migration` | Hermosillo Ruiz, Lau & Malhotra 2023, MNRAS |
| `gravitational_harmonics` (J2, J4) | standard; see the REBOUNDx documentation |
| `stochastic_forces` | Rein & Papaloizou 2009, A&A 497, 595 |
| `yarkovsky_effect` | Veras, Higuchi & Ida 2019, MNRAS 485, 708 |
| `lense_thirring` | Park et al. 2017, AJ 153, 121 |
| `gas_dynamical_friction` | Ostriker 1999, ApJ 513, 252; Kim & Kim 2007, ApJ 665, 432 |
| `gas_damping_timescale` | Kominami & Ida 2002, Icarus 157, 43 |
| `tides_dynamical` | Press & Teukolsky 1977, ApJ 213, 183 |
| `central_force`, `track_min_distance`, `modify_mass`, `integrate_force`, `interpolation` | Tamayo et al. 2019 |

BibTeX for all of these is collected upstream at
<https://github.com/hannorein/rebound#papers>,
<https://rebound.readthedocs.io/en/latest/citations/> and
<https://reboundx.readthedocs.io>.

### A note on the upstream projects' AI policy

REBOUND's README states:

> REBOUND is a labour of love, created by people. Please refrain from
> submitting issues or pull requests that have been generated by an LLM or
> other fully-automated tools. […] You may of course use AI assistants for your
> own work with REBOUND. Just don't submit any AI generated code.

This translation was produced with AI assistance, and it is **exactly the
"your own work" case that policy permits**: it is a separate derivative work in
its own repository. Accordingly, **nothing from this port has been or will be
submitted to the upstream REBOUND or REBOUNDx repositories** — no issues, no
pull requests. Please do not file bug reports upstream about this port; the
upstream maintainers did not write it and should not be asked to support it.
Report problems here instead, and verify against the original C before
reporting anything upstream.

---

## License

**GPL-3.0-or-later**, the same license as REBOUNDx, because this is a
derivative work of it. The full text is in `LICENSE`.

---

## Related repositories

| Repository | What it is |
|---|---|
| [`once-ere/rebound_rust`](https://github.com/once-ere/rebound_rust) | `rebound_rs` — the pure-Rust REBOUND 5.1.1 this crate requires |
| `rustSolveIt_macos-silicon_SUNDIALS_7_8_0` | **This repository** — the macOS / Apple Silicon edition of the physics simulator whose verification method these ports were built with; carries `rebound_rust/` and `reboundx_rust/` side by side at its top level |
| [`once-ere/rustSolveIt_Win11_SUNDIALS_7_8_0`](https://github.com/once-ere/rustSolveIt_Win11_SUNDIALS_7_8_0) | The Windows 11 edition, where these ports were first made and verified against Microsoft's C compiler |
| [`hannorein/rebound`](https://github.com/hannorein/rebound) | The original REBOUND, in C |
| [`dtamayo/reboundx`](https://github.com/dtamayo/reboundx) | The original REBOUNDx, in C |
