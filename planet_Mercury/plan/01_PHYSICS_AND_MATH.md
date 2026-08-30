# Planet Mercury Tidal-Locking Plan — Document 1 of 8: The Physics and the Math

**Project:** `planet_Mercury` — a rustSolveIt Jupyter-notebook simulation of how Mercury
became locked to the Sun.
**Audience:** written for a reader with U.S. high-school math and science. Every fact,
number, and equation needed to understand this document is *inside* this document.
**Status:** PLAN (awaiting approval). Nothing has been built yet. All arithmetic in
this document was re-derived by an independent audit agent before presentation.

Two recurring words, explained up front: **ODE** = ordinary differential equation, a
rule for how fast something changes; **secular** = astronomers' word for the slow,
steady, averaged-over-many-orbits part of a motion, as opposed to its fast wiggles.

---

## 1. The story we are going to simulate, in plain English

### 1.1 What is a "tide," and why does it slow a planet's spin?

When a small body (Mercury) orbits close to a big one (the Sun), the Sun's gravity pulls
harder on Mercury's near side than on its far side. That difference in pull stretches
Mercury slightly, raising two small bulges — one facing the Sun, one facing away. This
stretching is called a **tide** (the very same physics that makes Earth's oceans rise
and fall).

If the planet is spinning, the spin keeps dragging the bulges away from the Sun-line,
and the rock has to keep flexing to move the bulges back. Flexing rock wastes energy as
heat, exactly like a bent paperclip warming up. That wasted energy has to come from
somewhere: it comes out of the planet's spin. So the tide acts as a **brake on the
spin**. Over millions of years, the brake slows the spin down until the planet reaches
a stable, low-friction end state. This whole process is called **tidal locking**.

Our Moon is the most famous example: it was braked all the way to a **1:1 lock** — one
spin per orbit — which is why it always shows Earth the same face.

### 1.2 Mercury is locked too — but in a stranger way

For decades astronomers assumed Mercury was 1:1 locked to the Sun, like the Moon is to
Earth. In 1965, radar bounced off Mercury (Pettengill & Dyce) proved something odder:
Mercury spins **exactly three times for every two orbits** — a **3:2 spin-orbit
resonance**. That is Mercury's form of tidal locking, and it produces a wonderfully
weird calendar:

| Quantity | Value | Meaning |
|---|---|---|
| Year (orbital period) | **87.969 Earth days** | one trip around the Sun |
| Sidereal day (rotation period) | **58.646 Earth days** | one spin, measured against the stars — exactly 2/3 of the year |
| Solar day (sunrise to sunrise) | **175.938 Earth days** | exactly **two Mercury years** — on Mercury, one "day" is longer than one "year"! |

### 1.3 Why 3:2 and not 1:1? Two ingredients

1. **An oval orbit.** Mercury's orbit is noticeably stretched: its eccentricity is
   e ≈ 0.206 (0 would be a perfect circle; Earth's is only 0.017). Tidal braking pushes
   a planet toward spinning at roughly the speed it *moves past the Sun at closest
   approach* — and on an oval orbit that is faster than the average orbital speed. For
   Mercury's eccentricity, tides alone would park the spin at about **1.26×** the
   orbital rate. To get down to 1.26×, a fast-spinning planet must first pass **through**
   the 1.5× (= 3:2) mark.

2. **A slightly lopsided shape.** Mercury's equator is not a perfect circle — one axis
   is a tiny bit longer than the other (a "triaxial" body). The Sun's gravity can grab
   this permanent bulge like a **handle**. Whenever the spin rate passes exactly through
   a whole-or-half-number ratio (2:1, 3:2, 1:1, ...), the handle-torque can snap the
   spin into step and hold it — like a ball rolling down a bumpy slope that gets caught
   in a dimple.

**The history the simulation tells:** Mercury was likely born spinning fast (we start it
at one rotation every ≈ 11.6 hours). Solar tides braked the spin for a very long time.
As the spin rate fell toward 1.5× the orbital rate, the Sun's grip on the lopsided bulge
captured Mercury into the 3:2 resonance — where it has stayed ever since.

### 1.4 What "the entire existence of Mercury" means here

Mercury is about 4.5 billion years old. With the specification's tidal numbers, the
braking takes about 4.7 billion years — essentially the planet's whole existence
(Section 8 computes this from scratch). Simulating billions of years step-by-step is
not practical on any computer, so — exactly as the source specification intends — we
run a **time-compressed movie**: we make the tidal friction one thousand times stronger
than the specified value, so the entire multi-billion-year story plays out inside a
10-million-year simulation. The *sequence of events* (fast spin → long braking →
capture at 3:2 → locked forever) is preserved; only the clock is sped up. Section 8.3
defines the exact compression factor and the honest mapping between "movie time" and
"real time."

---

## 2. The model: exactly two bodies

- **The Sun** is treated as a **point mass** M☉ (all its mass concentrated at a point).
  Its own spin and tides raised *on the Sun* are ignored — they are far too small to
  matter for Mercury's story.
- **Mercury** is an **extended, deformable, nearly spherical body**: it has a size
  (radius R), a slightly lopsided permanent shape (its three principal moments of
  inertia obey A < B < C, with C the spin-axis moment — a moment of inertia measures a
  body's resistance to changes in its spin), and it can flex under the Sun's
  tidal pull (described by a "squishiness" number k₂ and a "sluggishness" time lag τ).
- **Deliberate exclusions** (each is a conscious modeling choice, to isolate the clean
  two-body tidal-capture mechanism):
  - **No Jupiter, no Venus, no other planets.** Their tugs are excluded. (In the fuller
    research literature, those tugs make Mercury's eccentricity wander, which actually
    *raises* the chance of 3:2 capture — Correia & Laskar 2004. That richer story is
    reserved for a planned second test.)
  - **No general relativity.** Einstein's famous 43-arcseconds-per-century perihelion
    correction is excluded here, and reserved for the planned second test.
  - **Zero obliquity.** Mercury's spin axis is held exactly perpendicular to its orbit
    plane (true for the real planet to within about 2 arcminutes). This makes the spin
    a single angle instead of a 3-D orientation.

---

## 3. The five numbers the computer tracks (the state vector)

The whole system is described at any instant by five numbers, collected in a vector

y(t) = [ a(t), e(t), M(t), θ(t), Ω(t) ]

| Symbol | Name | Plain meaning |
|---|---|---|
| a | semi-major axis | the size of the orbit (half the long diameter of the ellipse), in meters |
| e | eccentricity | how oval the orbit is (0 = circle, closer to 1 = more stretched) |
| M | mean anomaly | a "clock angle" (radians) that marks where Mercury is along its orbit; it increases at a perfectly steady rate |
| θ | rotation angle | which way Mercury's long axis points (radians) |
| Ω = dθ/dt | spin rate | how fast Mercury rotates (radians per second) |

The computer's job is to advance all five forward through time by solving five coupled
ODEs — equations that say "here is how fast each number is changing right now."

---

## 4. Every constant and starting value, with meanings

```toml
[constants]
G                 = 6.67430e-11   # Gravitational constant [m^3 kg^-1 s^-2]
M_sun             = 1.98847e30    # Sun mass [kg]
m_mercury         = 3.3011e23     # Mercury mass [kg]
R_mercury         = 2.4397e6      # Mercury mean radius [m]
C_factor          = 0.34          # Moment of inertia factor C / (m * R^2)
B_minus_A_over_C  = 1.0e-4        # Triaxial asymmetry ratio (B - A) / C
k2                = 0.12          # Secular Love number of degree 2
tau_lag           = 100.0         # Tidal constant time lag [s]  (see Section 8.3 for the
                                  # time-compressed value 1.0e5 s used in the main runs)

[initial_conditions]
a_0     = 5.790905e10   # Semi-major axis [m] (~0.387098 AU, today's value)
e_0     = 0.20563       # Orbital eccentricity (today's value)
M_0     = 0.0           # Mean anomaly [rad] (start at perihelion)
theta_0 = 0.0           # Spin angle [rad]   (see Section 9.3: this is the phase-sweep knob)
Omega_0 = 1.5e-4        # Fast primordial spin rate [rad/s] (one rotation per ~11.6 hours)
```

**What the less obvious ones mean:**

- **`C_factor = 0.34`** — how Mercury's mass is arranged inside. A uniform solid sphere
  would give 0.4; Mercury's dense iron core concentrates mass at the center and lowers
  it (the MESSENGER spacecraft measured ≈ 0.346, so 0.34 is realistic). The polar moment
  of inertia is
  C = 0.34 · m · R² = 0.34 × 3.3011×10²³ kg × (2.4397×10⁶ m)² ≈ **6.68×10³⁵ kg·m²**.
- **`B_minus_A_over_C = 1e-4`** — how lopsided Mercury's equator is: the "handle" the
  Sun grabs. The real measured value is of the same order (≈ 1.2×10⁻⁴). So
  (B − A) = 10⁻⁴ · C ≈ 6.68×10³¹ kg·m².
- **`k2 = 0.12` and `tau_lag = 100 s`** — how squishy Mercury is, and how sluggishly its
  tidal bulge answers the Sun's pull. Their **product k₂·τ** sets the *strength of the
  tidal brake*. These are **modeling choices**, not precise measurements (present-day
  measurements suggest k₂ ≈ 0.45–0.57, and the value billions of years ago is unknown).
  Section 8 shows that with k₂τ = 12 s the braking takes ~4.7 billion years, and defines
  the documented, boosted value used to compress the story into the simulation window.
- **`Omega_0 = 1.5e-4 rad/s`** — the assumed fast newborn spin: one full turn every
  2π / 1.5×10⁻⁴ = 41,888 s ≈ **11.6 hours**. That is about **181×** faster than the
  orbital rate. The end result is insensitive to this choice — any sufficiently fast
  starting spin funnels into the same braking story.

---

## 5. The governing equations (all of them, with plain-English readings)

### 5.1 Orbital clock (Kepler's third law)

dM/dt = n,  where  **n = √( G(M☉ + m) / a³ )**

n is the **mean motion**: the orbit's average angular speed. With the constants above,

n = √( 6.67430×10⁻¹¹ × (1.98847×10³⁰ + 3.3011×10²³) / (5.790905×10¹⁰)³ )
  = **8.2669×10⁻⁷ rad/s**,

which corresponds to an orbital period P_orb = 2π/n = 7.6004×10⁶ s = **87.97 days** —
the real Mercury year. (Because a changes slowly — Section 5.5 — n and P_orb also drift
slowly; the simulation recomputes n from a at every step.)

### 5.2 Where Mercury actually is on its oval (Kepler's equation)

The mean anomaly M is a steady clock, but on an oval orbit Mercury actually moves faster
near the Sun and slower far away. Two standard conversions recover the true position:

- **Kepler's equation** (solved numerically at every step):  M = E − e·sin E,
  where E is the "eccentric anomaly" (a geometric helper angle).
- **True anomaly** f (the actual Sun-to-Mercury direction angle):
  tan(f/2) = √( (1+e)/(1−e) ) · tan(E/2)
- **Sun–Mercury distance**:  r = a(1 − e²) / (1 + e·cos f)

At perihelion (closest approach, f = 0): r = a(1−e) ≈ 0.794 a. At aphelion (farthest,
f = π): r = a(1+e) ≈ 1.206 a. The r³ in the handle-torque below means the Sun's grip is
about (1.206/0.794)³ ≈ 3.5× stronger at perihelion than aphelion — this perihelion
dominance is precisely why the *3:2* state (which repeats its orientation at perihelion)
is special.

### 5.3 The "handle" torque (triaxial gravitational torque)

T_tri(f, θ) = −(3/2) · G · M☉ · (B − A) / r³ · sin( 2(θ − f) )

**Reading:** the Sun grabs Mercury's long axis and twists it toward the Sun-line. The
twist is strongest when the long axis is 45° off the Sun-line (the sin(2·…) factor), and
much stronger near perihelion (the 1/r³ factor). Averaged over many orbits this torque
nearly cancels — **except** when the spin is at (or crossing) a resonance ratio, where
it becomes the restoring force that locks the spin in step. Its peak size at r = a is

(3/2) × 6.67430×10⁻¹¹ × 1.98847×10³⁰ × 6.68×10³¹ / (5.790905×10¹⁰)³ ≈ 6.8×10¹⁹ N·m.

### 5.4 The tidal brake (secular tidal dissipation torque)

Using the **Constant Time Lag (CTL)** model of Hut (1981) — the same model implemented
by the REBOUNDx framework this specification cross-references — the orbit-averaged
tidal torque on Mercury's spin is:

⟨T_tidal⟩ = −K · [ Ω·f₁(e) − n·f₂(e) ],  with  **K = 3 G M☉² R⁵ k₂ τ / a⁶**

and the eccentricity functions

f₁(e) = (1 + 3e² + (3/8)e⁴) / (1 − e²)^(9/2)
f₂(e) = (1 + (15/2)e² + (45/8)e⁴ + (5/16)e⁶) / (1 − e²)⁶

**Reading:** the brake slows the spin whenever Ω·f₁ > n·f₂ and speeds it up otherwise.
The braking force vanishes at the **pseudo-synchronous rate**

Ω_eq = n · f₂(e)/f₁(e).

At e = 0.20563: f₁ = 1.3695, f₂ = 1.7200, so **Ω_eq ≈ 1.256 n**. Because 1.256 < 1.5,
tides alone would drag the spin *down through* the 3:2 mark toward 1.256 n — which is
exactly what gives the handle torque its chance to capture at 1.5 n. (These f₁, f₂ are
the **corrected** Hut-1981 forms; the uncorrected source had a typo that would have made
capture impossible — full errata in Section 11.)

With the spec constants, K = 3 × 6.67430×10⁻¹¹ × (1.98847×10³⁰)² × (2.4397×10⁶)⁵
× 0.12 × 100 / (5.790905×10¹⁰)⁶ ≈ **2.18×10¹⁹ kg·m²/s** (Section 8.1 uses this).

### 5.5 The orbit slowly answers back

The energy and angular momentum the tide moves around also reshape the orbit, slowly:

da/dt = (2K / (m·n·a)) · [ Ω·f₂(e) − n·f₃(e) ]

de/dt = (9K·e / (m·n·a²)) · [ (11/18)·Ω·f₄(e) − n·f₅(e) ]

with the remaining Hut (1981) eccentricity polynomials

f₃(e) = (1 + (31/2)e² + (255/8)e⁴ + (185/16)e⁶ + (25/64)e⁸) / (1 − e²)^(15/2)
f₄(e) = (1 + (3/2)e² + (1/8)e⁴) / (1 − e²)⁵
f₅(e) = (1 + (15/4)e² + (15/8)e⁴ + (5/64)e⁶) / (1 − e²)^(13/2)

**Reading:** while Mercury still spins fast (Ω large), spin angular momentum is pumped
into the orbit and a grows slightly (the same physics that makes our Moon recede from
Earth); once the spin is slow, the tide damps the eccentricity gently. Both effects are
*tiny* compared with the spin's evolution — the orbit is nearly a fixed stage on which
the spin drama plays out — but we integrate them anyway, because "how do the orbital
period and angular momenta change in time?" is one of the questions this project must
answer and display. (These two equations are the **corrected** forms; the source's
versions had wrong prefactors, a wrong sign, and swapped polynomials — Section 11. The
audit proved, with exact symbolic arithmetic, that the corrected trio of equations
conserves total angular momentum C·Ω + m·n·a²·√(1−e²) **exactly** — the uncorrected
source equations do not.)

### 5.6 Newton's law for the spin

dθ/dt = Ω
dΩ/dt = ( T_tri(f, θ) + ⟨T_tidal⟩ ) / C

The two torques, divided by the moment of inertia C, change the spin rate: the steady
tidal term slowly bleeds spin away; the oscillating handle term is what snaps the spin
into resonance when Ω/n drifts down through 1.5.

That completes the five ODEs: da/dt, de/dt, dM/dt, dθ/dt, dΩ/dt.

---

## 6. Quantities we compute FROM the state, for display and verification

| Quantity | Formula | Why we show it |
|---|---|---|
| Orbital period | P_orb = 2π/n, n = √(G(M☉+m)/a³) | "the year" — asked for explicitly; drifts as a drifts |
| Sidereal rotation period | P_rot = 2π/Ω | "the day" — the locked state must show P_rot = (2/3)·P_orb |
| Spin/orbit ratio | Ω/n | the star of the show: 181 → … → 2.0 → **1.5 forever** |
| Solar day | P_solar = 1 / \|1/P_rot − 1/P_orb\| | the "day longer than a year" punchline; = 2·P_orb when locked |
| Resonance angle | γ = 2θ − 3M | the "is it locked?" dial: γ circulates before capture, then **librates** (rocks gently around a fixed value) after capture |
| Spin angular momentum | L_spin = C·Ω | asked for explicitly ("angular momenta vs time") |
| Orbital angular momentum | L_orb = β·n·a²·√(1−e²), β = m·M☉/(M☉+m) ≈ m | the other half of the ledger |
| Total angular momentum | L_tot = L_spin + L_orb | the conservation audit |
| Spin kinetic energy | E_spin = ½·C·Ω² | energy bookkeeping |
| Orbital energy | E_orb = −G·M☉·m/(2a) | energy bookkeeping |
| Tidal heat released | Q_heat(t) = [E_spin+E_orb](0) − [E_spin+E_orb](t) | where the braking energy went: into warming Mercury |

---

## 7. What "success" looks like (verification targets)

When the spin ratio reaches 1.5 and the handle torque captures it, the resonance angle
γ = 2θ − 3M stops circulating and starts **librating** — rocking gently around a fixed
value, like a settling pendulum. The locked end state must reproduce Mercury's observed
clock:

| Parameter | Formula | Target |
|---|---|---|
| Orbital period P_orb | 2π/n | ≈ 87.969 Earth days (7.6005×10⁶ s) |
| Sidereal rotation period P_rot | 2π/Ω = (2/3)·P_orb | ≈ 58.646 Earth days (5.0670×10⁶ s) |
| Solar day P_solar | 1/\|1/P_rot − 1/P_orb\| = 2·P_orb | ≈ 175.938 Earth days (1.5201×10⁷ s) |

Additional pass criteria the notebook will check automatically:

1. **Capture:** after capture, the running average of Ω/n stays within 10⁻⁴ of 1.5000
   to the end of the run, and γ librates (bounded swing, decaying amplitude).
2. **Periods:** final P_orb, P_rot, P_solar match the table above (P_rot/P_orb = 2/3 to
   at least 5 significant figures).
3. **Angular-momentum ledger:** the drop in C·Ω is matched by the rise in L_orb, in
   TWO eras with different physics (a distinction the build's adversarial review
   sharpened). **Before capture** the handle torque averages to essentially zero over
   the circulating resonance angle, so the model conserves L_tot exactly (secular
   part): budget \|L_tot(t) − L_tot(0)\|/L_tot(0) < **1×10⁻⁹**. **After capture** the
   lock is held by a nonzero average handle torque ⟨T_tri⟩ = −⟨T_tidal⟩ =
   K·n·(1.5f₁ − f₂), and because the specification deliberately gives the handle
   torque no orbital back-reaction, the model itself leaks L_tot at that small
   secular rate (about 2×10⁻¹⁰ per compressed Myr); the locked era is therefore
   checked against the predicted leak, not against zero — and both eras' measured
   values are displayed honestly. (Numerically, the recorded locked-era leak is even
   smaller than predicted, because at lock the orbit's per-step change falls below
   what a 64-bit float can resolve in a — also documented.)
4. **Energy ledger:** E_spin + E_orb only decreases (heat is one-way), and the heat
   curve flattens after capture.

---

## 8. Timescale audit — done from scratch, because it changes the plan

This section is the arithmetic that every later planning decision rests on. One term
first: an **e-folding time** is how long an exponentially decaying quantity takes to
shrink by one factor of the number e ≈ 2.718; "k e-folds" means k of those shrink
steps in a row, a total shrink factor of e^k.

### 8.1 How long does the braking take with the given constants?

Away from the resonances the brake gives (using the linearized secular equation)

dΩ/dt ≈ −(K·f₁/C) · (Ω − Ω_eq),   Ω_eq = 1.256 n,

i.e. exponential decay of Ω **toward the pseudo-synchronous rate Ω_eq** (not toward
zero!) with e-folding time

t_brake = C / (K·f₁) = 6.68×10³⁵ / (2.18×10¹⁹ × 1.3695) ≈ 2.24×10¹⁶ s ≈ **710 million years**.

Getting from Ω₀ = 181.4 n down to the 1.5 n mark therefore needs

ln( (181.4 − 1.256) / (1.5 − 1.256) ) = ln(738) ≈ **6.6 e-folds** ≈ **4.7 billion years**

— essentially Mercury's whole 4.5-billion-year existence, which is astrophysically
realistic! But the specification's simulated window is
final_time = 3.15576×10¹⁴ s = **10 million years**, roughly **470× too short**. Run
literally as written, the spin would only fall from 181.4 n to about 179 n and nothing
else would happen.

> **Finding F1 (spec inconsistency, documented in the provenance):** the source
> specification's constants (k₂τ = 12 s) and its final_time (10 Myr) are mutually
> inconsistent — its own verification checklist verified periods and dimensions but
> never the braking timescale. The specification's *intent* is explicit ("dissipation is
> deliberately set strong so the whole story fits inside a tractable simulation"), so we
> honor the intent with a documented compression factor (Section 8.3) rather than the
> letter.

### 8.2 How expensive is the integration? (why we stage it)

The handle torque wiggles as sin(2(θ − f)). While Mercury spins fast, that wiggle has
period ≈ π/(Ω − n) ≈ π/1.5×10⁻⁴ s ≈ **6 hours**. A solver being asked for 12-digit
accuracy must take several steps per wiggle; over 10 million years that is on the order
of **10¹¹ steps** — hundreds of times more than the specification's own
max_steps = 5×10⁸, and days-to-weeks of compute. The wiggle it would be resolving is
physically irrelevant far from resonance (its average effect is essentially zero).

> **Finding F2 (spec infeasibility, documented in the provenance):** the specification's
> tolerance settings + its oscillatory torque + its 10-Myr window imply a step count
> that exceeds its own max_steps by ~2 orders of magnitude. The standard, physically
> justified fix (used throughout the research literature) is **staging**:
>
> - **Stage S (secular, far from resonance, Ω/n > 2.2):** integrate the smooth system
>   with the handle torque switched off (it averages to zero there). The solver then
>   takes giant steps; this stage costs almost nothing.
> - **Stage R (resonant, Ω/n ≤ 2.2):** integrate the full five-equation system with the
>   handle torque on, covering the 2:1 crossing, the 3:2 crossing, capture, and
>   libration. Feasible: ≈ 10⁸ steps, minutes of compute.
> - **Validation of the seam:** over a test window at Ω/n = 3, run both models and show
>   their secular drifts agree (a notebook cell + an automated test).
>
> Every stage still uses SUNDIALS CVODE — BDF (backward differentiation formulas) with
> Newton iteration (repeated refine-the-guess corrector steps) and a dense linear
> solver — at the specified tolerances. No hand-rolled integrator anywhere.

### 8.3 The time-compression factor

To fit the multi-billion-year story into the 10-Myr window, the main ("movie") run
multiplies the brake strength k₂τ by a documented factor **S = 1000**
(equivalently τ_lag: 100 s → 1.0×10⁵ s, keeping k₂ = 0.12):

t_brake(movie) = 710 Myr / 1000 ≈ **0.71 Myr**, so the milestones land at
6.6 × 0.71 ≈ **4.7 Myr** for the 3:2 crossing (with the stage handover at Ω/n = 2.2
around 3.7 Myr and the 2:1 crossing around 3.9 Myr) — leaving **more than 5 million
years** of window to display the locked state after capture. (The audit corrected an
earlier draft here: with the previously proposed S = 500, capture would have landed at
≈ 9.4 Myr — dangerously close to the end of the window.)

**Honest mapping:** the secular equations are exactly linear in k₂τ, so 1 year of movie
time ≈ 1000 years of real history for the braking part (the capture moment itself is a
fast event in both). Both the movie run and a spec-literal run (Section 10) are kept,
labeled, and stored in the database.

**Adiabaticity check (the compression does not break the physics):** capture is only
possible if the brake drifts the spin slowly compared with the resonance's own rocking
motion. The libration (rocking) frequency at 3:2 is

ω_lib = n·√( 3·(B−A)/C · H(e) ),  H(e) = (7/2)e − (123/16)e³ ≈ 0.653 at e = 0.20563,
so ω_lib ≈ 0.014·n ≈ 1.16×10⁻⁸ rad/s (a rocking period of ≈ 17 years).

The compressed brake drifts Ω at ≈ 9.0×10⁻²¹ rad/s² near the resonance, while the
adiabatic ceiling is ω_lib² ≈ 1.3×10⁻¹⁶ rad/s² — a safety factor of ≈ **15,000**
(over 4 orders of magnitude). The compression is safe. (The capture *probability* of
Section 9 is provably independent of S — the brake strength cancels out of the
formula — so the compression does not bias the coin-flip either.)

---

## 9. Resonance capture is a game of chance — and the plan embraces that

### 9.1 The crossings, in order

Starting at Ω = 181 n and braking downward, the spin crosses integer-and-half ratios in
order: … 3:1, 5:2, **2:1**, … then **3:2**. At each crossing the handle torque gets one
chance to capture. If every capture fails, the spin settles at the pseudo-synchronous
rate 1.256 n and stays there (it can never reach 1:1, because the brake turns off at
1.256 n).

### 9.2 The odds

Classic resonance-capture theory (Goldreich & Peale 1966) gives the capture probability
at 3:2 as

P ≈ 2 / (1 + π·V/(2W))

where V is the mean brake torque at the resonance and W is the part of the brake that
rocks back and forth with the libration. For this model's constants (spec-literal
strength; both V and W scale with the brake, so the ratio — and P — is the same at any
compression):

V = K·n·(1.5·f₁ − f₂) = 2.18×10¹⁹ × 8.2669×10⁻⁷ × 0.3343 ≈ **6.0×10¹² N·m**
W = K·f₁·ω_lib = 2.18×10¹⁹ × 1.3695 × 1.157×10⁻⁸ ≈ **3.5×10¹¹ N·m**
π·V/(2W) ≈ 27.4  →  P = 2/28.4 ≈ **7.0%** per crossing.

(At 2:1 the odds are smaller still — about 1%.) This is not a defect — it is the real
scientific puzzle of Mercury: the famous 1966 paper computed essentially this number,
and later work (Correia & Laskar 2004) showed how a wandering eccentricity raises the
lifetime odds. A single deterministic run with an arbitrary starting phase will
*usually sail through* the resonance.

### 9.3 The plan's deterministic answer: a documented phase sweep

Whether one particular run captures depends deterministically on the spin **phase** θ at
the moment of crossing — effectively a coin whose outcome is fixed by the starting angle
θ₀. So the build will:

1. Integrate once to just above the 3:2 crossing — the restart state is saved **when
   Ω/n first falls to 1.6** — and store that state in the database.
2. From the restart point, run a **sweep of Nₚ = 64 short branch runs** that differ only
   in a tiny, documented phase offset added to θ. Each branch is cheap (it only covers
   the crossing window). Every branch's outcome (captured / passed through) is stored.
3. The **measured capture fraction** is displayed next to the ≈ 7% theoretical estimate,
   with honest statistical error bars (64 coin flips at 7% odds typically give 2–8
   captures; there is about a 1% chance of zero captures, in which case one documented
   finer-grid re-sweep is run) — turning a nuisance into the notebook's best science
   lesson.
4. The first captured branch becomes the **canonical trajectory**: it is continued to
   the end of the window, and its final periods are checked against Section 7's targets.
   Its exact θ offset is recorded in the database and in the notebook, so the whole
   pipeline is reproducible bit-for-bit.

### 9.4 A guaranteed-capture encore (exploration, clearly labeled)

One extra display run raises the eccentricity to e = 0.285, where the pseudo-synchronous
rate f₂/f₁ equals exactly 1.5 — there the brake *delivers* the spin to the 3:2 doorstep
and capture is certain. This "what if Mercury's orbit were more oval?" run is clearly
labeled an exploration beyond the specification, and it beautifully demonstrates *why*
eccentricity is the secret ingredient of Mercury's lock.

---

## 10. The runs the project will produce (summary)

| Run | Purpose | Brake k₂τ | Window | Expected outcome |
|---|---|---|---|---|
| A. Spec-literal | fidelity to the source spec, run exactly as written | 12 s | 10 Myr (secular stage only — Section 8.2 makes the full literal run infeasible; a short full-model segment proves the equations) | spin falls 181.4 n → ≈ 179 n; no capture (documented Finding F1) |
| B. Movie | the braking story: fast spin → 2:1 pass → restart saved | 12,000 s (S = 1000) | 0 → restart at Ω/n = 1.6 (≈ 4.6 Myr) | restart state saved just above the 3:2 crossing |
| C. Phase sweep (64 branches) | measure the capture odds | 12,000 s | crossing window only, from the restart | capture fraction ≈ theory ≈ 7% (with stated error bars) |
| B-final. Canonical history | the first captured branch, continued | 12,000 s | restart → 10 Myr | **capture at 3:2 ≈ 4.7 Myr**, then locked; final periods = Section 7 targets |
| D. High-e encore | why eccentricity matters | 12,000 s | crossing window | guaranteed capture |
| E. Seam validation | secular-vs-full agreement at Ω/n = 3 | 12,000 s | short window | drift rates agree to tolerance |

All runs use SUNDIALS CVODE, BDF method (backward differentiation formulas), Newton
iteration, dense linear solver, relative tolerance 1.0×10⁻¹², absolute tolerances
[1.0×10⁻³ (a), 1.0×10⁻⁶ (e), 1.0×10⁻¹⁰ (M), 1.0×10⁻¹⁰ (θ), 1.0×10⁻¹⁴ (Ω)], max step
864,000 s (10 days), from the initial conditions of Section 4 (per-run overrides listed
in the table). The solver plumbing (exact SUNDIALS calls, output cadence, restart
mechanism) is specified in plan document 2, which restates all of these settings in
full.

---

## 11. Errata inherited from the source specification (already corrected)

The source specification (an email-derived provenance document) contained four errors,
found by checking against Hut (1981) — the very paper it cites — and by dimensional
analysis. This plan uses the corrected forms throughout (they appear in Section 5):

- **E1 — f₁ exponent:** the source wrote f₁'s denominator as (1−e²)^(15/2); correct is
  **(1−e²)^(9/2)**. With the wrong exponent, the pseudo-synchronous rate would sit
  *below* n at all eccentricities and 3:2 capture would be impossible — contradicting
  the document's own goal.
- **E2 — da/dt prefactor and sign:** the source wrote da/dt = −(2aK/(G·M☉·m))·[Ωf₂ − nf₃],
  which is dimensionally wrong (it does not even yield meters-per-second) and has the
  wrong sign. Correct: **da/dt = +(2K/(m·n·a))·[Ωf₂ − nf₃]** — a fast-spinning planet
  pushes its orbit *outward* (like the receding Moon), and at e = 0 this reduces to the
  classical da/dt = 2K(Ω−n)/(m·n·a).
- **E3 — de/dt scrambled:** the source misplaced the 11/18 factor (it belongs on the Ω
  term), attached the wrong polynomials to the wrong terms, used a polynomial that
  matches nothing in Hut (1981), and repeated E2's dimensional error. Correct:
  **de/dt = +(9Ke/(m·n·a²))·[(11/18)·Ω·f₄ − n·f₅]** with f₄, f₅ as in Section 5.5;
  this gives eccentricity *damping* for slow rotators, as physics requires.
- **E4 — a comment typo:** the source described Ω₀ = 1.5×10⁻⁴ rad/s as a "~1.16 day"
  period; the correct period is 2π/Ω₀ = 41,888 s ≈ **11.6 hours**. The numeric value is
  authoritative; the outcome is insensitive to it either way.

(The prefactors in E2/E3 use the approximation M☉ ≫ m — good here to 1 part in 10⁷.
The audit verified, symbolically and numerically, that the corrected trio E2/E3 + the
tidal torque conserves total angular momentum exactly.)

To these, this plan adds its own two findings — **F1** (braking timescale vs. window,
Section 8.1) and **F2** (step-count infeasibility, Section 8.2) — with their documented
resolutions. All findings also appear in the project's provenance document.

---

## 12. References (for the curious reader)

- Hut, P. (1981), "Tidal evolution in close binary systems," *Astronomy & Astrophysics*
  99, 126 — the Constant Time Lag tidal model and the f₁–f₅ polynomials.
- Goldreich, P. & Peale, S. (1966), "Spin-orbit coupling in the solar system," *The
  Astronomical Journal* 71, 425 — the theory of Mercury's 3:2 capture and its odds.
- Correia, A. C. M. & Laskar, J. (2004), "Mercury's capture into the 3/2 spin-orbit
  resonance as a result of its chaotic dynamics," *Nature* 429, 848.
- Pettengill, G. H. & Dyce, R. B. (1965), "A Radar Determination of the Rotation of the
  Planet Mercury," *Nature* 206, 1240.
- REBOUND / REBOUNDx — open-source N-body framework and effects library referenced by
  the specification (this workspace contains complete pure-Rust ports of both).
- SUNDIALS 7.8.0 (CVODE) — the ODE solver suite; this workspace's vendored pure-Rust
  port is the only integration backend used.
