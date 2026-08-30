# Historical Tidal Locking and 3:2 Spin-Orbit Capture of Planet Mercury

**Summary and Provenance Document**

| | |
|---|---|
| **Consolidated by** | Claude (Fable 5), 2026-08-25 |
| **Consolidated from** | 15 file fragments produced by a Gemini agent (inventory in [Part II](#part-ii--provenance-of-this-document)) |
| **Underlying source** | *"Mercury 3:2 Spin-Orbit Resonant Capture Provenance Specification"* (email from Patrick Nash to self, 2026-08-24) |
| **Status** | Duplicates removed; arithmetic verified; three equation errors in the source corrected and documented in [Part IV](#part-iv--verification-and-errata) |

---

## Part I — The Story in Plain English

*(Start here if you are new to the subject. The math comes later.)*

### What is "tidal locking"?

When a small body orbits close to a big one, the big body's gravity pulls harder on the near side of the small body than on its far side. This stretching force is called a **tide**. If the small body is spinning, the tide constantly kneads and flexes it, and that flexing wastes energy as heat — which acts like a brake on the spin. Over millions of years the spin slows down until it settles into a stable, low-friction end state. Our Moon is the most famous example: it got braked all the way down to a **1:1 lock**, spinning exactly once per orbit, which is why it always shows Earth the same face.

### Mercury is locked — but not the way the Moon is

For a long time astronomers assumed Mercury was also 1:1 locked to the Sun. In 1965, radar measurements showed something stranger: Mercury is caught in a **3:2 spin-orbit resonance**. It spins on its axis exactly **three** times for every **two** trips around the Sun.

The consequences are delightfully odd:

| Quantity | Value | Meaning |
|---|---|---|
| Year (orbital period) | **87.969 Earth days** | One trip around the Sun |
| Sidereal day (rotation period) | **58.646 Earth days** | One spin, measured against the stars — exactly ⅔ of the year |
| Solar day (sunrise to sunrise) | **175.938 Earth days** | Exactly **two Mercury years**. On Mercury, a single "day" lasts longer than its year! |

### Why 3:2 and not 1:1? Two key ingredients

1. **An eccentric (oval) orbit.** Mercury's orbit is noticeably stretched (eccentricity e ≈ 0.206, versus Earth's 0.017). Tidal braking pushes a body toward spinning at roughly the speed it moves *at closest approach* — and on an oval orbit that is faster than the average orbital rate. For Mercury's eccentricity, pure tidal friction alone would park the spin at about **1.26×** the orbital rate — meaning the spin must first slow down *through* the 1.5× (3:2) mark on its way there.

2. **A slightly lopsided shape.** Mercury is not a perfect sphere; one axis of its equator is a little longer than the other (a "triaxial" body). The Sun's gravity grabs this permanent bulge like a handle. Whenever the spin rate passes exactly through a whole- or half-number ratio (2:1, 3:2, 1:1, ...), this handle-torque can snap the spin into place and hold it there — like a ball rolling down a bumpy hill that gets caught in a dimple.

**The historical story the model simulates:** Mercury was likely born spinning fast (the model starts it at one rotation every ~12 hours). Over millions of years, solar tides braked the spin. As the spin rate fell and approached 1.5× the orbital rate, the gravitational grip on Mercury's lopsided bulge captured it into the 3:2 resonance — where it has remained ever since.

### What the model actually is

The source document is a **provenance specification**: a precise recipe (equations + constants + numerical-solver settings) for reproducing this capture on a computer. It describes:

- a **two-body system** (Sun + Mercury only),
- Mercury as a spinning, slightly lopsided, deformable body,
- **five numbers** tracked through time (orbit size, orbit shape, orbit position, spin angle, spin rate),
- two competing torques (tidal braking + the "handle" torque on the bulge),
- integrated for **10 million simulated years** with the SUNDIALS CVODE solver,
- with success declared if the simulation ends in the 3:2 state matching the observed periods in the table above.

---

## Part II — Provenance of This Document

### Fragment inventory (all 15 files, deduplicated)

All fifteen fragments turn out to be copies or partial screenshots of **one** underlying document. No fragment contains unique scientific content that is absent from the complete email text.

| # | Fragment | What it is | Unique content? |
|---|---|---|---|
| 1 | `Gmail - Model Provenance_ ... .pdf` | **Complete document** (email, 2026-08-24 08:30). Sections 1–5 in full. | ✅ Primary source |
| 2 | `Gmail - 2 Model Provenance_ ... .pdf` | Same document, emailed again 28 min later (08:59), pasted **twice** in one message | ❌ Duplicate of #1 |
| 3 | `Binder1-6.pdf` | Six browser screenshots of sections 3–5 bound into one PDF | ❌ Subset of #1 |
| 4 | `State Vector Definition.pdf` | Screenshot: §3.1–3.4 | ❌ Subset of #1 |
| 5 | `sys.pdf` | Screenshot: §3.1–3.4 (overlaps #4) | ❌ Subset of #1 |
| 6 | `Orbital Motion.pdf` | Screenshot: §3.2–3.4 | ❌ Subset of #1 |
| 7 | `Triaxial Gravitational Torque.pdf` | Screenshot: §3.3–3.5 | ❌ Subset of #1 |
| 8 | `Secular Tidal Dissipation TorqueA.pdf` | Screenshot: §3.4–3.5 | ❌ Subset of #1 |
| 9 | `Secular Tidal Dissipation Torque.pdf` | Screenshot: §3.5–3.6 | ❌ Subset of #1 |
| 10 | `metadata.txt` | The `[metadata]` TOML block | ❌ Subset of #1 |
| 11 | `constants.txt` | The `[constants]` + `[initial_conditions]` TOML blocks | ❌ Subset of #1 |
| 12 | `solver.txt` | The `[solver]` TOML block | ❌ Subset of #1 |
| 13 | `gemini-code-1787584900782.txt` | Byte-identical copy of `metadata.txt` | ❌ Exact duplicate |
| 14 | `gemini-code-1787584879145.txt` | Byte-identical copy of `constants.txt` | ❌ Exact duplicate |
| 15 | `gemini-code-1787584924422.txt` | Byte-identical copy of `solver.txt` | ❌ Exact duplicate |

### Authorship chain

1. **Physics lineage:** Analytical spin-orbit mechanics — tidal model of **Hut (1981)** (Constant Time Lag), resonance-capture mechanics of **Goldreich & Peale (1966)**; formulation cross-referenced to the **REBOUND / REBOUNDx** open-source frameworks (`spin`, `tidal_dissipation` effects).
2. **Specification:** Composed as *"Mercury 3:2 Spin-Orbit Resonant Capture Provenance Specification"*, targeting the solver build `rustSolveIt_Win11_SUNDIALS_7_8_0` (SUNDIALS 7.8.0 on Windows 11).
3. **Transmission:** Emailed by Patrick Nash to himself twice on 2026-08-24 (08:30 and 08:59).
4. **Fragmentation:** A Gemini agent attempted to consolidate the material and instead emitted the 15 fragments above.
5. **This document:** Consolidation, deduplication, verification, and correction of the fragments (see errata, Part IV).

---

## Part III — The Model Specification (cleaned and corrected)

> Equations marked **[corrected]** differ from the source fragments; the exact discrepancies and the reasons for each correction are documented in Part IV.

### 1. Metadata

```toml
[metadata]
title                 = "Mercury 3:2 Spin-Orbit Resonant Capture Provenance Specification"
system                = "Two-Body Sun-Mercury"
solver_target         = "rustSolveIt_Win11_SUNDIALS_7_8_0"
framework_references  = ["REBOUND (N-body)", "REBOUNDx (spin, tidal_dissipation)"]
authorship_provenance = "Analytical Spin-Orbit Mechanics Formulation"
```

### 2. System Definition and Model Scope

The model simulates the historical rotational deceleration and 3:2 resonance capture of Mercury by integrating a system of ordinary differential equations (ODEs) with SUNDIALS CVODE.

- **Primary body (Sun):** a point mass $M_\odot$.
- **Secondary body (Mercury):** an extended, non-spherical body with mass $m$, mean radius $R$, principal moments of inertia $A < B < C$, tidal Love number $k_2$, and tidal time lag $\tau$.
- **Exclusions (deliberate simplifications):**
  - Obliquity fixed at $\epsilon = 0°$ — Mercury's spin axis stays perpendicular to its orbit (true to within ~2 arcminutes for the real planet).
  - No third bodies (Jupiter, Venus, ...) — their perturbations are excluded to isolate the two-body tidal-capture mechanism. *(Caveat for the curious: in the fuller literature, planetary perturbations that make Mercury's eccentricity wander actually **raise** the probability of 3:2 capture — see Correia & Laskar 2004. This model deliberately studies the clean two-body case.)*

### 3. Physical Constants and Initial Conditions

```toml
[constants]
G                 = 6.67430e-11   # Gravitational constant [m^3 kg^-1 s^-2]
M_sun             = 1.98847e30    # Sun mass [kg]
m_mercury         = 3.3011e23     # Mercury mass [kg]
R_mercury         = 2.4397e6      # Mercury mean radius [m]
C_factor          = 0.34          # Moment of inertia factor C / (m * R^2)
B_minus_A_over_C  = 1.0e-4        # Triaxial asymmetry ratio (B - A) / C
k2                = 0.12          # Secular Love number of degree 2
tau_lag           = 100.0         # Tidal constant time lag [s]

[initial_conditions]
a_0     = 5.790905e10   # Semi-major axis [m] (~0.387098 AU)
e_0     = 0.20563       # Initial orbital eccentricity
M_0     = 0.0           # Initial mean anomaly [rad]
theta_0 = 0.0           # Initial spin angle [rad]
Omega_0 = 1.5e-4        # Initial fast rotation rate [rad/s] (Period ~11.6 hours)  [corrected comment]
```

**What these numbers mean, in words:**

- `C_factor = 0.34`: how Mercury's mass is spread internally (a uniform sphere would be 0.4; Mercury's dense iron core lowers it — MESSENGER measured ≈ 0.346, so 0.34 is realistic). The polar moment of inertia is $C = 0.34\,mR^2 \approx 6.68 \times 10^{35}\ \mathrm{kg\,m^2}$.
- `B_minus_A_over_C = 1e-4`: how lopsided Mercury's equator is — the "handle" the Sun grabs. The real measured value is of this same order ($\approx 1.2\times10^{-4}$).
- `k2 = 0.12`, `tau_lag = 100 s`: how squishy Mercury is and how sluggishly its tidal bulge responds. These two set the *strength* of tidal braking. They are **modeling choices**, not precise measurements: the real capture took far longer than 10 Myr, and dissipation is deliberately set strong here so the whole story fits inside a tractable simulation. (For reference, MESSENGER-era estimates give a present-day $k_2 \approx 0.45$–0.57; the historical effective value is unknown.)
- `a_0`, `e_0`: today's orbit size and shape, used as the starting orbit.
- `Omega_0 = 1.5e-4 rad/s`: the assumed fast primordial spin — one rotation every **≈ 11.6 hours** ($2\pi/\Omega_0 = 41{,}888$ s), about **181×** faster than the orbital rate. *(The source comment said "~1.16 days"; that is arithmetically inconsistent with the given value — see errata E4.)*

### 4. Mathematical Formulation — the Governing ODEs

#### 4.1 State vector

Five quantities are evolved through time:

$$
\mathbf{y}(t) = \begin{bmatrix} a(t) \\ e(t) \\ M(t) \\ \theta(t) \\ \Omega(t) \end{bmatrix}
\in \mathbb{R}^5
$$

| Symbol | Name | Plain meaning |
|---|---|---|
| $a$ | semi-major axis | size of the orbit |
| $e$ | eccentricity | how oval the orbit is (0 = circle) |
| $M$ | mean anomaly | "clock angle" marking where Mercury is along its orbit |
| $\theta$ | rotation angle | which way Mercury's long axis points |
| $\Omega = d\theta/dt$ | spin rate | how fast Mercury rotates |

#### 4.2 Orbital motion

$$
\frac{dM}{dt} = n = \sqrt{\frac{G(M_\odot + m)}{a^3}}
$$

$n$ is the **mean motion** — the average angular speed of the orbit (Kepler's third law). For the given constants, $n = 8.267\times10^{-7}$ rad/s, i.e. one orbit per 87.97 days.

#### 4.3 Triaxial gravitational torque $T_{\mathrm{tri}}$ (the "handle" torque)

The permanent mass asymmetry $(B - A)$ lets the Sun torque Mercury's orientation:

$$
T_{\mathrm{tri}}(f, \theta) = -\frac{3}{2}\, G M_\odot\, \frac{B - A}{r^3}\, \sin\!\big(2(\theta - f)\big)
$$

where the true anomaly $f$ (actual angle of Mercury as seen from the Sun) and distance $r$ follow from the standard Keplerian relations:

$$
r = \frac{a(1 - e^2)}{1 + e \cos f}, \qquad
\tan\frac{f}{2} = \sqrt{\frac{1+e}{1-e}}\, \tan\frac{E}{2}, \qquad
M = E - e \sin E
$$

($E$ is the eccentric anomaly; the last relation is Kepler's equation, solved numerically at each step.) This torque averages to nearly zero — *except* near a resonance, where it becomes the restoring force that locks the spin.

#### 4.4 Secular tidal dissipation torque $\langle T_{\mathrm{tidal}} \rangle$ (the brake)

Using the **Constant Time Lag (CTL)** model (Hut 1981; as implemented in REBOUNDx `tidal_dissipation`), the orbit-averaged tidal torque on Mercury's spin is:

$$
\left\langle T_{\mathrm{tidal}} \right\rangle
= -\,\frac{3\, G M_\odot^2 R^5 k_2 \tau}{a^6}\,
\Big[\, \Omega\, f_1(e) - n\, f_2(e) \,\Big]
$$

with eccentricity functions **[corrected — see erratum E1]**:

$$
f_1(e) = \frac{1 + 3e^2 + \tfrac{3}{8}e^4}{(1 - e^2)^{9/2}}, \qquad
f_2(e) = \frac{1 + \tfrac{15}{2}e^2 + \tfrac{45}{8}e^4 + \tfrac{5}{16}e^6}{(1 - e^2)^6}
$$

**Reading it in words:** the torque brakes the spin whenever $\Omega f_1 > n f_2$ and spins it up otherwise. The braking vanishes at the **pseudo-synchronous rate** $\Omega_{\mathrm{eq}} = n\, f_2(e)/f_1(e)$; at Mercury's $e = 0.20563$ this is $\Omega_{\mathrm{eq}} \approx 1.256\,n$. Because $1.256 < 1.5$, tides alone would drag the spin *down through* the 3:2 mark — which is exactly what gives the triaxial torque its chance to capture.

#### 4.5 Orbital evolution (back-reaction on the orbit)

Tidal dissipation also slowly reshapes the orbit **[corrected — see errata E2, E3]**:

$$
\frac{da}{dt} = \frac{2K}{m\, n\, a}\, \Big[\, \Omega\, f_2(e) - n\, f_3(e) \,\Big]
$$

$$
\frac{de}{dt} = \frac{9\,K e}{m\, n\, a^2}\, \Big[\, \tfrac{11}{18}\,\Omega\, f_4(e) - n\, f_5(e) \,\Big]
$$

where

$$
K = \frac{3\, G M_\odot^2 R^5 k_2 \tau}{a^6}
$$

and the remaining Hut (1981) eccentricity polynomials are:

$$
f_3(e) = \frac{1 + \tfrac{31}{2}e^2 + \tfrac{255}{8}e^4 + \tfrac{185}{16}e^6 + \tfrac{25}{64}e^8}{(1 - e^2)^{15/2}}
$$

$$
f_4(e) = \frac{1 + \tfrac{3}{2}e^2 + \tfrac{1}{8}e^4}{(1 - e^2)^{5}}, \qquad
f_5(e) = \frac{1 + \tfrac{15}{4}e^2 + \tfrac{15}{8}e^4 + \tfrac{5}{64}e^6}{(1 - e^2)^{13/2}}
$$

**Reading it in words:** while the planet still spins fast ($\Omega$ large), spin angular momentum is pumped into the orbit and $a$ grows slightly; eccentricity is (for slow spins) damped. Both effects are tiny compared with the spin evolution — the orbit is nearly a fixed stage on which the spin drama plays out.

#### 4.6 Complete rotational dynamics

$$
\frac{d\theta}{dt} = \Omega, \qquad
\frac{d\Omega}{dt} = \frac{1}{C}\,\Big( T_{\mathrm{tri}}(f, \theta) + \big\langle T_{\mathrm{tidal}} \big\rangle \Big)
$$

Newton's law for rotation: the two torques, divided by the moment of inertia $C = 0.34\,mR^2$, change the spin rate. The steady tidal term slowly bleeds spin; the oscillating triaxial term is what snaps the spin into resonance when the ratio $\Omega/n$ drifts through 1.5.

### 5. SUNDIALS Integrator Configuration (`rustSolveIt_Win11_SUNDIALS_7_8_0`)

```toml
[solver]
method        = "CVODE_BDF"     # Backward Differentiation Formulas for stiff spin-orbit evolution
iteration     = "NEWTON"        # Newton iteration with direct dense linear solver
rel_tolerance = 1.0e-12         # Relative tolerance for angular precision
abs_tolerance = [1.0e-3, 1.0e-6, 1.0e-10, 1.0e-10, 1.0e-14]   # per state variable [a, e, M, theta, Omega]
initial_time  = 0.0             # [s]
final_time    = 3.15576e14      # [s] = 10 million years
max_steps     = 500000000
max_step_size = 864000.0        # [s] = 10 days
```

**Why these choices:** the problem is *stiff* — it mixes fast oscillations (the orbit, ~88 days) with glacial drifts (the spin-down, millions of years) — which is exactly what the implicit BDF method in CVODE is designed for. The five absolute tolerances match the five state variables in order; the spin rate $\Omega$ gets the tightest tolerance ($10^{-14}$) because resonance capture hinges on resolving tiny changes in it.

### 6. Resonant Capture Verification Targets

When the spin ratio reaches the 3:2 attractor ($\Omega/n = 1.5$), the restoring torque $T_{\mathrm{tri}}$ captures the spin into **libration** (a gentle rocking) around the resonance angle $\gamma = 2\theta - 3M$. A successful simulation must reproduce Mercury's observed clock:

| Parameter | Formula | Target value | Verified ✓ |
|---|---|---|---|
| Orbital period $P_{\mathrm{orb}}$ | $2\pi/n$ | ≈ 87.969 Earth days $(7.6005\times10^6\ \mathrm{s})$ | 87.968 d ✓ |
| Sidereal rotation period $P_{\mathrm{rot}}$ | $2\pi/\Omega = \tfrac{2}{3}P_{\mathrm{orb}}$ | ≈ 58.646 Earth days $(5.0670\times10^6\ \mathrm{s})$ | 58.645 d ✓ |
| Solar day $P_{\mathrm{solar}}$ | $\left\lvert \tfrac{1}{P_{\mathrm{rot}}} - \tfrac{1}{P_{\mathrm{orb}}} \right\rvert^{-1} = 2P_{\mathrm{orb}}$ | ≈ 175.938 Earth days $(1.5201\times10^7\ \mathrm{s})$ | 175.935 d = 2·P_orb exactly ✓ |

*("Verified" column: recomputed independently from `G`, `M_sun`, `m_mercury`, `a_0` during consolidation — see Part IV.)*

---

## Part IV — Verification and Errata

### A. Checks that PASSED (source values confirmed)

| # | Check | Result |
|---|---|---|
| V1 | $n = \sqrt{G(M_\odot+m)/a_0^3}$ | $8.2669\times10^{-7}$ rad/s → $P_{\mathrm{orb}} = 87.968$ d ✓ matches target 87.969 d and the real Mercury |
| V2 | $P_{\mathrm{rot}} = \tfrac{2}{3}P_{\mathrm{orb}}$ | 58.645 d ✓ matches target 58.646 d |
| V3 | Solar-day identity $\lvert 1/P_{\mathrm{rot}} - 1/P_{\mathrm{orb}}\rvert^{-1} = 2P_{\mathrm{orb}}$ | Exact algebraically for a 3:2 lock; numerically 175.935 d ✓ |
| V4 | Second-to-day conversions in the targets table | $7.6005\times10^6$ s, $5.0670\times10^6$ s, $1.5201\times10^7$ s ✓ all consistent |
| V5 | `a_0` ↔ "~0.387098 AU" | $5.790905\times10^{10}$ m = 0.3870981 AU ✓ |
| V6 | `final_time` ↔ "~10 Million Years" | $3.15576\times10^{14}/3.15576\times10^{7} = 10^7$ yr ✓ exact |
| V7 | `max_step_size` ↔ "10 days" | $864{,}000$ s = 10.0 d ✓ |
| V8 | Tolerance vector length (5) matches state dimension (5) | ✓ |
| V9 | Physical plausibility of constants | $M_\odot$, $m$, $R$ standard; $C/mR^2=0.34$ ≈ measured 0.346; $(B\!-\!A)/C = 10^{-4}$ ≈ measured $1.2\times10^{-4}$ ✓ |
| V10 | Capture logic self-consistency | Pseudo-synchronous rate $f_2/f_1 = 1.256\,n < 1.5\,n$ at $e=0.20563$ → tidal spin-down necessarily crosses the 3:2 resonance, enabling capture ✓ |

### B. Errata — errors found in the source fragments and CORRECTED above

The following were checked against Hut (1981) — the very paper the source cites — and against dimensional analysis. All four errors appear identically in every fragment that contains the affected passage (so they originate in the source document, not in the fragmentation).

| # | Location | Source fragments said | Corrected to | Why |
|---|---|---|---|---|
| **E1** | §3.4, $f_1(e)$ | denominator $(1-e^2)^{15/2}$ | $(1-e^2)^{9/2}$ | Hut (1981) eq. (11): the spin-torque $\Omega$-term carries $(1-e^2)^{3/2} f_5^{\mathrm{Hut}} / (1-e^2)^6 = f_5^{\mathrm{Hut}}/(1-e^2)^{9/2}$. With $15/2$ the pseudo-synchronous spin rate comes out *below* $n$ at all eccentricities, which would make 3:2 capture impossible — contradicting the document's own §5. |
| **E2** | §3.5, $da/dt$ | $-\dfrac{2aK}{GM_\odot m}\big[\Omega f_2 - n f_3\big]$ | $+\dfrac{2K}{m n a}\big[\Omega f_2 - n f_3\big]$ | (i) **Dimensions:** the source prefactor makes $da/dt$ dimensionless instead of m/s (off by a factor $n a$). (ii) **Sign:** a fast-spinning planet feeds angular momentum *into* its orbit, so $da/dt > 0$ when $\Omega f_2 > n f_3$ (cf. the Moon receding from Earth). The corrected form reduces at $e=0$ to the classical $da/dt = 2K(\Omega - n)/(m n a)$ and satisfies conservation of total angular momentum together with §4.4/§4.6. |
| **E3** | §3.5, $de/dt$ | $-\dfrac{9eK}{GM_\odot m}\big[\Omega f_4 - \tfrac{11}{18} n f_5\big]$ with $f_4 = \frac{1+\frac{15}{4}e^2+\frac{15}{8}e^4+\frac{5}{64}e^6}{(1-e^2)^{13/2}}$, $f_5 = \frac{1+\frac{9}{2}e^2+\frac{5}{8}e^4+\frac{1}{16}e^6}{(1-e^2)^{5}}$ | $+\dfrac{9Ke}{m n a^2}\big[\tfrac{11}{18}\Omega f_4 - n f_5\big]$ with $f_4 = \frac{1+\frac{3}{2}e^2+\frac{1}{8}e^4}{(1-e^2)^{5}}$, $f_5 = \frac{1+\frac{15}{4}e^2+\frac{15}{8}e^4+\frac{5}{64}e^6}{(1-e^2)^{13/2}}$ | The source scrambled Hut (1981) eq. (10): the factor $\tfrac{11}{18}$ belongs on the **spin** ($\Omega$) term, not the orbital one; the two polynomials were attached to the wrong terms; the numerator $1+\frac92 e^2+\frac58 e^4+\frac1{16}e^6$ matches nothing in Hut 1981 (his $f_4 = 1+\frac32 e^2+\frac18 e^4$); the prefactor has the same dimensional problem as E2. The corrected form gives $de/dt < 0$ (eccentricity damping) for slow rotators — the physically required behavior. |
| **E4** | `[initial_conditions]`, comment on `Omega_0` | "Period ~1.16 days" | "Period ~11.6 hours" | $2\pi / (1.5\times10^{-4}\ \mathrm{rad/s}) = 41{,}888$ s $= 11.64$ h $= 0.485$ d. The numeric parameter value was kept authoritative and the comment corrected. (If instead a ~1.16-day period was intended, `Omega_0` would have to be $6.27\times10^{-5}$ rad/s; the simulation outcome is insensitive to this choice — any sufficiently fast initial spin funnels into the same tidal spin-down.) |

*Note on E2/E3 prefactors: they are written for $M_\odot \gg m$ (reduced mass $\beta = mM_\odot/(M_\odot+m) \approx m$), an approximation good to 1 part in $10^7$ here.*

### C. What was NOT changed

- All physical constants, initial conditions, and solver settings — reproduced verbatim (they all passed checks V1–V9).
- The structure of the tidal torque (§4.4), the triaxial torque (§4.3), Kepler relations, state vector, rotational dynamics (§4.6), and all verification targets (§6) — reproduced as sourced; all verified correct.

---

## References

- **Hut, P. (1981)**, "Tidal evolution in close binary systems", *Astronomy & Astrophysics* **99**, 126 — source of the Constant Time Lag tidal model and the eccentricity polynomials $f_1$–$f_5$.
- **Goldreich, P. & Peale, S. (1966)**, "Spin-orbit coupling in the solar system", *The Astronomical Journal* **71**, 425 — classic theory of Mercury's 3:2 resonance capture.
- **Correia, A. C. M. & Laskar, J. (2004)**, "Mercury's capture into the 3/2 spin-orbit resonance as a result of its chaotic dynamics", *Nature* **429**, 848 — modern capture-probability study (context for the two-body simplification).
- **Pettengill, G. H. & Dyce, R. B. (1965)**, "A Radar Determination of the Rotation of the Planet Mercury", *Nature* **206**, 1240 — the discovery that Mercury is not 1:1 locked.
- **REBOUND / REBOUNDx** — open-source N-body framework and effects library (`spin`, `tidal_dissipation`) referenced by the specification. https://rebound.readthedocs.io
- **SUNDIALS 7.8.0 (CVODE)** — LLNL suite of nonlinear/differential-algebraic solvers; the BDF integrator targeted by the specification. https://computing.llnl.gov/projects/sundials
