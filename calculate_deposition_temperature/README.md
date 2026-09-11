# Hg/Au(111) Deposition-Temperature Calculation

## Purpose

This directory contains a compact numerical implementation of the
Eichler–Zvara thermochromatographic deposition-temperature relation for the
Hg/Au(111) system.

The present `solve_for_Ta.py` program takes a set of thermodynamic and
thermochromatographic parameters, converts the carrier-gas flow from sccm to
SI units, evaluates the compact Eichler–Zvara equation as a scalar residual,
and locates the first root in a prescribed temperature interval using a
bracketed Brent solver.

The numerical results documented below are taken directly from:

- `solve_for_Ta.py`;
- `input_params.txt`; and
- the three supplied solver output logs for `Q = 0.002`, `0.005`, and
  `0.010 sccm`.

No numerical values in this README are taken from an older version of the
README.

---

## 1. Scope and model definition

The source code implements the following compact deposition-temperature
relation:

```math
\frac{\Delta H_{\mathrm{des}}}{R T_a}
=
\frac{\Delta S_{\mathrm{des}}}{R}
+
\ln\!\left[
\frac{s_0\nu_B\sqrt{2\pi M}\,Q\,g}
{s\,T_a^{3/2}}
\right].
```

<div align="right">(1)</div>

Here:

| Symbol | Meaning | Value/unit convention in the implementation |
|---|---|---|
| $T_a$ | deposition temperature | K |
| $\Delta H_{\mathrm{des}}$ | molar desorption enthalpy | $\mathrm{J\,mol^{-1}}$ |
| $\Delta S_{\mathrm{des}}$ | molar desorption entropy | $\mathrm{J\,mol^{-1}\,K^{-1}}$ |
| $R$ | gas constant | $\mathrm{J\,mol^{-1}\,K^{-1}}$ |
| $s_0$ | standard surface area | $\mathrm{m^2}$ |
| $\nu_B$ | characteristic frequency | $\mathrm{s^{-1}}$ |
| $M$ | molar mass | **entered numerically in $\mathrm{g\,mol^{-1}}$** in Eq. (1), as required by the project implementation |
| $Q$ | carrier-gas volumetric flow | $\mathrm{m^3\,s^{-1}}$ after conversion |
| $g$ | temperature gradient | $\mathrm{K\,m^{-1}}$ |
| $s$ | column cross-sectional area | $\mathrm{m^2}$ |

### Important unit convention for $M$

The compact equation in the source code explicitly uses the numerical molar
mass in **$\mathrm{g\,mol^{-1}}$**. Thus the value

```math
M=200.59\ {\rm g\,mol^{-1}}
```

is inserted into the square-root factor as `200.59`, not as `0.20059`.

This is a convention of the implemented compact equation. Converting $M$ to
k$\mathrm{g\,mol^{-1}}$ inside Eq. (1) without simultaneously re-deriving the
prefactor would change the numerical result.

---

## 2. Thermodynamic sign convention

The input file defines adsorption quantities with the convention

```math
\Delta H_{\mathrm{ads}}<0
```

for exothermic adsorption and

```math
\Delta S_{\mathrm{ads}}<0
```

for the adsorption process.

The solver converts these quantities to desorption quantities by

```math
\Delta H_{\mathrm{des}}=-\Delta H_{\mathrm{ads}},
\qquad
\Delta S_{\mathrm{des}}=-\Delta S_{\mathrm{ads}}.
```

<div align="right">(2)</div>

For the present Hg/Au(111) input:

```math
\Delta H_{\mathrm{ads}}
=
-47.0765871755\ {\rm kJ\,mol^{-1}},
```

```math
\Delta S_{\mathrm{ads}}
=
-104.7019236234\ {\rm J\,mol^{-1}\,K^{-1}},
```

and therefore

```math
\Delta H_{\mathrm{des}}
=
+47.0765871755\ {\rm kJ\,mol^{-1}},
```

```math
\Delta S_{\mathrm{des}}
=
+104.7019236234\ {\rm J\,mol^{-1}\,K^{-1}}.
```

<div align="right">(3)</div>

The thermodynamic quantities are those specified in `input_params.txt` as
Hg/Au(111) thermochemistry evaluated at 298.15 K and 1 bar.

---

## 3. Flow conversion

The source code defines

```math
1\ {\rm sccm}
=
1\ {\rm cm^3\,min^{-1}}
=
\frac{10^{-6}}{60}\ {\rm m^3\,s^{-1}}.
```

<div align="right">(4)</div>

Therefore,

```math
Q_{\mathrm{SI}}
=
Q_{\mathrm{sccm}}\frac{10^{-6}}{60}.
```

<div align="right">(5)</div>

The three supplied calculations use:

| $Q$ (sccm) | $Q$ ($\mathrm{m^3\,s^{-1}}$) |
|---:|---:|
| 0.002 | $3.333333333\times10^{-11}$ |
| 0.005 | $8.333333333\times10^{-11}$ |
| 0.010 | $1.666666667\times10^{-10}$ |

The source code also documents an historical conversion error: `2e-9 m³/s`
corresponds to **0.12 sccm**, not 0.002 sccm. The corrected conversion in
Eq. (5) is used for all results reported here.

---

## 4. Column geometry

The column cross-sectional area is calculated directly from the specified
diameter:

```math
s=\pi\left(\frac{d}{2}\right)^2.
```

<div align="right">(6)</div>

With

```math
d=0.020000\ {\rm m},
```

the program obtains

```math
s=3.141592654\times10^{-4}\ {\rm m^2}.
```

The standard surface area is

```math
s_0=1.000000\times10^{-4}\ {\rm m^2},
```

corresponding to the 1 c$\mathrm{m^2}$ value stated in the input file.

---

## 5. Characteristic Hg frequency

The input file retains two Hg-dominated vibrational frequencies for
traceability:

```math
\tilde\nu_1=19.312\ {\rm cm^{-1}},
\qquad
\tilde\nu_2=19.5168\ {\rm cm^{-1}}.
```

The characteristic frequency used by the compact deposition-temperature
solver is their geometric mean:

```math
\tilde\nu_B
=
\sqrt{\tilde\nu_1\tilde\nu_2}
=
19.414130\ {\rm cm^{-1}},
```

which is supplied to the solver as

```math
\nu_B=5.8202097368\times10^{11}\ {\rm s^{-1}}.
```

<div align="right">(7)</div>

The README does not reinterpret these modes; it records the frequency
definition already specified in `input_params.txt`.

---

## 6. Compact prefactor

For numerical convenience the implementation forms

```math
A=
\frac{s_0\nu_B\sqrt{2\pi M}\,Q\,g}{s}.
```

<div align="right">(8)</div>

Equation (1) can then be written as

```math
f(T)
=
\frac{\Delta H_{\mathrm{des}}}{RT}
-\frac{\Delta S_{\mathrm{des}}}{R}
-\ln\!\left(\frac{A}{T^{3/2}}\right).
```

<div align="right">(9)</div>

The deposition temperature is the numerical root

```math
f(T_a)=0.
```

<div align="right">(10)</div>

For the three supplied flow rates, the resulting prefactors are:

| $Q$ (sccm) | $A$ |
|---:|---:|
| 0.002 | $1.096180534\times10^4$ |
| 0.005 | $2.740451334\times10^4$ |
| 0.010 | $5.480902668\times10^4$ |

Because $A\propto Q$, increasing the flow increases the compact
prefactor and, for the present thermodynamic parameter set, decreases the
computed deposition temperature.

---

## 7. Numerical root-finding

The source code does not rely on a single arbitrary initial guess.

The temperature interval is read from:

```text
T_min_K = 50.0
T_max_K = 1000.0
```

The interval is divided into 2000 equal subintervals, giving 2001 grid
points. The program searches sequentially for the first interval in which
the residual changes sign.

Once a sign-changing interval $[a,b]$ is found, the root is refined using
`scipy.optimize.brentq` with:

```text
xtol = 1e-10 K
rtol = 1e-12
```

The first root encountered is reported.

If no sign-changing interval exists in the requested temperature range, the
program raises a `RuntimeError`.

The reported residuals in the supplied logs are effectively zero:

| $Q$ (sccm) | $T_a$ (K) | $f(T_a)$ |
|---:|---:|---:|
| 0.002 | 444.043288 | $1.126\times10^{-12}$ |
| 0.005 | 410.761265 | $0.000\times10^{0}$ |
| 0.010 | 388.890165 | $-1.332\times10^{-15}$ |

These residuals demonstrate numerical convergence of the implemented
scalar equation; they do **not** establish physical accuracy of the
underlying thermodynamic model.

---

## 8. Deposition-temperature results

### 8.1 Summary of the supplied calculations

| Flow $Q$ (sccm) | Flow $Q$ ($\mathrm{m^3\,s^{-1}}$) | $T_a$ (K) | $T_a$ (°C) | Mobile $\Delta S_{\rm ads}(T_a)$ ($\mathrm{J\,mol^{-1}\,K^{-1}}$) |
|---:|---:|---:|---:|---:|
| 0.002 | $3.333333333\times10^{-11}$ | **444.043288** | **170.893288** | -149.604914 |
| 0.005 | $8.333333333\times10^{-11}$ | **410.761265** | **137.611265** | -149.928803 |
| 0.010 | $1.666666667\times10^{-10}$ | **388.890165** | **115.740165** | -150.156267 |

Thus, within the compact model and the supplied parameter range,

```math
Q\uparrow
\quad\Longrightarrow\quad
T_a\downarrow.
```

For the supplied three-point series, increasing the flow by a factor of five,
from 0.002 to 0.010 sccm, lowers the predicted deposition temperature by

```math
444.043288-388.890165
=
55.153123\ {\rm K}.
```

<div align="right">(11)</div>


---

## 9. Constant-$H,S$ Gibbs-energy reference temperature

The solver additionally reports

```math
T_{\mathrm{H/S}}
=
\frac{\Delta H_{\mathrm{ads}}}{\Delta S_{\mathrm{ads}}}.
```

<div align="right">(12)</div>

For the present input,

```math
T_{\mathrm{H/S}}
=
449.624855\ {\rm K}
=
176.474855^\circ{\rm C}.
```

This quantity is a useful algebraic reference for a simplified
constant-$\Delta H$, constant-$\Delta S$ relation

```math
\Delta G_{\mathrm{ads}}(T)
=
\Delta H_{\mathrm{ads}}
-
T\Delta S_{\mathrm{ads}},
```

but it is **not** the same quantity as the deposition temperature obtained
from Eq. (1). The latter also contains the thermochromatographic transport
and frequency/geometry factor.

The value $449.624855$ K is identical in all three supplied solver logs
because it depends only on the fixed $\Delta H_{\mathrm{ads}}$ and
$\Delta S_{\mathrm{ads}}$, not on $Q$.

---

## 10. `T_cross_K` in the input file

`input_params.txt` also contains

```text
T_cross_K = 450.921106
```

The source program does **not** read this parameter when solving for
`T_a`. It is therefore not the solver's deposition-temperature result and
must not be substituted for the values in Section 8.

The input file labels this value as a comparison value from the same
DFT/phonon calculation. The two distinct temperatures should therefore be
kept separate:

| Quantity | Value | Role |
|---|---:|---|
| Solver $T_a$, at 0.002 sccm | 444.043288 K | compact thermochromatographic deposition root |
| Solver $T_a$, at 0.005 sccm | 410.761265 K | compact thermochromatographic deposition root |
| Solver $T_a$, at 0.010 sccm | 388.890165 K | compact thermochromatographic deposition root |
| Constant-$H,S$ $T_{\mathrm{H/S}}$ | 449.624855 K | algebraic reference |
| Input `T_cross_K` | 450.921106 K | external/comparison value; unused by `solve_for_Ta.py` |

---

## 11. Mobile-adsorption entropy diagnostic

The source code also provides a separate diagnostic function:

```math
\Delta S_{\mathrm{mobile}}
=
R\ln\!\left[
\frac{1}{L_0}
\frac{1}{\nu_B}
\sqrt{\frac{k_BT}{2\pi m}}
\right]
+\frac{R}{2},
```

<div align="right">(13)</div>

where

```math
m=\frac{M}{N_A}
```

is the mass of one Hg particle and the default standard length is

```math
L_0=0.01\ {\rm m}.
```

The implementation converts the input molar mass from $\mathrm{g\,mol^{-1}}$ to kg per
particle before evaluating Eq. (13).

For the three deposition temperatures, the diagnostic gives:

| $Q$ (sccm) | $T_a$ (K) | $\Delta S_{\mathrm{mobile}}(T_a)$ ($\mathrm{J\,mol^{-1}\,K^{-1}}$) |
|---:|---:|---:|
| 0.002 | 444.043288 | -149.604914 |
| 0.005 | 410.761265 | -149.928803 |
| 0.010 | 388.890165 | -150.156267 |

This diagnostic is printed by the program but is **not** inserted back into
the compact residual of Eq. (9). It should therefore be interpreted as a
reported thermodynamic diagnostic rather than an additional term in the
root equation.

---

## 12. Complete active parameter set

The following table collects the parameters actually used by
`solve_for_Ta.py` for the supplied calculations.

| Parameter | Active value | Unit | Role |
|---|---:|---|---|
| $M$ | 200.59 | $\mathrm{g\,mol^{-1}}$ | Hg molar mass |
| $\nu_B$ | $5.8202097368\times10^{11}$ | $\mathrm{s^{-1}}$ | characteristic Hg frequency |
| $s_0$ | $1.0\times10^{-4}$ | $\mathrm{m^2}$ | standard surface area |
| $Q$ | 0.002 / 0.005 / 0.010 | $\mathrm{sccm}$ | carrier-gas flow |
| $g$ | 50.0 | $\mathrm{K\,m^{-1}}$ | temperature gradient |
| $d$ | 0.020 | m | column diameter |
| $s$ | $3.141592654\times10^{-4}$ | $\mathrm{m^2}$ | calculated column area |
| $\Delta H_{\mathrm{ads}}$ | -47.0765871755 | $\mathrm{kJ\,mol^{-1}}$ | adsorption enthalpy |
| $\Delta S_{\mathrm{ads}}$ | -104.7019236234 | $\mathrm{J\,mol^{-1}\,K^{-1}}$ | adsorption entropy |
| $R$ | 8.31446261815324 | $\mathrm{J\,mol^{-1}\,K^{-1}}$ | gas constant |
| $T_{\min}$ | 50 | K | root-search lower bound |
| $T_{\max}$ | 1000 | K | root-search upper bound |
| $L_0$ | 0.01 | m | mobile-entropy diagnostic standard length |
| $N_A$ | $6.02214076\times10^{23}$ | $\mathrm{mol^{-1}}$ | Avogadro constant |
| $k_B$ | $1.380649\times10^{-23}$ | $\mathrm{J\,K^{-1}}$ | Boltzmann constant |

### Parameters retained only for traceability

The input file also contains:

| Parameter | Value | Status in `solve_for_Ta.py` |
|---|---:|---|
| `hg_frequency_1_cm1` | 19.312 $\mathrm{cm^{-1}}$ | traceability only |
| `hg_frequency_2_cm1` | 19.5168 $\mathrm{cm^{-1}}$ | traceability only |
| `T_cross_K` | 450.921106 K | comparison value; **unused by solver** |

---

## 13. Relationship to the Hg/Au(111) atomistic calculation

The thermodynamic input used by the deposition-temperature solver is associated
with the Hg/Au(111) atomistic workflow documented separately in the project.

The supplied computational specification identifies:

| Quantity | Specification |
|---|---|
| Surface | Au(111) |
| Adsorbate | one Hg atom |
| Surface construction | ASE `fcc111` |
| Surface cell | $4\times4$ |
| Slab thickness | 6 Au layers |
| Au atoms | 96 |
| Vacuum | 18 Å |
| Lattice constant | 4.08 Å |
| Adsorption site | fcc |
| Initial Hg height | 3.0 Å |
| Force/energy model | MACE |
| MACE model | `mace-mp-0b3-medium.model` |
| ASE version | 3.29.0 |
| MACE version | 0.3.16 |
| Vibrational method | finite-displacement Hessian |
| Default displacement | 0.005 Å |
| Frequency cutoff | 1 $\mathrm{cm^{-1}}$ |
| Representative thermodynamic pressure | 1 bar in the compact-equation input |
| Thermochemical reference temperature | 298.15 K |

The atomistic calculation is a constrained-coordinate, Gamma-point harmonic
model. It should not be described as a Brillouin-zone-converged phonon
calculation.

For the constrained Hessian workflow, frozen atoms remain present in the
force calculations; only the selected mobile coordinates are displaced and
included in the Hessian coordinate space.

---

## 14. Reproducibility

### Run with the default input file

```bash
python solve_for_Ta.py
```

### Run with an explicit parameter file

```bash
python solve_for_Ta.py input_params.txt
```

The script requires Python together with SciPy because it imports

```python
from scipy.optimize import brentq
```

The calculation is deterministic for a fixed input file and software
environment.

The output explicitly reports the input file used, thermodynamic parameters,
transport parameters, calculated column area, compact prefactor, deposition
temperature, root residual, constant-$H,S$ reference temperature, mobile
entropy diagnostic, and flow-conversion check.

---

## 15. Numerical and physical interpretation

The three supplied calculations show a monotonic decrease in $T_a$ with
increasing flow:

```math
0.002\rightarrow0.005\rightarrow0.010\ {\rm sccm}
```

corresponds to

```math
444.0433\rightarrow410.7613\rightarrow388.8902\ {\rm K}.
```

<div align="right">(14)</div>

This trend is a property of the implemented compact equation under the
specified parameter set. It should not be generalized beyond the model
without checking the assumptions behind the Eichler–Zvara relation.

The calculation should be viewed as a model-based deposition-temperature
prediction. It is not, by itself, an experimental calibration or a proof
that the physical deposition temperature is known to the numerical
precision of the root solver.

In particular, the very small residuals reported in Section 7 measure only
the numerical solution of the implemented equation. They do not quantify
uncertainty in:

- the MACE potential;
- the optimized Hg/Au(111) structure;
- the harmonic approximation;
- the selected Hg frequencies;
- the thermodynamic reference state;
- the compact Eichler–Zvara approximation;
- the flow or temperature-gradient measurements; or
- the representation of the thermochromatographic column.

---

## 16. Important implementation notes

1. **Flow conversion is corrected.**  
   The program uses $1\ {\rm sccm}=10^{-6}/60\ {\rm m^3\,s^{-1}}$. The
   historical `2e-9 m³/s` value corresponds to 0.12 sccm.

2. **The active flow in `input_params.txt` is 0.010 sccm.**  
   The 0.002 and 0.005 sccm cases are represented by the supplied output
   logs, not by the currently uncommented `Q_sccm` line.

3. **The compact equation uses $M$ in $\mathrm{g\,mol^{-1}}$.**  
   This is intentional and source-defined.

4. **Adsorption/desorption signs are explicitly reversed in the solver.**  
   The input stores adsorption values; Eq. (1) uses the corresponding
   desorption quantities.

5. **`T_cross_K` is not used by the solver.**  
   It is retained in the input for comparison/traceability only.

6. **The mobile entropy function is diagnostic only.**  
   Its output is printed after solving but does not alter the deposition root.

7. **The root search uses the first sign-changing interval.**  
   If a future parameter set produces multiple roots inside the requested
   temperature range, this implementation will return the first one found
   on the increasing temperature grid.

8. **The input parser is intentionally simple.**  
   It reads numeric `name=value` lines, ignores blank lines and comment lines,
   and raises an error for a non-numeric value.

9. **Positive-value checks are applied to $M$, $\nu_B$, $s_0$, $Q$,
   $g$, and the calculated column area $s$.**

---

## 17. Compact results for direct reporting

For the present Hg/Au(111) parameter set:

> **$Q=0.002\ {\rm sccm}$: $T_a=444.043288\ {\rm K}=170.893288^\circ{\rm C}$.**

> **$Q=0.005\ {\rm sccm}$: $T_a=410.761265\ {\rm K}=137.611265^\circ{\rm C}$.**

> **$Q=0.010\ {\rm sccm}$: $T_a=388.890165\ {\rm K}=115.740165^\circ{\rm C}$.**

The currently active `input_params.txt` case is therefore:

```math
\boxed{
Q=0.010\ {\rm sccm},
\qquad
T_a=388.890165\ {\rm K}
=
115.740165^\circ{\rm C}
}
```

<div align="right">(15)</div>

with a numerical equation residual of

```math
-1.332\times10^{-15}.
```

---

## 18. Provenance

This README is regenerated from the source and outputs supplied for the
deposition-temperature calculation. The numerical result tables in Sections
7–10 reproduce the values printed by the three attached `solve_for_Ta.py`
logs.

The associated Hg/Au(111) atomistic parameters in Section 13 are included
only to document the thermochemical provenance of the input quantities; they
are not re-computed by `solve_for_Ta.py`.

---

## 19. Summary

The implemented workflow is:

```math
\boxed{
\text{Hg/Au(111) thermochemistry}
\rightarrow
(\Delta H_{\rm ads},\Delta S_{\rm ads},\nu_B)
\rightarrow
\text{flow/geometry conversion}
\rightarrow
\text{compact Eichler--Zvara residual}
\rightarrow
\text{bracketed Brent root}
\rightarrow
T_a
}
```

<div align="right">(16)</div>

Using the supplied parameters, the compact model predicts decreasing
deposition temperature with increasing carrier-gas flow, from

```math
444.043288\ {\rm K}
```

at 0.002 sccm to

```math
388.890165\ {\rm K}
```

at 0.010 sccm.

The README deliberately distinguishes this transport-model deposition
temperature from the constant-$H,S$ Gibbs-energy reference temperature
($449.624855$ K) and from the separate `T_cross_K` value stored in the
input file ($450.921106$ K).
