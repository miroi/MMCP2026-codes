# Eichler–Zvara Deposition-Temperature Calculation

## Overview

This package contains a corrected version of the supplied
`solve_for_Ta.py` and a corrected `input_params.txt`.

Two concrete errors in the original files are corrected:

- the conversion of sccm to m³/s was wrong by a factor of 60;
- the molar-mass convention in the compact equation was implicit, so the
  input is now explicit as g mol⁻¹.

The code retains the compact deposition-temperature equation used by the
project. It does **not** claim that this compact relation is the complete
Eichler–Zvara thermochromatographic transport equation.

## Primary reference

B. Eichler and I. Zvara, “Evaluation of the Enthalpy of Adsorption from
Thermochromatographical Data”, Radiochimica Acta 30(4), 233–238 (1982),
DOI 10.1524/ract.1982.30.4.233.

The paper describes determination of adsorption thermodynamic quantities from
thermochromatographic deposition temperatures and flow-rate dependence, and
also gives a model-based route for estimating adsorption entropy.

The later *Handbook of Nuclear Chemistry*, chapter “Radiochemical Separations
by Thermochromatography”, reproduces the Eichler–Zvara transport formalism
and the mobile-adsorption entropy model.

## Full thermochromatographic model

The handbook gives the general Eichler–Zvara transport relation as its
Eq. 53.2. In schematic form it contains

- a gas-phase transport contribution;
- a surface-residence contribution;
- the carrier-gas velocity under standard conditions;
- the temperature gradient;
- column geometry;
- a standard-state V/A ratio;
- ΔH°_ads and ΔS°_ads;
- the starting temperature T_s;
- the deposition temperature T_a;
- and an integral over the temperature interval.

Thus deposition temperature is an operating-condition-dependent quantity,
not simply the temperature at which ΔG_ads becomes zero.

The full equation requires experimental quantities that are not present in
the supplied `input_params.txt` (notably the starting temperature and
transport/exposure information). Therefore the present solver deliberately
implements the compact algebraic relation that was already used by the
project rather than inventing missing inputs.

## Compact equation implemented

The solver uses

# ΔH_des/(R T_a)
```math
\Delta S_{\mathrm{des}}/R
  +
\ln\left[
    s0 \nu_B \sqrt{2\pi M} Q g
    /
    (s T_a^(3/2))
\right]
```

<p align="right">(1)</p>
with:

`ΔH_des`
```text
positive desorption enthalpy, J mol⁻¹.
```
`ΔS_des`
```text
positive desorption entropy, J mol⁻¹ K⁻¹.
```
`R`
```text
gas constant, J mol⁻¹ K⁻¹.
```
`T_a`
```text
deposition temperature, K.
```
`s0`
```text
standard surface area.
```
`ν_B`
```text
characteristic adsorbent frequency, s⁻¹.
```
`M`
```text
molar mass used by this compact project equation, explicitly supplied
here in g mol⁻¹.
```
`Q`
```text
volumetric carrier-gas flow, m³ s⁻¹.
```
`g`
```text
temperature gradient, K m⁻¹.
```
`s`
```text
column cross-sectional area, m².
```
The solver uses a bracketed Brent root finder instead of `fsolve`. This
makes the solution independent of an arbitrary initial guess and provides
an explicit failure if no root exists in the requested temperature interval.

## Mass convention

The original input stated

```text
M_Hg = 0.20059  # kg/mol
```
and inserted that value directly into `sqrt(2πM)`.

The corrected input instead states

```text
M_adsorbate_g_mol = 200.59
```
and the code documents that this is the convention used by the compact
equation.

This should not be confused with the single-particle mass in the
statistical-mechanical mobile-adsorption model. For that model,

```math
m = M / N_{\mathrm A}
```

<p align="right">(2)</p>
after conversion of M from g mol⁻¹ to kg mol⁻¹.

## Mobile-adsorption entropy

The handbook gives the mobile-adsorption standard entropy for the conventional
V/A = 1 cm standard state as

# ΔS_a°
```math
R \ln\left[
    (1/(1 cm))
    (1/\nu_B)
    \sqrt{\frac{k_{\mathrm B}T_a}{2\pi m}}
]
+ R/2.
```

<p align="right">(3)</p>
Here `m` is the mass of one adsorbate particle, not the molar mass.

The corrected program evaluates this expression only as a diagnostic. It does
not silently replace the supplied DFT/phonon entropy with this model entropy.

For the characteristic substrate frequency, the handbook quotes
ν_B = 5 × 10^12 s⁻¹ for quartz. The project uses the same numerical
frequency as its model input.

## Flow-rate correction

The original file claimed

```text
0.002 sccm = 2e-9 m^3/s.
```
This is incorrect.

By definition,

```text
1 sccm = 1 cm^3/min
       = 1e-6 m^3 / 60 s
       = 1.6666666667e-8 m^3/s.
```
Therefore,

```text
0.002 sccm = 3.3333333333e-11 m^3/s.
```
Conversely,

```text
2e-9 m^3/s = 0.12 sccm.
```
The original calculation therefore used a flow 60 times larger than the
stated 0.002-sccm experiment.

The old output of about 291.3 K is reproducible with the numerical value
`Q = 2e-9 m^3/s`, but that value corresponds to 0.12 sccm.

The corrected input uses

```text
Q_sccm = 0.002
```
and the program performs the conversion internally.

## Corrected input

The corrected parameters are:

```text
M_adsorbate_g_mol = 200.59
nu_phonon_s_inv = 5.0e12
s0_m2 = 1.0e-4
Q_sccm = 0.002
temperature_gradient_K_m = 50.0
column_diameter_m = 0.020
delta_H_ads_J_mol = -44100.0
delta_S_ads_J_mol_K = -121.6
R_J_mol_K = 8.31446261815324
```
The column area is calculated as

```math
s = \pi(d/2)^2
```

<p align="right">(4)</p>
giving

```text
s = 3.1415926536e-4 m²
```
for a 20-mm column.

## Corrected result

With the compact equation and the corrected 0.002-sccm flow, the solver
gives approximately

```math
Q = 3.333333333e-11 m³/s
T_a ≈ 387.6 K
T_a ≈ 114.5 °C.
```

<p align="right">(5)</p>
This is intentionally different from the original 291.3-K result.

The 291.3-K result corresponds to the old numerical flow
`2e-9 m³/s = 0.12 sccm`.

## Thermodynamic ΔG=0 reference

The program separately reports

```math
T_{\mathrm{eq}} = \Delta H_{\mathrm{ads}} / \Delta S_{\mathrm{ads}}.
```

<p align="right">(6)</p>
For

```math
\Delta H_{\mathrm{ads}} = -44.1 kJ mol⁻¹
\Delta S_{\mathrm{ads}} = -121.6 J mol⁻¹ K⁻¹
```

<p align="right">(7)</p>
this gives approximately

```math
T_{\mathrm{eq}} = 362.7 K.
```

<p align="right">(8)</p>
This is only the temperature satisfying

```math
\Delta G_{\mathrm{ads}} = \Delta H_{\mathrm{ads}} - T \Delta S_{\mathrm{ads}} = 0
```

<p align="right">(9)</p>
under the assumption that the supplied H and S are temperature independent.

It is **not** the Eichler–Zvara deposition temperature.

The deposition temperature depends on transport and experimental conditions;
the thermodynamic ΔG=0 temperature does not contain the carrier-flow,
gradient, or column-geometry terms.

## Numerical implementation

The corrected solver:

1. reads a simple `name = value` input file;
2. converts the flow rate using the exact sccm definition;
3. calculates the column cross-sectional area;
4. constructs the compact Eichler–Zvara prefactor;
5. brackets the temperature root over a specified interval;
6. solves it with Brent's method;
7. reports the residual;
8. reports the independent ΔG=0 temperature;
9. evaluates the mobile-adsorption entropy as a diagnostic.

Run with

```text
python solve_for_Ta.py
```
or

```text
python solve_for_Ta.py input_params.txt
```
## Scientific scope and limitations

This correction fixes the coding/input errors identified in the supplied
files. It does not establish that the compact equation is sufficient for a
publication-level Hg/Au(111) thermochromatography prediction.

Before using the result quantitatively, verify:

- the actual experimental carrier-gas flow;
- whether Q is reported at standard conditions;
- the precise V/A standard state;
- compatibility of the supplied DFT/phonon ΔS with that standard state;
- the appropriate characteristic frequency for Au(111);
- whether Hg is a mobile 2-D adsorbate or a localized adsorbate;
- starting temperature and experiment duration;
- whether the full Eichler–Zvara transport integral should be evaluated;
- temperature dependence of ΔH and ΔS;
- coverage/site effects and experimental uncertainties.

## References

Eichler, B.; Zvara, I. (1982).
Evaluation of the Enthalpy of Adsorption from Thermochromatographical Data.
Radiochimica Acta 30(4), 233–238.
DOI: 10.1524/ract.1982.30.4.233

Novgorodov, A. F.; Rösch, F.; Korolev, N. A. (2011).
Radiochemical Separations by Thermochromatography.
Handbook of Nuclear Chemistry, Chapter 53.
DOI: 10.1007/978-1-4419-0720-2_53
