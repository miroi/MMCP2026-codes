# Hg adsorption on Au(111): Applied Theory and Code Description

## Purpose

This document describes the physical model, equations, numerical
procedures, configuration controls, file outputs, and program
architecture implemented in the modular Hg-on-Au(111) workflow supplied
with this project. It is intended as a stand-alone
Supplementary-Information-style description of the calculation rather
than as a line-by-line programming manual.

## Equation numbering and cross-references

All displayed mathematical equations are assigned a stable reference
label. Within an RST/Sphinx document, an equation can be referenced with
`` :eq:`eq-001 ``<span class="title-ref"> when the
</span><span class="title-ref">sphinx.ext.mathjax</span>\` extension is
enabled. The explicit labels are intended to remain stable even if
surrounding prose is edited.

The workflow combines:

- construction of a periodic Au(111) slab with ASE;
- placement of one Hg atom at a selected adsorption site;
- constrained geometry optimization using a MACE interatomic potential;
- Gamma-point finite-displacement vibrational analysis;
- two selectable vibrational backends, `custom` and `ase`;
- two selectable thermochemistry backends, `custom` and `ase`;
- monatomic ideal-gas thermodynamics for the gas-phase Hg reference;
- temperature-dependent adsorption thermodynamics;
- a temperature scan; and
- a robust bisection search for the temperature at which
  $`\Delta G_{\text{ads}}(T) = 0`$.

The present implementation is a constrained-coordinate, Gamma-point
harmonic model. It is not a Brillouin-zone-converged phonon calculation.

## 0. Complete computational specification

The following parameter set is the explicit computational specification
documented in the supplied supplementary-information files. It should be
treated as the reference setup when reproducing the representative
numerical results unless a different configuration is stated explicitly.

- surface: Au(111);
- adsorbate: one Hg atom;
- surface construction: ASE `fcc111`;
- surface cell: 4 x 4;
- slab thickness: 6 Au layers;
- vacuum thickness: 18 Angstrom;
- lattice constant: 4.08 Angstrom;
- adsorption site: fcc;
- initial Hg height: 3.0 Angstrom;
- electronic-structure/force model: MACE;
- MACE model: `mace-mp-0b3-medium.model`;
- ASE version: 3.29.0;
- MACE version: 0.3.16;
- phonon method: finite-displacement Hessian;
- default finite-displacement amplitude: 0.005 Angstrom;
- thermodynamic frequency cutoff: 1 cm<sup>-1</sup>;
- representative Hg pressure: 1e-5 bar;
- temperature scan: 50--600 K;
- representative scan temperatures: 50, 100, 200, 298.15, 400, 500, and
  600 K.

The project configuration is expressed in bar for pressure, Angstrom for
geometric displacements, eV for energies, cm<sup>-1</sup> for
vibrational frequencies, and K for temperature. Conversion to SI
pressure is performed internally where required.

The representative numerical values in this document are therefore
model-specific. They are not intended to establish convergence of the
physical problem without the convergence studies described later.

## 1. Physical system and structural model

### 1.1 Au(111) slab

The clean surface is generated with ASE's `fcc111` construction. The
standard configuration supplied with the project uses

- surface cell: `4 x 4`;
- slab thickness: 6 Au layers;
- lattice constant: 4.08 Angstrom;
- vacuum thickness: 18.0 Angstrom.

The resulting slab contains 96 Au atoms, i.e. 16 Au atoms per layer for
six layers. Periodic boundary conditions are inherited from the ASE
surface construction and are retained during optimization and force
evaluation.

The bottom Au layers can be constrained during geometry relaxation
through `n_frozen_layers`. In the standard setup this value is 2.

### 1.2 Hg adsorption geometry

One Hg atom is added to the optimized clean slab. The adsorption site is
selected with `adsorption_site`; the supplied configuration uses the fcc
site. `hg_height` specifies the initial vertical height used when
placing Hg. The standard initial height is 3.0 Angstrom.

The Hg/Au system therefore contains 97 atoms before vibrational
coordinate selection: 96 Au atoms plus one Hg atom.

### 1.3 Separation of relaxation and vibrational constraints

A central feature of the implementation is that the constraints used for
geometry relaxation are not automatically the constraints used to define
the vibrational coordinate space.

`n_frozen_layers` controls the geometry optimization. A separate
parameter, `n_frozen_layers_vibrations`, controls which bottom Au layers
are excluded from the vibrational coordinate space. This second setting
does not alter the optimized structure.

The vibrational freezing can therefore be increased after optimization.
It may be equal to the total number of Au layers. Hg is not frozen by
this Au layer setting.

For example, with a six-layer slab and `n_frozen_layers_vibrations = 6`:

- all 96 Au atoms are frozen in the vibrational coordinate space for
  clean Au;
- no Au vibrational coordinates remain for clean Au;
- the single Hg atom remains mobile for Hg/Au;
- the Hg/Au vibrational Hessian consequently contains 3 coordinates and
  3 modes.

The frozen atoms are not deleted from the physical system. They remain
in all force evaluations. Only the selected mobile coordinates are
displaced and included in the Hessian.

## 2. Electronic-energy model and adsorption energy

The atomic energies are evaluated with a MACE calculator. The default
model is `mace-mp-0b3-medium.model`. The calculator is configured for
CUDA when a CUDA-capable PyTorch installation is available; otherwise
CPU execution is used.

The adsorption energy is defined as

<span id="eq-008">
``` math
\Delta E_{\text{ads}} = E_{\text{HgAu}} - E_{\text{Au}} - E_{\text{Hg}}
```
</span>

where

- $`E_{\text{HgAu}}`$ is the potential energy of the optimized Hg/Au
  slab;
- $`E_{\text{Au}}`$ is the potential energy of the optimized clean Au
  slab; and
- $`E_{\text{Hg}}`$ is the potential energy of an isolated Hg atom
  evaluated with the same MACE calculator.

Negative $`\Delta E_{\text{ads}}`$ therefore denotes energetically
favorable adsorption at the electronic-energy level.

For the representative calculation documented in the project outputs:

    E_Au    = -304.8265897028 eV
    E_HgAu  = -305.4524593272 eV
    E_Hg    =   -0.1245467400 eV

and hence $`\Delta E_{\text{ads}} = -0.5013228844`$ eV. Using
96.48533212 kJ mol<sup>-1</sup> per eV gives approximately
$`\Delta E_{\text{ads}} = -48.3703`$ kJ mol<sup>-1</sup>.

The electronic adsorption energy is kept separate from vibrational and
gas thermodynamic corrections throughout the implementation.

## 3. Geometry optimization

The clean Au slab and Hg/Au system are optimized independently. The
optimizer is configurable; the supplied configuration uses FIRE.

The force convergence thresholds are independently configurable:

- `fmax_au` for the clean slab;
- `fmax_hg` for the Hg/Au system.

The geometry constraints are applied to the specified bottom Au layers.
The reported maximum force is evaluated on the mobile atoms.

The optimization stage produces optimized atomic structures and optional
optimizer trajectories. The final potential energies are then used in
the adsorption-energy expression above.

The Hg--surface distance is reported before and after optimization,
together with its change.

## 4. Harmonic vibrational theory

### 4.1 Constrained Hessian

Let the mobile vibrational Cartesian coordinates be collected into a
vector $`q`$ of dimension $`3N_{\text{mobile}}`$. The harmonic potential
around the optimized structure is

<span id="eq-001">
``` math
V(q) = V_0 + \tfrac{1}{2} q^{T} H q
```
</span>

where $`H`$ is the mass-unweighted Cartesian Hessian in the selected
coordinate space.

The key constrained-system construction is:

- all atoms remain present in the force calculation;
- only mobile atoms are displaced;
- forces on the complete system are calculated after every displacement;
  and
- only the Hessian block associated with mobile coordinates is
  diagonalized.

This preserves the interaction of mobile atoms with the vibrationally
frozen atoms.

### 4.2 Finite differences

The custom backend constructs the Hessian from finite differences of
forces. For a central two-point scheme, the derivative is represented
schematically by

<span id="eq-002">
``` math
H_{ij} = -\frac{F_i(q_j + d) - F_i(q_j - d)}{2d}
```
</span>

where $`d`$ is the configured displacement in Angstrom. The minus sign
appears because forces are the negative gradient of the potential.

The implementation also supports the ASE finite-displacement route. For
the ASE backend, ASE's `Vibrations` machinery is used with the
explicitly selected mobile atom indices. The full `Atoms` object and its
calculator remain available to the force calculations.

The displacement is controlled by `phonon_displacement = 0.005` in
Angstrom. For ASE, `ase_vibration_nfree` controls the finite-difference
scheme exposed by ASE.

For a four-point central finite-difference force derivative, when the
selected implementation uses `nfree = 4`, the corresponding first
derivative can be formed from the forces at $`-2d`$, $`-d`$, $`+d`$, and
$`+2d`$. With the sign convention $`F = -dV/dq`$, the Hessian is
obtained from the negative force derivative. The two-point central
formula remains

<span id="eq-003">
``` math
H_{ij} = -\frac{F_i(q_j+d) - F_i(q_j-d)}{2d}.
```
</span>

For a structure with $`N_{\text{mobile}}`$ mobile atoms, the central
two-point Hessian requires $`2 \times (3 N_{\text{mobile}})`$ complete
force evaluations, because each of the $`3 N_{\text{mobile}}`$ Cartesian
coordinates is displaced in both directions.

For the default 4 x 4 x 6 slab:

- clean Au: 64 mobile atoms -\> 192 Hessian coordinates -\> 384 force
  evaluations;
- Hg/Au: 65 mobile atoms -\> 195 Hessian coordinates -\> 390 force
  evaluations.

These counts refer to force evaluations on the full physical systems.
The frozen atoms are retained in every such calculation.

### 4.3 Mass weighting and normal modes

After construction of the Cartesian Hessian, the Hessian is mass
weighted. For Cartesian coordinates $`i`$ and $`j`$ associated with
atoms $`a`$ and $`b`$,

<span id="eq-004">
``` math
D_{ij} = \frac{H_{ij}}{\sqrt{m_a m_b}}
```
</span>

with the appropriate atomic masses in atomic mass units.

Diagonalization gives eigenvalues $`\lambda_k`$ and normal-mode
eigenvectors. The vibrational angular frequencies are related to the
eigenvalues by

<span id="eq-005">
``` math
\omega_k = \sqrt{\lambda_k}
```
</span>

for positive $`\lambda_k`$. Frequencies are converted to wavenumbers in
cm<sup>-1</sup>.

The code retains signed frequencies for diagnostics, so imaginary modes
can be identified explicitly.

### 4.4 Imaginary-mode handling

The configured `frequency_cutoff_cm` determines which positive
frequencies are retained for thermochemistry. The supplied default is
1.0 cm<sup>-1</sup>.

Significant imaginary modes are not silently converted into
thermodynamic modes. With `fail_on_imaginary = true` the calculation
stops when a significant imaginary instability is detected. This is
preferable to masking an unstable optimized structure by simply removing
negative frequencies.

### 4.5 Acoustic modes and constrained slabs

For a freely translating isolated system, translational acoustic modes
can be expected near zero frequency. Here, however, bottom Au atoms are
constrained and the vibrational problem is a constrained-coordinate
problem. Consequently a global translation of all mobile atoms is not
necessarily a zero mode.

The implementation therefore does not apply an acoustic sum-rule
correction by default. This is particularly important when the bottom
layers are frozen.

## 5. Hg mode projection

For Hg/Au, the code identifies which normal modes are dominated by Hg
motion. The Hg atom index in the mobile-coordinate space is determined
and the mode projection onto the Hg Cartesian subspace is evaluated.

A mode is classified as Hg-dominated when its Hg projection is greater
than or equal to `hg_mode_projection_threshold` with the supplied
default value 0.20.

This diagnostic is not a replacement for the complete vibrational
partition function. Thermochemistry uses the selected positive
vibrational frequencies of the complete mobile-coordinate system.
Hg-dominated frequencies are reported separately to help interpret the
low-frequency adsorbate motion.

## 6. Thermodynamic state-function construction

The calculation distinguishes the electronic potential-energy
contribution from thermal and vibrational contributions. The adsorbed
state is represented by the optimized Hg/Au slab and its selected
mobile-coordinate harmonic modes; the reference state is the optimized
clean Au slab plus one gas-phase Hg atom.

For the adsorption reaction

<span id="eq-006">
``` math
\text{Hg}(g,T,P) + \text{Au}(111) \rightarrow \text{Hg}^{*}/\text{Au}(111),
```
</span>

the working thermodynamic convention is

<span id="eq-007">
``` math
\Delta G_{\text{ads}}(T,P)
  = \Delta E_{\text{ads}} + \Delta F_{\text{vib}}(T)
    - \mu_{\text{Hg,thermal}}(T,P).
```
</span>

The same quantity is evaluated independently through

<span id="eq-009">
``` math
\Delta G_{\text{ads}} = \Delta H_{\text{ads}} - T \Delta S_{\text{ads}}.
```
</span>

The implementation checks that these two expressions agree numerically.

The electronic term is

<span id="eq-010">
``` math
\Delta E_{\text{ads}} = E_{\text{HgAu}} - E_{\text{Au}} - E_{\text{Hg}}.
```
</span>

The solid-state vibrational free-energy difference is

<span id="eq-011">
``` math
\Delta F_{\text{vib}} = F_{\text{vib,HgAu}} - F_{\text{vib,Au}}.
```
</span>

The gas chemical potential is defined relative to the isolated-Hg
electronic energy already included in $`\Delta E_{\text{ads}}`$; hence
the gas term used here is the thermal ideal-gas contribution rather than
a second electronic atomic-energy term.

### 6.1 Harmonic oscillator energy spectrum

For a normal mode of frequency $`\nu`$,

<span id="eq-012">
``` math
\epsilon = h \nu = h c \tilde{\nu},
```
</span>

where $`\tilde{\nu}`$ is the wavenumber in cm<sup>-1</sup>.

The dimensionless temperature parameter is

<span id="eq-013">
``` math
x = \frac{\epsilon}{k_B T} = \frac{h \nu}{k_B T}.
```
</span>

The quantum harmonic-oscillator energy levels are

``` math
E_n = h \nu \left(n + \tfrac{1}{2}\right), \qquad n = 0, 1, 2, \ldots
```

The single-mode vibrational partition function is

<span id="eq-014">
``` math
q_{\text{vib}}
  = \sum_n \exp\left[-\beta h \nu \left(n + \tfrac{1}{2}\right)\right]
  = \frac{\exp(-x/2)}{1 - \exp(-x)},
```
</span>

with $`\beta = 1/(k_B T)`$.

The corresponding Helmholtz free energy is

<span id="eq-015">
``` math
F_{\text{vib,mode}}
  = -k_B T \ln(q_{\text{vib}})
  = \frac{h \nu}{2} + k_B T \ln\left(1 - e^{-x}\right).
```
</span>

The mean energy is

<span id="eq-016">
``` math
U_{\text{vib,mode}}
  = -\frac{d \ln(q_{\text{vib}})}{d\beta}
  = \frac{h \nu}{2} + \frac{h \nu}{e^{x} - 1}.
```
</span>

Thus the zero-point term and the finite-temperature thermal term are
naturally separated as

<span id="eq-017">
``` math
\text{ZPE}_{\text{mode}} = \frac{h \nu}{2}
```
</span>

and

<span id="eq-018">
``` math
U_{\text{vib,thermal,mode}} = \frac{h \nu}{e^{x} - 1}.
```
</span>

### 6.2 Vibrational entropy from the partition function

For a mode,

<span id="eq-019">
``` math
S_{\text{vib,mode}}
  = \frac{U_{\text{vib,mode}} - F_{\text{vib,mode}}}{T}.
```
</span>

Substitution gives

<span id="eq-020">
``` math
S_{\text{vib,mode}}
  = k_B \left[
      \frac{x}{e^{x} - 1} - \ln\left(1 - e^{-x}\right)
    \right].
```
</span>

For independent harmonic modes, the total quantities are sums over all
retained positive modes.

### 6.3 Heat capacity

Differentiation of the thermal vibrational energy gives

<span id="eq-021">
``` math
C_{V,\text{vib,mode}}
  = \frac{d U_{\text{vib,thermal}}}{dT}
  = k_B \frac{x^2 e^{x}}{(e^{x} - 1)^2}.
```
</span>

The total harmonic vibrational heat capacity is the sum over modes.

The implementation does not require $`C_V`$ explicitly for the
adsorption free energy, but this relation is useful for interpreting the
temperature dependence of the vibrational internal energy and enthalpy.

### 6.4 Thermodynamic identities

For each harmonic state,

<span id="eq-022">
``` math
F_{\text{vib}} = U_{\text{vib,total}} - T S_{\text{vib}},
```
</span>

and therefore

<span id="eq-023">
``` math
\Delta F_{\text{vib}}
  = \Delta U_{\text{vib,total}} - T \Delta S_{\text{vib}}.
```
</span>

When ZPE is included,

<span id="eq-024">
``` math
\Delta U_{\text{vib,total}}
  = \Delta \text{ZPE} + \Delta U_{\text{vib,thermal}},
```
</span>

so that

<span id="eq-025">
``` math
\Delta F_{\text{vib}}
  = \Delta \text{ZPE}
    + \Delta U_{\text{vib,thermal}}
    - T \Delta S_{\text{vib}}.
```
</span>

This identity is the algebraic bridge between the two
adsorption-free-energy routes implemented in the code.

## 7. Vibrational statistical mechanics

For each retained positive frequency $`\nu`$ in cm<sup>-1</sup>, the
quantum of vibrational energy is

<span id="eq-026">
``` math
\epsilon = h c \nu.
```
</span>

In electron-volt units the implementation uses

<span id="eq-027">
``` math
\epsilon (\text{eV}) = 1.2398419843320026 \times 10^{-4} \,
\nu (\text{cm}^{-1}).
```
</span>

### 7.1 Zero-point energy

For a set of harmonic modes,

``` math
\text{ZPE} = \tfrac{1}{2} \sum_k \epsilon_k.
```

The adsorption ZPE correction is

<span id="eq-028">
``` math
\Delta \text{ZPE} = \text{ZPE}_{\text{HgAu}} - \text{ZPE}_{\text{Au}}.
```
</span>

The code retains the raw ZPE difference and applies the `include_zpe`
switch explicitly.

### 7.2 Thermal vibrational internal energy

The thermal excitation energy of one harmonic mode is

<span id="eq-029">
``` math
U_{k,\text{thermal}}
  = \frac{\epsilon_k}{\exp(\epsilon_k / (k_B T)) - 1}.
```
</span>

Therefore

<span id="eq-030">
``` math
U_{\text{vib,thermal}}(T) = \sum_k U_{k,\text{thermal}}.
```
</span>

The adsorption vibrational thermal-energy correction is

<span id="eq-031">
``` math
\Delta U_{\text{vib,thermal}}
  = U_{\text{HgAu,thermal}} - U_{\text{Au,thermal}}.
```
</span>

If ZPE is included, the total vibrational internal energy is

<span id="eq-032">
``` math
U_{\text{vib,total}} = \text{ZPE} + U_{\text{vib,thermal}}.
```
</span>

### 7.3 Vibrational Helmholtz free energy

For a harmonic oscillator, including zero-point energy,

``` math
F_k = \tfrac{1}{2} \epsilon_k
      + k_B T \ln\left[1 - \exp(-\epsilon_k / (k_B T))\right].
```

Thus

<span id="eq-033">
``` math
F_{\text{vib}}(T) = \sum_k F_k
```
</span>

and

<span id="eq-034">
``` math
\Delta F_{\text{vib}}
  = F_{\text{HgAu,vib}} - F_{\text{Au,vib}}.
```
</span>

When ZPE is disabled, the zero-point term is removed consistently from
the reported vibrational free-energy correction.

### 7.4 Vibrational entropy

The vibrational entropy can be obtained from

<span id="eq-035">
``` math
S_{\text{vib}}
  = \frac{U_{\text{vib,total}} - F_{\text{vib}}}{T}.
```
</span>

The adsorption vibrational entropy change is

<span id="eq-036">
``` math
\Delta S_{\text{vib}}
  = S_{\text{HgAu,vib}} - S_{\text{Au,vib}}.
```
</span>

The implementation evaluates this quantity from the same harmonic mode
sets used for the vibrational energy and free energy.

### 7.5 Complete harmonic-thermodynamic formula set

This section collects the complete set of thermodynamic equations used
by the workflow in a form suitable for reproducing the implementation
independently of the Python code. The surface and adsorbate are treated
as harmonic vibrational systems, while gas-phase Hg is treated as a
monatomic ideal gas. All thermodynamic quantities below are per adsorbed
Hg atom unless explicitly stated otherwise.

#### 7.5.1 Energy units and frequency conversion

The code uses eV for energies and eV/K for entropies. Vibrational
frequencies are stored as wavenumbers in cm<sup>-1</sup>. The
photon-like quantum associated with a vibrational wavenumber
$`\bar{\nu}`$ is

<span id="eq-037">
``` math
\epsilon = h c \bar{\nu}.
```
</span>

In the project's eV units,

<span id="eq-038">
``` math
\epsilon\,[\text{eV}] = (hc)_{\text{eV}\cdot\text{cm}} \,
\bar{\nu}\,[\text{cm}^{-1}],
```
</span>

where `HC_EV_CM` is the corresponding conversion constant. The
dimensionless harmonic-oscillator variable is

<span id="eq-039">
``` math
x = \frac{\epsilon}{k_B T}.
```
</span>

The same $`x`$ is used in the energy, entropy, and free-energy
expressions.

#### 7.5.2 One harmonic vibrational mode

For one mode of quantum energy $`\epsilon = h \nu`$, the harmonic-
oscillator energy levels are

``` math
E_n = \left(n + \tfrac{1}{2}\right) \epsilon,
\qquad n = 0, 1, 2, \ldots
```

The canonical partition function is

``` math
q_{\text{vib}}
  = \sum_{n=0}^{\infty} \exp\left[-\beta \left(n + \tfrac{1}{2}\right)
    \epsilon\right]
  = \frac{\exp(-\beta \epsilon / 2)}{1 - \exp(-\beta \epsilon)},
```

with $`\beta = 1/(k_B T)`$.

For a set of independent modes, the total vibrational partition function
is the product

<span id="eq-040">
``` math
Q_{\text{vib}} = \prod_k q_k.
```
</span>

Equivalently,

<span id="eq-041">
``` math
\ln Q_{\text{vib}} = \sum_k \ln q_k.
```
</span>

#### 7.5.3 Vibrational zero-point energy

The zero-point energy of one harmonic mode is

``` math
\text{ZPE}_k = \tfrac{1}{2} \epsilon_k.
```

For all retained positive modes,

``` math
\text{ZPE}
  = \tfrac{1}{2} \sum_k \epsilon_k
  = \tfrac{1}{2} \sum_k h c \bar{\nu}_k.
```

For adsorption, the vibrational zero-point correction is the difference
between the adsorbed system and the clean surface:

<span id="eq-042">
``` math
\Delta \text{ZPE} = \text{ZPE}_{\text{HgAu}} - \text{ZPE}_{\text{Au}}.
```
</span>

If `include_zpe = false`, the implementation sets this adsorption
correction to zero in the thermochemical combination. The raw ZPE is
still available as a diagnostic quantity, and the ASE free-energy output
is adjusted so that the switch has the same meaning in both backends.

#### 7.5.4 Thermal vibrational internal energy

For one harmonic mode, excluding its zero-point contribution,

<span id="eq-043">
``` math
U_{\text{vib,thermal},k}
  = \frac{\epsilon_k}{\exp(\epsilon_k / (k_B T)) - 1}
  = \frac{\epsilon_k}{e^{x_k} - 1}.
```
</span>

Thus

<span id="eq-044">
``` math
U_{\text{vib,thermal}}(T)
  = \sum_k \frac{\epsilon_k}{e^{x_k} - 1}.
```
</span>

The total harmonic vibrational internal energy, when ZPE is included, is

<span id="eq-045">
``` math
U_{\text{vib,total}}(T)
  = \text{ZPE} + U_{\text{vib,thermal}}(T).
```
</span>

For an adsorption process,

<span id="eq-046">
``` math
\Delta U_{\text{vib,thermal}}
  = U_{\text{vib,thermal,HgAu}} - U_{\text{vib,thermal,Au}},
```
</span>

and

<span id="eq-047">
``` math
\Delta U_{\text{vib,total}}
  = \Delta \text{ZPE} + \Delta U_{\text{vib,thermal}}
```
</span>

when ZPE is enabled.

The implementation evaluates the Bose-Einstein denominator numerically
as `expm1(x)`. For very large $`x`$, the thermal contribution is
numerically set to zero once it is below floating-point relevance.

#### 7.5.5 Vibrational Helmholtz free energy

For one harmonic mode,

<span id="eq-048">
``` math
F_{\text{vib},k} = -k_B T \ln(q_k).
```
</span>

Using the harmonic partition function gives

<span id="eq-049">
``` math
F_{\text{vib},k}
  = \tfrac{1}{2} \epsilon_k
    + k_B T \ln\left[1 - \exp(-\epsilon_k / (k_B T))\right].
```
</span>

Therefore,

<span id="eq-050">
``` math
F_{\text{vib}}(T)
  = \sum_k \left[
      \tfrac{1}{2} \epsilon_k + k_B T \ln(1 - e^{-x_k})
    \right].
```
</span>

The free energy can also be written as

<span id="eq-051">
``` math
F_{\text{vib}}
  = \text{ZPE} + k_B T \sum_k \ln(1 - e^{-x_k}).
```
</span>

The adsorption vibrational free-energy correction is

<span id="eq-052">
``` math
\Delta F_{\text{vib}}(T)
  = F_{\text{vib,HgAu}}(T) - F_{\text{vib,Au}}(T).
```
</span>

If ZPE is excluded by the project option, the corresponding ZPE term is
removed consistently from this expression. With ZPE included, the
identity

<span id="eq-053">
``` math
F_{\text{vib}} = U_{\text{vib,total}} - T S_{\text{vib}}
```
</span>

holds for the harmonic oscillator.

#### 7.5.6 Vibrational entropy

For one harmonic mode, the entropy is

<span id="eq-054">
``` math
S_{\text{vib},k}
  = k_B \left[
      \frac{x_k}{e^{x_k} - 1} - \ln(1 - e^{-x_k})
    \right].
```
</span>

Hence

<span id="eq-055">
``` math
S_{\text{vib}}(T)
  = k_B \sum_k \left[
      \frac{x_k}{e^{x_k} - 1} - \ln(1 - e^{-x_k})
    \right].
```
</span>

The equivalent thermodynamic identity used as a check is

<span id="eq-056">
``` math
S_{\text{vib}}
  = \frac{U_{\text{vib,total}} - F_{\text{vib}}}{T}.
```
</span>

Because the zero-point term cancels between $`U`$ and $`F`$, the same
entropy results whether $`U`$ is written as ZPE plus thermal energy or
in the corresponding partition-function form.

The adsorption vibrational entropy change is

<span id="eq-057">
``` math
\Delta S_{\text{vib}}(T)
  = S_{\text{vib,HgAu}}(T) - S_{\text{vib,Au}}(T).
```
</span>

At fixed positive frequency, the low-temperature limit is
$`S_{\text{vib}} \to 0`$. For a very soft mode, however, the harmonic
entropy becomes large and therefore the numerical treatment of low
frequencies is particularly important.

#### 7.5.7 Vibrational heat capacity (derived quantity)

Although heat capacity is not required as an independent input to the
final $`\Delta G`$ calculation, it follows from the same harmonic model
and is useful for interpreting the temperature dependence:

<span id="eq-058">
``` math
C_{V,\text{vib}}
  = \left(\frac{dU_{\text{vib,total}}}{dT}\right)_V
  = k_B \sum_k \frac{x_k^2 e^{x_k}}{(e^{x_k} - 1)^2}.
```
</span>

The zero-point term has zero temperature derivative and therefore does
not contribute to $`C_V`$.

#### 7.5.8 Thermodynamic identities used by the implementation

For the harmonic surface model,

<span id="eq-059">
``` math
F_{\text{vib}} = U_{\text{vib,total}} - T S_{\text{vib}},
```
</span>

and therefore

<span id="eq-060">
``` math
\Delta F_{\text{vib}}
  = \Delta U_{\text{vib,total}} - T \Delta S_{\text{vib}}.
```
</span>

When ZPE is included,

<span id="eq-061">
``` math
\Delta U_{\text{vib,total}}
  = \Delta \text{ZPE} + \Delta U_{\text{vib,thermal}}.
```
</span>

Consequently,

<span id="eq-062">
``` math
\Delta F_{\text{vib}}
  = \Delta \text{ZPE} + \Delta U_{\text{vib,thermal}}
    - T \Delta S_{\text{vib}}.
```
</span>

This identity is the algebraic bridge between the explicit
$`\Delta H - T \Delta S`$ route and the compact
$`\Delta E + \Delta F - \mu`$ route used in the code.

#### 7.5.9 Clean-surface and adsorbed-system definitions

The two solid-state vibrational states are treated separately:

``` math
\text{Au}: \{\nu_k^{\text{Au}}\}, \qquad
\text{HgAu}: \{\nu_k^{\text{HgAu}}\}.
```

For any vibrational function $`X`$,

<span id="eq-063">
``` math
\Delta X_{\text{vib}}
  = X_{\text{vib,HgAu}} - X_{\text{vib,Au}}.
```
</span>

Thus, explicitly,

<span id="eq-064">
``` math
\Delta \text{ZPE}
  = \tfrac{1}{2} \sum_k \epsilon_k^{\text{HgAu}}
    - \tfrac{1}{2} \sum_j \epsilon_j^{\text{Au}},
```
</span>

<span id="eq-065">
``` math
\Delta U_{\text{vib,thermal}}(T)
  = \sum_k \frac{\epsilon_k^{\text{HgAu}}}
                 {\exp(\epsilon_k^{\text{HgAu}} / (k_B T)) - 1}
    - \sum_j \frac{\epsilon_j^{\text{Au}}}
                 {\exp(\epsilon_j^{\text{Au}} / (k_B T)) - 1},
```
</span>

<span id="eq-066">
``` math
\Delta F_{\text{vib}}(T)
  = \sum_k \left[
      \frac{\epsilon_k^{\text{HgAu}}}{2}
      + k_B T \ln\left(1 - e^{-\epsilon_k^{\text{HgAu}} / (k_B T)}\right)
    \right]
    - \sum_j \left[
      \frac{\epsilon_j^{\text{Au}}}{2}
      + k_B T \ln\left(1 - e^{-\epsilon_j^{\text{Au}} / (k_B T)}\right)
    \right],
```
</span>

<span id="eq-067">
``` math
\Delta S_{\text{vib}}(T)
  = k_B \sum_k \left[
      \frac{x_k^{\text{HgAu}}}{e^{x_k^{\text{HgAu}}} - 1}
      - \ln\left(1 - e^{-x_k^{\text{HgAu}}}\right)
    \right]
    - k_B \sum_j \left[
      \frac{x_j^{\text{Au}}}{e^{x_j^{\text{Au}}} - 1}
      - \ln\left(1 - e^{-x_j^{\text{Au}}}\right)
    \right].
```
</span>

Only frequencies passing the positive-frequency selection are included
in these thermodynamic sums.

## 8. Complete ideal-gas Hg statistical mechanics

The gas reference is one monatomic Hg atom. There are no molecular
rotational or vibrational degrees of freedom. The thermal gas model
therefore contains translational motion only.

### 8.1 Translational partition function

For $`N`$ identical monatomic particles in volume $`V`$, the
translational partition function is

<span id="eq-068">
``` math
Q_{\text{trans}} = \frac{q_{\text{trans}}^{N}}{N!},
```
</span>

where the one-particle partition function is

<span id="eq-069">
``` math
q_{\text{trans}}
  = \frac{V}{\Lambda^{3}}
  = V \left(\frac{2 \pi m k_B T}{h^{2}}\right)^{3/2},
```
</span>

and the thermal de Broglie wavelength is

<span id="eq-070">
``` math
\Lambda = \frac{h}{\sqrt{2 \pi m k_B T}}.
```
</span>

The Hg atomic mass used by the implementation is 200.59 u, converted to
kg per atom where required.

### 8.2 Ideal-gas entropy

Using Stirling's approximation for $`N!`$ gives the Sackur--Tetrode form

<span id="eq-071">
``` math
S_{\text{trans}}
  = N k_B \left[
      \ln\left(\frac{V}{N \Lambda^{3}}\right) + \tfrac{5}{2}
    \right].
```
</span>

Per Hg atom,

<span id="eq-072">
``` math
S_{\text{Hg}}
  = k_B \left[
      \ln\left(\frac{V}{N \Lambda^{3}}\right) + \tfrac{5}{2}
    \right].
```
</span>

Using the ideal-gas equation

<span id="eq-073">
``` math
P V = N k_B T
```
</span>

gives

<span id="eq-074">
``` math
S_{\text{Hg}}(T,P)
  = k_B \left[
      \ln\left(\frac{k_B T}{P \Lambda^{3}}\right) + \tfrac{5}{2}
    \right].
```
</span>

Relative to a reference pressure $`P_0`$,

<span id="eq-075">
``` math
S_{\text{Hg}}(T,P)
  = S_{\text{Hg}}(T,P_0) - k_B \ln(P / P_0).
```
</span>

Thus lowering the gas pressure increases the gas entropy and makes
adsorption less favorable through the chemical-potential term.

### 8.3 Ideal-gas internal energy and enthalpy

A monatomic ideal gas has three translational quadratic degrees of
freedom. Therefore, per atom,

``` math
U_{\text{gas}} = \tfrac{3}{2} k_B T.
```

The ideal-gas equation gives

<span id="eq-076">
``` math
P V = k_B T
```
</span>

per atom, and hence

<span id="eq-077">
``` math
H_{\text{gas}} = U_{\text{gas}} + P V
               = \tfrac{5}{2} k_B T.
```
</span>

The implementation calls this quantity `H_gas_thermal`.

### 8.4 Gas chemical potential

The thermal chemical-potential contribution is evaluated as

<span id="eq-078">
``` math
\mu_{\text{Hg,thermal}}(T,P)
  = H_{\text{gas,thermal}}(T) - T S_{\text{Hg}}(T,P).
```
</span>

Equivalently,

<span id="eq-079">
``` math
\mu_{\text{Hg,thermal}}
  = -k_B T \ln(q_{\text{trans}} / N)
```
</span>

for the corresponding classical ideal-gas reference, with the same
pressure, temperature, and standard-state convention.

The pressure dependence can be written as

<span id="eq-080">
``` math
\mu_{\text{Hg,thermal}}(T,P)
  = \mu_{\text{Hg,thermal}}(T,P_0)
    + k_B T \ln(P / P_0).
```
</span>

Consequently, increasing Hg pressure makes the gas chemical potential
less negative and adsorption thermodynamically more favorable.

### 8.5 Unit conversion

The project uses

<span id="eq-081">
``` math
1\,\text{bar} = 10^{5}\,\text{Pa}
```
</span>

and

<span id="eq-082">
``` math
1\,\text{eV} = 1.602176634 \times 10^{-19}\,\text{J}.
```
</span>

For vibrational wavenumbers,

<span id="eq-083">
``` math
h c = 1.2398419843320026 \times 10^{-4}\,\text{eV}\,\text{cm},
```
</span>

so

<span id="eq-084">
``` math
\epsilon\,(\text{eV})
  = 1.2398419843320026 \times 10^{-4} \,
    \tilde{\nu}\,(\text{cm}^{-1}).
```
</span>

## 9. Gas-phase Hg thermodynamics

The reference state for Hg is a monatomic ideal gas. The gas reference
is not represented by a harmonic molecular vibration calculation because
an isolated atom has no molecular vibrational or rotational degrees of
freedom.

### 9.1 Ideal-gas PV contribution

For one mole of ideal gas,

<span id="eq-085">
``` math
P V = R T.
```
</span>

Per particle,

<span id="eq-086">
``` math
P V = k_B T.
```
</span>

Therefore the thermal enthalpy contribution of a monatomic ideal gas is

``` math
H_{\text{gas,thermal}} = \tfrac{5}{2} k_B T.
```

This term is subtracted from the adsorbed-state enthalpy correction
because the gas reference contains translational PV work whereas the
adsorbed atom does not.

### 9.2 Translational entropy

The translational partition function is

<span id="eq-087">
``` math
q_{\text{trans}}
  = \left(\frac{2 \pi m k_B T}{h^{2}}\right)^{3/2} V.
```
</span>

The Sackur--Tetrode expression used by the custom backend can be written
as

<span id="eq-088">
``` math
S_{\text{trans}}
  = k_B \left[\ln(q_{\text{trans}} / N) + \tfrac{5}{2}\right].
```
</span>

Using the ideal-gas relation

<span id="eq-089">
``` math
V = \frac{N k_B T}{P}
```
</span>

introduces the pressure dependence. Relative to a reference pressure
$`P_0`$,

<span id="eq-090">
``` math
S(T,P) = S(T,P_0) - k_B \ln(P / P_0).
```
</span>

The gas chemical potential contribution is represented thermodynamically
as

<span id="eq-091">
``` math
\mu_{\text{Hg,thermal}}(T,P)
  = H_{\text{gas,thermal}} - T S_{\text{gas}}(T,P).
```
</span>

The pressure used by the project is entered in bar and converted
internally to SI pressure where required.

### 9.3 Complete ideal-gas Hg statistical mechanics

The gas reference is a single monatomic Hg atom. Internal molecular
rotational and vibrational contributions are absent because Hg is
monatomic. The model therefore contains translational motion only.

#### 9.3.1 Translational partition function

For one classical particle in a volume $`V`$, the translational
partition function is

<span id="eq-092">
``` math
q_{\text{trans}} = \frac{V}{\Lambda^{3}},
```
</span>

where the thermal de Broglie wavelength is

<span id="eq-093">
``` math
\Lambda = \frac{h}{\sqrt{2 \pi m k_B T}}.
```
</span>

For $`N`$ indistinguishable particles,

<span id="eq-094">
``` math
Q_{\text{trans}} = \frac{q_{\text{trans}}^{N}}{N!}.
```
</span>

Using Stirling's approximation for large $`N`$ and the thermodynamic
limit gives the standard ideal-gas expressions used by the
Sackur--Tetrode formulation.

For a monatomic ideal gas,

``` math
U_{\text{trans}} = \tfrac{3}{2} N k_B T,
```

and

``` math
H_{\text{trans}} = U_{\text{trans}} + P V
                 = \tfrac{5}{2} N k_B T,
```

because $`P V = N k_B T`$.

Per Hg atom, the thermal enthalpy used by the code is therefore

``` math
H_{\text{gas,thermal}}(T) = \tfrac{5}{2} k_B T.
```

The isolated-Hg electronic energy is already included separately through
$`\Delta E_{\text{ads}}`$; consequently this gas correction contains
only the thermal ideal-gas contribution.

#### 9.3.2 Sackur--Tetrode entropy

For a monatomic ideal gas, the entropy per particle can be written

<span id="eq-095">
``` math
\frac{S}{N}
  = k_B \left[
      \ln\left(\frac{V}{N \Lambda^{3}}\right) + \tfrac{5}{2}
    \right].
```
</span>

Using the ideal-gas equation

<span id="eq-096">
``` math
P V = N k_B T
```
</span>

gives

<span id="eq-097">
``` math
\frac{V}{N} = \frac{k_B T}{P},
```
</span>

and hence

<span id="eq-098">
``` math
S_{\text{gas}}(T,P)
  = k_B \left[
      \ln\left(\frac{k_B T}{P \Lambda^{3}}\right) + \tfrac{5}{2}
    \right].
```
</span>

The pressure dependence at fixed $`T`$ is therefore

<span id="eq-099">
``` math
S_{\text{gas}}(T,P)
  = S_{\text{gas}}(T,P_0) - k_B \ln(P / P_0),
```
</span>

where $`P_0`$ is the configured reference pressure.

This is the explicit pressure correction implemented by the custom
backend. The pressure is supplied in bar in the input file and converted
to pascal for this SI statistical-mechanical expression.

#### 9.3.3 Chemical potential

For an ideal gas, the chemical potential can be written as

<span id="eq-100">
``` math
\mu = -k_B T \ln(q_{\text{trans}} / N)
```
</span>

up to the standard indistinguishability formulation, or equivalently
through thermodynamic identities. The implementation uses the
enthalpy-entropy form for the thermal contribution:

<span id="eq-101">
``` math
\mu_{\text{Hg,thermal}}(T,P)
  = H_{\text{gas,thermal}}(T) - T S_{\text{gas}}(T,P).
```
</span>

Substitution of the monatomic ideal-gas enthalpy gives

<span id="eq-102">
``` math
\mu_{\text{Hg,thermal}}(T,P)
  = \tfrac{5}{2} k_B T - T S_{\text{gas}}(T,P).
```
</span>

Because

<span id="eq-103">
``` math
S_{\text{gas}}(T,P)
  = S_{\text{gas}}(T,P_0) - k_B \ln(P / P_0),
```
</span>

one obtains the explicit pressure dependence

<span id="eq-104">
``` math
\mu_{\text{Hg,thermal}}(T,P)
  = \mu_{\text{Hg,thermal}}(T,P_0)
    + k_B T \ln(P / P_0).
```
</span>

Thus increasing the gas pressure increases the gas chemical potential.
Since the adsorption free energy contains $`-\mu_{\text{Hg,thermal}}`$,
increasing Hg pressure makes adsorption thermodynamically more favorable
in this convention.

#### 9.3.4 Standard pressure versus actual pressure

The code distinguishes the actual Hg pressure $`P`$ from a reference
pressure $`P_0`$. The standard-pressure gas entropy is

<span id="eq-105">
``` math
S_{\text{gas,standard}}(T) = S_{\text{gas}}(T,P_0),
```
</span>

and the actual-pressure entropy is

<span id="eq-106">
``` math
S_{\text{gas}}(T,P)
  = S_{\text{gas,standard}}(T) - k_B \ln(P / P_0).
```
</span>

The thermodynamic result therefore depends on both the configured
pressure and the configured reference pressure. The latter fixes the
zero of the gas chemical-potential convention; the former represents the
physical gas condition being evaluated.

#### 9.3.5 Pressure derivative of the adsorption free energy

At fixed temperature, the ideal-gas contribution implies

<span id="eq-107">
``` math
\frac{d \mu_{\text{Hg,thermal}}}{d \ln P} = k_B T.
```
</span>

Therefore, for the adsorption free energy

<span id="eq-108">
``` math
\Delta G_{\text{ads}}
  = \Delta E_{\text{ads}} + \Delta F_{\text{vib}}
    - \mu_{\text{Hg,thermal}},
```
</span>

one has

<span id="eq-109">
``` math
\frac{d \Delta G_{\text{ads}}}{d \ln P} = -k_B T.
```
</span>

This relation is a useful analytical check on pressure-dependent scans.

## 10. Adsorption enthalpy, entropy and Gibbs free energy

### 10.1 Reaction definition and thermodynamic state functions

The modeled adsorption reaction is

<span id="eq-110">
``` math
\text{Hg}(g) + \text{Au(surface)}
  \rightarrow \text{Hg}^{*}/\text{Au(surface)},
```
</span>

where the asterisk denotes the adsorbed Hg atom in the selected
optimized adsorption state. The electronic-energy contribution is

<span id="eq-111">
``` math
\Delta E_{\text{ads}}
  = E_{\text{HgAu}} - E_{\text{Au}} - E_{\text{Hg}}.
```
</span>

The thermodynamic adsorption free energy is constructed by adding the
vibrational correction for the solid states and subtracting the
gas-phase Hg chemical-potential contribution.

The central state-function definitions are

<span id="eq-112">
``` math
H = U + P V,
```
</span>

<span id="eq-113">
``` math
F = U - T S,
```
</span>

<span id="eq-114">
``` math
G = H - T S = F + P V,
```
</span>

and for a gas-phase species the relevant reservoir quantity is its
chemical potential $`\mu`$ (Gibbs free energy per particle in the
ideal-gas limit).

For the adsorption reaction, the code uses

<span id="eq-115">
``` math
\Delta H_{\text{ads}}
  = H_{\text{HgAu}} - H_{\text{Au}} - H_{\text{Hg}}(g),
```
</span>

<span id="eq-116">
``` math
\Delta S_{\text{ads}}
  = S_{\text{HgAu}} - S_{\text{Au}} - S_{\text{Hg}}(g),
```
</span>

<span id="eq-117">
``` math
\Delta G_{\text{ads}}
  = \Delta H_{\text{ads}} - T \Delta S_{\text{ads}}.
```
</span>

Because the electronic energy is separated from the thermal gas
contribution, these become the explicit equations given below.

### 10.2 Explicit enthalpy equation

The surface and adsorbate electronic energies are supplied by MACE. The
vibrational correction to their internal energies is harmonic. The gas
phase contributes its thermal ideal-gas enthalpy. Thus

<span id="eq-118">
``` math
\Delta H_{\text{ads}}(T)
  = \Delta E_{\text{ads}}
    + \Delta \text{ZPE}
    + \Delta U_{\text{vib,thermal}}(T)
    - H_{\text{gas,thermal}}(T),
```
</span>

with

``` math
H_{\text{gas,thermal}}(T) = \tfrac{5}{2} k_B T.
```

If ZPE is disabled,

<span id="eq-119">
``` math
\Delta H_{\text{ads}}(T)
  = \Delta E_{\text{ads}}
    + \Delta U_{\text{vib,thermal}}(T)
    - H_{\text{gas,thermal}}(T).
```
</span>

### 10.3 Explicit entropy equation

The solid-state entropy change is

<span id="eq-120">
``` math
\Delta S_{\text{vib}}(T)
  = S_{\text{vib,HgAu}}(T) - S_{\text{vib,Au}}(T).
```
</span>

The reaction entropy is

<span id="eq-121">
``` math
\Delta S_{\text{ads}}(T,P)
  = \Delta S_{\text{vib}}(T) - S_{\text{gas}}(T,P).
```
</span>

With the ideal-gas expression,

<span id="eq-122">
``` math
\Delta S_{\text{ads}}(T,P)
  = \Delta S_{\text{vib}}(T)
    - k_B \left[
        \ln\left(\frac{k_B T}{P \Lambda^{3}}\right) + \tfrac{5}{2}
      \right].
```
</span>

Equivalently, relative to $`P_0`$,

<span id="eq-123">
``` math
\Delta S_{\text{ads}}(T,P)
  = \Delta S_{\text{ads}}(T,P_0) + k_B \ln(P / P_0).
```
</span>

### 10.4 Explicit Gibbs-energy equation

Combining the preceding expressions gives

<span id="eq-124">
``` math
\Delta G_{\text{ads}}(T,P)
  = \Delta E_{\text{ads}}
    + \Delta \text{ZPE}
    + \Delta U_{\text{vib,thermal}}(T)
    - H_{\text{gas,thermal}}(T)
    - T \Delta S_{\text{vib}}(T)
    + T S_{\text{gas}}(T,P).
```
</span>

Using

<span id="eq-125">
``` math
\Delta F_{\text{vib}}
  = \Delta \text{ZPE} + \Delta U_{\text{vib,thermal}}
    - T \Delta S_{\text{vib}}
```
</span>

and

<span id="eq-126">
``` math
\mu_{\text{Hg,thermal}}
  = H_{\text{gas,thermal}} - T S_{\text{gas}},
```
</span>

this reduces exactly to

<span id="eq-127">
``` math
\Delta G_{\text{ads}}(T,P)
  = \Delta E_{\text{ads}} + \Delta F_{\text{vib}}(T)
    - \mu_{\text{Hg,thermal}}(T,P).
```
</span>

This is the compact expression used for the internal cross-check.

### 10.5 Separation into electronic, vibrational and gas terms

It is useful to display the free energy as

<span id="eq-128">
``` math
\Delta G_{\text{ads}}
  = \Delta E_{\text{ads}}
    + \left[\Delta \text{ZPE}
      + \Delta U_{\text{vib,thermal}}
      - T \Delta S_{\text{vib}}\right]
    - \left[H_{\text{gas,thermal}} - T S_{\text{gas}}\right].
```
</span>

The three physically distinct contributions are therefore:

- electronic adsorption energy: $`\Delta E_{\text{ads}}`$;
- surface/adsorbate vibrational free-energy correction:
  $`\Delta F_{\text{vib}}`$; and
- gas-phase Hg thermal chemical potential: $`\mu_{\text{Hg,thermal}}`$.

No gas-phase Hg electronic energy is added a second time:
$`E_{\text{Hg}}`$ is already part of $`\Delta E_{\text{ads}}`$.

### 10.6 Pressure dependence of Delta G_ads

At fixed $`T`$, the pressure dependence follows directly from the
ideal-gas chemical potential:

<span id="eq-129">
``` math
\mu_{\text{Hg}}(T,P)
  = \mu_{\text{Hg}}(T,P_0) + k_B T \ln(P / P_0).
```
</span>

Therefore

<span id="eq-130">
``` math
\Delta G_{\text{ads}}(T,P)
  = \Delta G_{\text{ads}}(T,P_0) - k_B T \ln(P / P_0).
```
</span>

This is the analytical pressure dependence of the present model.

### 10.7 Alternative free-energy bookkeeping

The code's preferred compact route is

<span id="eq-131">
``` math
\Delta G_{\text{ads}}
  = \Delta E_{\text{ads}} + \Delta F_{\text{vib}}
    - \mu_{\text{Hg,thermal}}.
```
</span>

The expanded route is

<span id="eq-132">
``` math
\Delta G_{\text{ads}}
  = \Delta H_{\text{ads}} - T \Delta S_{\text{ads}}.
```
</span>

The implementation evaluates both and raises an error if they disagree
beyond the configured numerical tolerance. Agreement is expected because

<span id="eq-133">
``` math
\Delta F_{\text{vib}}
  = \Delta \text{ZPE} + \Delta U_{\text{vib,thermal}}
    - T \Delta S_{\text{vib}}
```
</span>

and

<span id="eq-134">
``` math
\mu_{\text{Hg,thermal}}
  = H_{\text{gas,thermal}} - T S_{\text{gas}}.
```
</span>

### 10.8 ZPE switch semantics

When `include_zpe = true`:

<span id="eq-135">
``` math
\Delta H_{\text{ads}}
  = \Delta E_{\text{ads}} + \Delta \text{ZPE}
    + \Delta U_{\text{vib,thermal}} - H_{\text{gas}},
```
</span>

<span id="eq-136">
``` math
\Delta F_{\text{vib}}
  = \Delta \text{ZPE} + \Delta U_{\text{vib,thermal}}
    - T \Delta S_{\text{vib}}.
```
</span>

When `include_zpe = false`:

<span id="eq-137">
``` math
\Delta H_{\text{ads}}
  = \Delta E_{\text{ads}} + \Delta U_{\text{vib,thermal}}
    - H_{\text{gas}},
```
</span>

and the vibrational free-energy correction is correspondingly evaluated
with the ZPE term removed. Entropy is unchanged because ZPE is
temperature independent.

The implemented adsorption enthalpy is

<span id="eq-138">
``` math
\Delta H_{\text{ads}}(T)
  = \Delta E_{\text{ads}}
    + \Delta \text{ZPE}
    + \Delta U_{\text{vib,thermal}}(T)
    - H_{\text{gas,thermal}}(T).
```
</span>

The adsorption entropy is

<span id="eq-139">
``` math
\Delta S_{\text{ads}}(T,P)
  = \Delta S_{\text{vib}}(T) - S_{\text{gas}}(T,P).
```
</span>

The adsorption Gibbs free energy is then

<span id="eq-140">
``` math
\Delta G_{\text{ads}}(T,P)
  = \Delta H_{\text{ads}}(T,P) - T \Delta S_{\text{ads}}(T,P).
```
</span>

An algebraically equivalent and useful implementation check is

<span id="eq-141">
``` math
\Delta G_{\text{ads}}
  = \Delta E_{\text{ads}} + \Delta F_{\text{vib}}
    - \mu_{\text{Hg,thermal}}.
```
</span>

The program evaluates both routes and checks that they agree within
numerical tolerance. This provides an internal thermodynamic-consistency
test.

The sign convention is therefore:

- $`\Delta G_{\text{ads}} < 0`$: adsorption is thermodynamically
  favorable under the stated $`T`$ and $`P`$ conditions;
- $`\Delta G_{\text{ads}} > 0`$: desorption is thermodynamically favored
  relative to the chosen gas reference;
- $`\Delta G_{\text{ads}} = 0`$: adsorption and desorption are at the
  thermodynamic crossover in this model.

## 11. Temperature scan

The temperature scan does not repeat the electronic-structure
calculation and does not recompute the Hessian. Once the optimized
structures, adsorption energy, and vibrational frequencies are
available, the scan repeatedly re-evaluates the analytical thermodynamic
expressions at each configured temperature.

For each $`T`$ the program evaluates and records at least:

- $`\Delta E_{\text{ads}}`$;
- $`\Delta \text{ZPE}`$;
- $`\Delta U_{\text{vib,thermal}}`$;
- $`\Delta U_{\text{vib,total}}`$;
- $`\Delta F_{\text{vib}}`$;
- $`S_{\text{gas}}`$;
- $`\Delta S_{\text{vib}}`$;
- $`\Delta S_{\text{ads}}`$;
- $`H_{\text{gas,thermal}}`$;
- $`\Delta H_{\text{ads}}`$; and
- $`\Delta G_{\text{ads}}`$.

The scan is written to `hg_au111_thermodynamics_scan.csv` when the
output prefix is the supplied default. In `mode = scan` the values are
also printed in a compact terminal table.

For the representative scan in the project documentation, the values
are:

<table style="width:94%;">
<colgroup>
<col style="width: 15%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads</th>
<th>Delta S_ads</th>
<th>Delta G_ads</th>
</tr>
</thead>
<tbody>
<tr>
<td>50.00</td>
<td>-0.49620164 eV</td>
<td>-2.12731440e-3 eV/K</td>
<td>-0.38983592 eV</td>
</tr>
<tr>
<td>100.00</td>
<td>-0.49539024 eV</td>
<td>-2.11729969e-3 eV/K</td>
<td>-0.28366027 eV</td>
</tr>
<tr>
<td>200.00</td>
<td>-0.49187176 eV</td>
<td>-2.09333473e-3 eV/K</td>
<td>-0.07320481 eV</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48791444 eV</td>
<td>-2.07726519e-3 eV/K</td>
<td><blockquote>
<p>0.13142218 eV</p>
</blockquote></td>
</tr>
<tr>
<td>400.00</td>
<td>-0.48366842 eV</td>
<td>-2.06502003e-3 eV/K</td>
<td><blockquote>
<p>0.34233959 eV</p>
</blockquote></td>
</tr>
<tr>
<td>500.00</td>
<td>-0.47944343 eV</td>
<td>-2.05559378e-3 eV/K</td>
<td><blockquote>
<p>0.54835346 eV</p>
</blockquote></td>
</tr>
<tr>
<td>600.00</td>
<td>-0.47519066 eV</td>
<td>-2.04784062e-3 eV/K</td>
<td><blockquote>
<p>0.75351372 eV</p>
</blockquote></td>
</tr>
</tbody>
</table>

## 12. Determination of T_cross

The code can determine the temperature at which
$`\Delta G_{\text{ads}}(T) = 0`$.

The configured `temperatures` list defines the search interval. The code
sorts the positive temperatures, evaluates $`\Delta G`$ on that grid,
and searches for the first adjacent sign change.

If a grid point is numerically zero, it is accepted directly. Otherwise,
if a sign change is present, a deterministic bisection method is
applied.

The bisection uses a temperature tolerance of 1e-4 K and a Gibbs-energy
tolerance of 1e-10 eV, with a maximum of 100 iterations. Each bisection
step only evaluates the thermodynamic equations; no MACE force
calculation and no phonon calculation is performed.

For the representative calculation,
$`T_{\text{cross}} \approx 235.025029`$ K.

At the reported solution the residual was approximately
$`\Delta G_{\text{ads}}(T_{\text{cross}}) = -3.18 \times 10^{-8}`$ eV,
with a final numerical bracket of approximately
$`[235.02498236,\, 235.02507596]`$ K.

The root should be interpreted as a result of the stated harmonic,
Gamma-point, slab, gas-reference, pressure, and electronic-structure
model. It is not an experimentally universal desorption temperature.

### 12.1 Mathematical definition and numerical root criterion

The crossover temperature is the positive solution of

``` math
f(T) = \Delta G_{\text{ads}}(T,P) = 0.
```

The workflow does not solve this equation symbolically. It evaluates the
thermodynamic function numerically using the already computed electronic
energy and frequency sets.

Given a bracket $`[T_{\text{low}},\, T_{\text{high}}]`$ satisfying

<span id="eq-142">
``` math
f(T_{\text{low}}) \, f(T_{\text{high}}) < 0,
```
</span>

the midpoint is

``` math
T_{\text{mid}} = \tfrac{1}{2}\left(T_{\text{low}} + T_{\text{high}}\right).
```

If

``` math
f(T_{\text{low}}) \, f(T_{\text{mid}}) \le 0,
```

the new bracket is $`[T_{\text{low}},\, T_{\text{mid}}]`$; otherwise it
is $`[T_{\text{mid}},\, T_{\text{high}}]`$. The process is repeated
until either the temperature bracket is narrower than the configured
temperature tolerance or the absolute Gibbs-energy residual is below the
configured Gibbs-energy tolerance.

For a bracket containing a single continuous crossing, bisection is
guaranteed to converge. The reported `T_cross` is therefore a numerical
root of the specific model function $`\Delta G_{\text{ads}}(T,P)`$, not
an independently optimized structure at `T_cross`.

### 12.2 Analytical pressure shift of the crossover condition

Because

<span id="eq-143">
``` math
\Delta G_{\text{ads}}(T,P)
  = \Delta G_{\text{ads}}(T,P_0) - k_B T \ln(P / P_0),
```
</span>

the crossover condition can also be written

<span id="eq-144">
``` math
\Delta G_{\text{ads}}(T_{\text{cross}},P_0)
  = k_B T_{\text{cross}} \ln(P / P_0).
```
</span>

This relation does not by itself give a closed-form $`T_{\text{cross}}`$
because the vibrational terms are temperature dependent, but it provides
a useful check on pressure-dependent crossover calculations.

## 13. Selectable computational backends

### 13.1 Vibrational backend

The configuration keyword is

``` ini
[phonons]
method = custom
```

or

``` ini
[phonons]
method = ase
```

`custom` invokes the project's direct finite-difference Hessian
implementation. `ase` invokes ASE's `Vibrations` implementation. Both
operate on the same optimized structure and the same selected mobile
coordinates; they differ in the Hessian/mode implementation.

The terminal output explicitly reports the requested and actual
vibrational backend so that accidental backend substitution is visible.

### 13.2 Thermodynamics backend

The configuration keyword is

``` ini
[thermodynamics]
method = custom
```

or

``` ini
[thermodynamics]
method = ase
```

The `custom` branch evaluates the project's explicit harmonic and
ideal-gas equations described above.

The `ase` branch uses ASE's `HarmonicThermo` for the vibrational systems
and ASE's `IdealGasThermo` for monatomic Hg(g). The project retains the
same adsorption-energy convention and explicitly handles the
`include_zpe` switch so that the two backends are thermodynamically
comparable.

The thermodynamics backend is reported once at the beginning of the
thermochemistry section. It is deliberately not printed inside the
repeated temperature-scan or bisection evaluations.

## 14. Program architecture

The main execution path is

```
run_hg_au111.py
    |
    +-- config.py
    |
    +-- calculator.py
    |
    +-- workflow.py
            |
            +-- surface.py
            +-- optimization.py
            +-- energies.py
            +-- phonons.py
            +-- thermodynamics.py
            +-- constants.py
```

### 14.1 `run_hg_au111.py`

This is the command-line entry point. It loads the configuration,
applies command-line overrides, creates the MACE calculator, executes
the workflow, and prints final results.

Typical execution is

``` bash
python run_hg_au111.py -c hg_au111_config.ini
```

### 14.2 `config.py`

`Config` stores defaults, reads the INI file, converts strings to
numerical values, and validates incompatible settings.

Important controls include:

- slab dimensions and lattice constant;
- geometry-relaxation frozen layers;
- vibrationally frozen layers;
- adsorption site and Hg height;
- optimization thresholds and optimizer;
- phonon backend and finite-difference settings;
- thermodynamics backend;
- temperature and pressure;
- scan temperatures;
- $`\Delta G = 0`$ search activation; and
- MACE model path.

The validation layer ensures, among other things, that the vibrationally
frozen layer count is not smaller than the geometry-relaxation frozen
layer count and does not exceed the number of Au layers.

### 14.3 `surface.py`

Builds the Au(111) slab, identifies Au layers, applies bottom-layer
geometry constraints, determines adsorption-site coordinates, adds Hg,
and reports surface details.

### 14.4 `optimization.py`

Performs constrained structural optimization using ASE optimizers and
the MACE calculator. The relaxation constraints remain separate from the
later vibrational coordinate selection.

### 14.5 `energies.py`

Contains adsorption-energy bookkeeping and isolated-Hg energy evaluation
and prints the electronic-energy summary.

### 14.6 `phonons.py`

Contains the constrained finite-displacement vibrational machinery. It:

- identifies the Au layers;
- constructs the vibrational frozen mask;
- retains all atoms in force evaluations;
- identifies mobile vibrational coordinates;
- dispatches to `custom` or `ase`;
- builds/reads the mode information used by the workflow;
- identifies significant imaginary modes;
- applies the positive-frequency cutoff;
- performs Hg-mode projection diagnostics; and
- writes frequency and diagnostic files.

The important conceptual distinction is that `phonon_atoms` represent
the coordinate space, not the physical system passed to the calculator.

### 14.7 `thermodynamics.py`

Implements the harmonic oscillator formulas, ideal-gas Hg
thermodynamics, the custom thermochemistry backend, and the ASE
thermochemistry backend. It also contains internal consistency checks
between the H-TS and F-mu routes.

### 14.8 `workflow.py`

Coordinates the complete calculation. It builds and optimizes both
systems, calculates the adsorption energy, obtains phonons, evaluates
thermochemistry, performs the scan, optionally searches for `T_cross`,
and writes the final JSON results.

The temperature scan and `T_cross` search intentionally reuse the
already calculated vibrational frequencies.

### 14.9 `calculator.py`

Creates the MACE calculator and provides conservative CUDA-memory
cleanup functions around energy and force evaluations.

### 14.10 `constants.py`

Contains the physical constants and unit-conversion factors needed for
the vibrational frequency and thermodynamic calculations.

## 15. Input configuration

A representative configuration is structured as follows:

``` ini
[system]
vacuum_size = 18.0
surface_size = 4,4,6
lattice_constant = 4.08
n_frozen_layers = 2
n_frozen_layers_vibrations = 6
adsorption_site = fcc
hg_height = 3.0

[calculation]
mode = scan
single_temperature = 298.15
temperatures = 50,100,200,298.15,400,500,600
calculate_g_ads_zero_temperature = true
pressure = 1.0
output_prefix = hg_au111   verbose = true
save_trajectories = true
optimizer = FIRE

[phonons]
method = custom
ase_vibration_nfree = 2
phonon_displacement = 0.005
frequency_cutoff_cm = 1.0
symmetrize = true
fail_on_imaginary = true
use_hg_projected_modes = true
hg_mode_projection_threshold = 0.20

[thermodynamics]
method = custom
include_zpe = true
gas_reference_pressure = 1.0
calculate_temperature_scan = false

[runtime]
cuda_clear_between_phonons = true

[model]
path = /path/to/mace-mp-0b3-medium.model
```

The pressure is specified in bar in the project configuration. Internal
conversions to pascal are made where required by ASE or the statistical-
mechanical formulas.

## 16. Output files and provenance

The workflow writes files using the configured `output_prefix`.
Depending on the selected options and backend, the principal outputs
include:

- optimizer trajectory files when `save_trajectories = true`;
- `Au_phonons_frequencies_cm-1.dat`;
- `Au_phonons_signed_frequencies_cm-1.dat`;
- `Au_phonons_gamma_eigenvalues.dat`;
- `Au_phonons_mode_diagnostics.dat`;
- corresponding `HgAu_phonons_*` files;
- phonon metadata JSON files;
- `hg_au111_thermodynamics_scan.csv`; and
- `hg_au111_results.json`.

The frequency file without `signed` contains the positive modes retained
for thermochemistry. The signed-frequency file is intended for stability
and imaginary-mode diagnostics.

The final JSON combines the electronic adsorption energy, thermodynamic
corrections, backend-independent summary quantities, vibrational mode
counts, Hg-dominated frequencies, and the optional `T_cross` result.

## 17. Corrected constrained-phonon validation

The critical correction in the present workflow is that frozen Au atoms
are retained in every finite-difference force calculation. They are not
removed from the physical system. Only their coordinates are excluded
from the vibrational coordinate vector.

For the default 96-atom clean slab with 32 geometry-frozen atoms, force
evaluations still contain all 96 Au atoms, while the mobile-coordinate
Hessian is 192 x 192.

For Hg/Au, force evaluations contain all 97 atoms and the
mobile-coordinate Hessian is 195 x 195.

The expected metadata for the corrected custom finite-difference
calculation contains the conceptual flags

```
frozen_atoms_retained_in_force_calculations: true
```

and

```
method: central_finite_difference_full_system_mobile_hessian
```

These provide a direct provenance check that the old, incorrect
`remove-frozen-atoms` construction has not been used.

### 17.1 Historical phonon error and numerical consequences

An earlier version of the workflow removed the frozen Au atoms before
calculating displaced forces. That procedure changed the physical
mechanical response and produced a different vibrational spectrum.

The previously reported values from that incorrect calculation included:

- minimum frequency: about 10.7 cm<sup>-1</sup>;
- ZPE correction: about 0.00891 eV;
- $`\Delta G_{\text{ads}}(298\,\text{K})`$: about -0.16293 eV;
- $`\Delta G_{\text{ads}} = 0`$ crossover: about 448.887 K.

These values must not be mixed with results from the corrected
constrained finite-difference construction.

The corrected workflow must be rerun from the optimized structures. A
change relative to the historical values is expected and is not, by
itself, evidence of a regression: it reflects the correction of the
force-evaluation system.

### 17.2 Backend validation

The code provides two thermochemistry backends:

- `custom`: explicit project equations for harmonic vibrations and
  ideal-gas Hg;
- `ase`: ASE `HarmonicThermo` plus `IdealGasThermo`.

The supplied validation record reports that:

- all Python modules compile successfully with `python -m py_compile`;
- custom backend dispatch was exercised;
- the ASE backend was exercised against a local test double implementing
  the documented ASE thermochemistry equations;
- the resulting thermodynamic quantities agree with the custom backend
  to numerical precision;
- a full ASE/MACE calculation should nevertheless be run in the intended
  `mace_env` environment for end-to-end validation.

The implementation also checks internally that

``` math
\Delta G_{\text{ads}}
  = \Delta H_{\text{ads}} - T \Delta S_{\text{ads}}
```

and

``` math
\Delta G_{\text{ads}}
  = \Delta E_{\text{ads}} + \Delta F_{\text{vib}}
    - \mu_{\text{Hg,thermal}}
```

agree within numerical tolerance.

### 17.3 Explicit convergence protocol

Before treating $`\Delta G`$ or `T_cross` as publication-level
quantitative results, repeat the phonon calculation with at least

```
phonon_displacement = 0.003 Angstrom
phonon_displacement = 0.005 Angstrom
phonon_displacement = 0.008 Angstrom
```

and examine the stability of the frequencies, ZPE,
$`\Delta F_{\text{vib}}`$, $`\Delta G`$, and `T_cross`.

Also test:

- a larger surface cell than 4 x 4;
- increased slab thickness beyond six layers;
- the number of vibrationally frozen layers;
- adsorption sites fcc, hcp, bridge, and top;
- finite-temperature structural effects;
- explicit coverage dependence;
- and, where quantitative periodic vibrational free energies are
  required, Brillouin-zone/q-point sampling.

### 17.4 Expected mode counts and stability checks

For the default model:

```
clean Au: 64 mobile atoms -> 192 modes
Hg/Au:    65 mobile atoms -> 195 modes
```

The positive-frequency thermochemistry list is obtained only after the
signed spectrum has been examined.

With `fail_on_imaginary = true`, a significant imaginary mode terminates
the thermochemical stage rather than being silently discarded. A
low-frequency or imaginary result should instead trigger inspection of
the optimized geometry, displacement amplitude, slab model, frozen-layer
definition, and convergence.

Because the bottom Au layers are fixed, the constrained system does not
have the same translational zero modes as a freely translating isolated
system. Accordingly, the acoustic sum-rule correction is disabled.

## 18. Computational scope and limitations

### 18.1 Gamma-point approximation

The vibrational calculation is a Gamma-point calculation for the finite
surface supercell. It does not integrate phonon free energies over the
Brillouin zone. For a periodic surface, quantitative vibrational
free-energy convergence may require explicit q-point sampling or an
appropriate periodic phonon treatment.

### 18.2 Harmonic approximation

The thermochemistry is harmonic. Anharmonicity, temperature-dependent
structural relaxation, phonon-phonon interactions, and thermal expansion
are not included.

### 18.3 Frozen-layer approximation

The vibrational free energy corresponds to the selected constrained
coordinate space. Changing `n_frozen_layers_vibrations` changes the
vibrational model without changing the optimized geometry. This
parameter therefore represents a deliberate physical/modeling choice and
should be tested for convergence.

### 18.4 Surface and slab convergence

The 4x4x6 model should not automatically be regarded as converged.
Before using quantitative $`\Delta G`$ or `T_cross` values for
publication, test at least surface-cell size, slab thickness, number of
frozen layers, finite-difference displacement, and the treatment of
low-frequency modes.

### 18.5 Adsorption-site comparison

The supplied calculation uses an fcc starting site. Assigning fcc as the
thermodynamic equilibrium site requires comparison with competing sites
such as top, bridge, and hcp, with consistent geometry and vibrational
treatment.

### 18.6 Gas-reference pressure

The calculated $`\Delta G_{\text{ads}}`$ and `T_cross` are pressure
dependent because the gas entropy and chemical potential depend on Hg
pressure. Report the pressure explicitly with any thermodynamic result.

### 18.7 Comparison of vibrational constraint models

Two vibrational treatments were documented for the same electronic
adsorption energy. The purpose of this comparison is to show the effect
of the selected vibrational coordinate space, not to change the relaxed
geometry.

Standard treatment:

```
Delta ZPE             = 0.00881388 eV
Delta H_ads           = -0.48791444 eV
Delta S_ads           = -2.07727 x 10^-3 eV/K
Delta G_ads(298 K)    =  0.13142 eV
Delta G_ads = 0 T     = about 235 K
```

Frozen-slab vibrational treatment:

```
Delta ZPE             = 0.00637762 eV
Delta H_ads           = -0.48825772 eV
Delta S_ads           = -2.05303 x 10^-3 eV/K
Delta G_ads(298 K)    =  0.12385 eV
Delta G_ads = 0 T     = about 238 K
```

The electronic adsorption energy is unchanged between these treatments.
The difference comes from the modified vibrational spectrum and
therefore from $`\Delta F_{\text{vib}}`$, $`\Delta U_{\text{vib}}`$, and
$`\Delta S_{\text{vib}}`$.

## 19. Representative numerical result

For the documented custom-backend calculation with two geometry-frozen
Au layers and two vibrationally frozen Au layers, the representative
values at 298.15 K and 1e-5 bar were approximately:

```
Delta E_ads       = -0.50132288 eV
Delta ZPE         =  0.00881388 eV
Delta U_vib       =  0.06882601 eV
H_gas,thermal     =  0.06423145 eV
Delta H_ads       = -0.48791444 eV
Delta S_ads       = -2.077265e-3 eV/K
Delta F_vib       = -0.13949947 eV
mu_Hg,thermal     = -0.77224453 eV
Delta G_ads       =  0.13142218 eV
T_cross           =  235.025029 K
```

These numbers are model-specific and should be reproduced only with the
same structure, calculator/model, vibrational coordinate definition,
frequency cutoff, pressure, and thermochemistry conventions.

### 19.1 Explicit representative temperature scan

For the representative custom-backend calculation at 1e-5 bar, the
documented scan is:

<table style="width:94%;">
<colgroup>
<col style="width: 15%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads(eV)</th>
<th>Delta S_ads(eV/K)</th>
<th>Delta G_ads(eV)</th>
</tr>
</thead>
<tbody>
<tr>
<td>50.00</td>
<td>-0.49620164</td>
<td>-2.12731440e-3</td>
<td>-0.38983592</td>
</tr>
<tr>
<td>100.00</td>
<td>-0.49539024</td>
<td>-2.11729969e-3</td>
<td>-0.28366027</td>
</tr>
<tr>
<td>200.00</td>
<td>-0.49187176</td>
<td>-2.09333473e-3</td>
<td>-0.07320481</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48791444</td>
<td>-2.07726519e-3</td>
<td><blockquote>
<p>0.13142218</p>
</blockquote></td>
</tr>
<tr>
<td>400.00</td>
<td>-0.48366842</td>
<td>-2.06502003e-3</td>
<td><blockquote>
<p>0.34233959</p>
</blockquote></td>
</tr>
<tr>
<td>500.00</td>
<td>-0.47944343</td>
<td>-2.05559378e-3</td>
<td><blockquote>
<p>0.54835346</p>
</blockquote></td>
</tr>
<tr>
<td>600.00</td>
<td>-0.47519066</td>
<td>-2.04784062e-3</td>
<td><blockquote>
<p>0.75351372</p>
</blockquote></td>
</tr>
</tbody>
</table>

The corresponding thermodynamic scan is written to
`hg_au111_thermodynamics_scan.csv`.

Publication-oriented plots documented by the supplied supplementary
material are:

- `DeltaH_ads_vs_T.png`
- `DeltaS_ads_vs_T.png`
- `DeltaG_ads_vs_T.png`

These figures visualize the temperature dependence of adsorption
enthalpy, entropy, and Gibbs free energy.

### 19.2 Physical interpretation of the temperature dependence

The electronic interaction is favorable because
$`\Delta E_{\text{ads}} < 0`$.

Nevertheless, adsorption removes the large translational entropy
associated with a gas-phase Hg atom. Thus $`\Delta S_{\text{ads}} < 0`$,
and the term $`-T \Delta S_{\text{ads}}`$ becomes increasingly positive
as temperature rises.

The calculated $`\Delta G_{\text{ads}}`$ therefore becomes less
favorable with increasing temperature in the documented pressure range.
The crossover temperature is the point where the favorable electronic
and vibrational contributions are balanced by the gas-phase
thermodynamic penalty.

## 20. Recommended reproducibility protocol

For a publication-quality calculation:

1.  Record the exact MACE model file and software versions.
2.  Record the complete INI configuration.
3.  Save the optimized Au and Hg/Au structures.
4.  Save the signed and positive vibrational frequencies and metadata.
5.  Verify that significant imaginary modes are absent.
6.  Check convergence with respect to finite-difference displacement.
7.  Check convergence with respect to slab thickness and surface cell.
8.  Check convergence with respect to the number of vibrationally frozen
    Au layers.
9.  Repeat with both `custom` and `ase` backends as an implementation
    cross-check.
10. Report temperature, Hg pressure, gas reference pressure, and ZPE
    convention together with $`\Delta G`$ and `T_cross`.

The resulting $`\Delta G(T)`$ and `T_cross` should be regarded as
properties of the specified computational model, not as standalone
experimental observables.

### 20.1 Results from the supplied ASE/MACE production logs

The uploaded production logs provide a consistent set of electronic
energies and ASE vibrational/thermochemical results for several Hg
pressures. The electronic energies are pressure-independent:

```
E_Au     = -304.8265897028 eV
E_HgAu   = -305.4524593272 eV
E_Hg     = -0.1245467400 eV
```

and therefore $`\Delta E_{\text{ads}} = -0.5013228844`$ eV.

For the standard vibrational treatment (all 192 clean-Au mobile modes
and all 195 Hg/Au mobile modes retained), the spectra reported in the
logs contain no significant imaginary modes. The clean Au frequencies
span 4.309--140.294 cm<sup>-1</sup> and the Hg/Au frequencies span
4.281--140.090 cm<sup>-1</sup>. Two Hg-dominated Hg/Au modes are
identified at approximately 19.312 and 19.517 cm<sup>-1</sup> using the
reported projection threshold of 0.200.

At 298.15 K and 1e-5 bar, the standard treatment gives:

```
Delta ZPE          =  0.00881388 eV
Delta U_vib        =  0.06882601 eV
H_Hg,thermal       =  0.06423145 eV
Delta H_ads        = -0.48791444 eV
S_Hg               =  2.80555417e-3 eV/K
Delta S_vib        =  7.28288986e-4 eV/K
Delta S_ads        = -2.07726519e-3 eV/K
Delta F_vib        = -0.13949947 eV
mu_Hg,thermal      = -0.77224453 eV
Delta G_ads        =  0.13142218 eV
```

The corresponding zero of the Gibbs free energy is
$`T_{\text{cross}} = 235.025029`$ K, with the bisection residual
reported as
$`\Delta G_{\text{ads}}(T_{\text{cross}}) = -3.16 \times 10^{-8}`$ eV.

The same pressure-independent electronic adsorption energy is combined
with the pressure-dependent gas chemical potential. The
standard-vibrational ASE logs give the following pressure dependence at
298.15 K:

<table style="width:69%;">
<colgroup>
<col style="width: 19%" />
<col style="width: 26%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>P(Hg) [bar]</th>
<th>Delta G_ads [eV]</th>
<th>T_cross [K]</th>
</tr>
</thead>
<tbody>
<tr>
<td>1</td>
<td>-0.16437430</td>
<td>450.921106</td>
</tr>
<tr>
<td>1e-1</td>
<td>-0.10521497</td>
<td>380.457339</td>
</tr>
<tr>
<td>1e-2</td>
<td>-0.04605564</td>
<td>329.270472</td>
</tr>
<tr>
<td>1e-3</td>
<td><blockquote>
<p>0.01310369</p>
</blockquote></td>
<td>290.354684</td>
</tr>
<tr>
<td>1e-4</td>
<td><blockquote>
<p>0.07226302</p>
</blockquote></td>
<td>259.745055</td>
</tr>
<tr>
<td>1e-5</td>
<td><blockquote>
<p>0.13142218</p>
</blockquote></td>
<td>235.025029</td>
</tr>
</tbody>
</table>

This pressure series demonstrates directly the expected logarithmic
ideal-gas pressure dependence. At fixed temperature, decreasing Hg
pressure makes $`\Delta G_{\text{ads}}`$ less favorable. Over this
pressure interval, the crossover temperature correspondingly decreases
from about 451 K at 1 bar to 235 K at 1e-5 bar.

### 20.2 Standard-vibrational temperature scan at multiple pressures

For the standard vibrational treatment, $`\Delta H_{\text{ads}}`$ is
essentially pressure-independent because the pressure enters through the
gas entropy and chemical potential rather than through the electronic or
solid vibrational terms.

**P = 1 bar**

<table style="width:93%;">
<colgroup>
<col style="width: 13%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads(eV)</th>
<th>Delta S_ads(eV/K)</th>
<th>Delta G_ads(eV)</th>
</tr>
</thead>
<tbody>
<tr>
<td>50</td>
<td>-0.49620164</td>
<td>-1.13520814e-3</td>
<td>-0.43944123</td>
</tr>
<tr>
<td>100</td>
<td>-0.49539024</td>
<td>-1.12519345e-3</td>
<td>-0.38287090</td>
</tr>
<tr>
<td>200</td>
<td>-0.49187176</td>
<td>-1.10122850e-3</td>
<td>-0.27162606</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48791444</td>
<td>-1.08515897e-3</td>
<td>-0.16437430</td>
</tr>
<tr>
<td>400</td>
<td>-0.48366843</td>
<td>-1.07291381e-3</td>
<td>-0.05450290</td>
</tr>
<tr>
<td>500</td>
<td>-0.47944344</td>
<td>-1.06348757e-3</td>
<td><blockquote>
<p>0.05230035</p>
</blockquote></td>
</tr>
<tr>
<td>600</td>
<td>-0.47519067</td>
<td>-1.05573441e-3</td>
<td><blockquote>
<p>0.15824998</p>
</blockquote></td>
</tr>
</tbody>
</table>

**P = 1e-1 bar**

<table style="width:93%;">
<colgroup>
<col style="width: 13%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads(eV)</th>
<th>Delta S_ads(eV/K)</th>
<th>Delta G_ads(eV)</th>
</tr>
</thead>
<tbody>
<tr>
<td>50</td>
<td>-0.49620164</td>
<td>-1.33362951e-3</td>
<td>-0.42952016</td>
</tr>
<tr>
<td>100</td>
<td>-0.49539024</td>
<td>-1.32361482e-3</td>
<td>-0.36302876</td>
</tr>
<tr>
<td>200</td>
<td>-0.49187176</td>
<td>-1.29964987e-3</td>
<td>-0.23194179</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48791444</td>
<td>-1.28358033e-3</td>
<td>-0.10521497</td>
</tr>
<tr>
<td>400</td>
<td>-0.48366843</td>
<td>-1.27133518e-3</td>
<td><blockquote>
<p>0.02486564</p>
</blockquote></td>
</tr>
<tr>
<td>500</td>
<td>-0.47944344</td>
<td>-1.26190893e-3</td>
<td><blockquote>
<p>0.15151103</p>
</blockquote></td>
</tr>
<tr>
<td>600</td>
<td>-0.47519067</td>
<td>-1.25415577e-3</td>
<td><blockquote>
<p>0.27730280</p>
</blockquote></td>
</tr>
</tbody>
</table>

**P = 1e-2 bar**

<table style="width:93%;">
<colgroup>
<col style="width: 13%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads(eV)</th>
<th>Delta S_ads(eV/K)</th>
<th>Delta G_ads(eV)</th>
</tr>
</thead>
<tbody>
<tr>
<td>50</td>
<td>-0.49620164</td>
<td>-1.53205087e-3</td>
<td>-0.41959909</td>
</tr>
<tr>
<td>100</td>
<td>-0.49539024</td>
<td>-1.52203618e-3</td>
<td>-0.34318662</td>
</tr>
<tr>
<td>200</td>
<td>-0.49187176</td>
<td>-1.49807123e-3</td>
<td>-0.19225751</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48791444</td>
<td>-1.48200169e-3</td>
<td>-0.04605564</td>
</tr>
<tr>
<td>400</td>
<td>-0.48366843</td>
<td>-1.46975654e-3</td>
<td><blockquote>
<p>0.10423419</p>
</blockquote></td>
</tr>
<tr>
<td>500</td>
<td>-0.47944344</td>
<td>-1.46033030e-3</td>
<td><blockquote>
<p>0.25072171</p>
</blockquote></td>
</tr>
<tr>
<td>600</td>
<td>-0.47519067</td>
<td>-1.45257714e-3</td>
<td><blockquote>
<p>0.39635562</p>
</blockquote></td>
</tr>
</tbody>
</table>

**P = 1e-3 bar**

<table style="width:93%;">
<colgroup>
<col style="width: 13%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads(eV)</th>
<th>Delta S_ads(eV/K)</th>
<th>Delta G_ads(eV)</th>
</tr>
</thead>
<tbody>
<tr>
<td>50</td>
<td>-0.49620164</td>
<td>-1.73047223e-3</td>
<td>-0.40967802</td>
</tr>
<tr>
<td>100</td>
<td>-0.49539024</td>
<td>-1.72045755e-3</td>
<td>-0.32334449</td>
</tr>
<tr>
<td>200</td>
<td>-0.49187176</td>
<td>-1.69649259e-3</td>
<td>-0.15257324</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48791444</td>
<td>-1.68042306e-3</td>
<td><blockquote>
<p>0.01310369</p>
</blockquote></td>
</tr>
<tr>
<td>400</td>
<td>-0.48366843</td>
<td>-1.66817790e-3</td>
<td><blockquote>
<p>0.18360273</p>
</blockquote></td>
</tr>
<tr>
<td>500</td>
<td>-0.47944344</td>
<td>-1.65875166e-3</td>
<td><blockquote>
<p>0.34993239</p>
</blockquote></td>
</tr>
<tr>
<td>600</td>
<td>-0.47519067</td>
<td>-1.65099850e-3</td>
<td><blockquote>
<p>0.51540844</p>
</blockquote></td>
</tr>
</tbody>
</table>

**P = 1e-4 bar**

<table style="width:93%;">
<colgroup>
<col style="width: 13%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads(eV)</th>
<th>Delta S_ads(eV/K)</th>
<th>Delta G_ads(eV)</th>
</tr>
</thead>
<tbody>
<tr>
<td>50</td>
<td>-0.49620164</td>
<td>-1.92889360e-3</td>
<td>-0.39975696</td>
</tr>
<tr>
<td>100</td>
<td>-0.49539024</td>
<td>-1.91887891e-3</td>
<td>-0.30350235</td>
</tr>
<tr>
<td>200</td>
<td>-0.49187176</td>
<td>-1.89491396e-3</td>
<td>-0.11288897</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48791444</td>
<td>-1.87884442e-3</td>
<td><blockquote>
<p>0.07226302</p>
</blockquote></td>
</tr>
<tr>
<td>400</td>
<td>-0.48366843</td>
<td>-1.86659927e-3</td>
<td><blockquote>
<p>0.26297128</p>
</blockquote></td>
</tr>
<tr>
<td>500</td>
<td>-0.47944344</td>
<td>-1.85717302e-3</td>
<td><blockquote>
<p>0.44914307</p>
</blockquote></td>
</tr>
<tr>
<td>600</td>
<td>-0.47519067</td>
<td>-1.84941986e-3</td>
<td><blockquote>
<p>0.63446125</p>
</blockquote></td>
</tr>
</tbody>
</table>

### 20.3 Frozen-slab vibrational control from the supplied logs

A separate control calculation freezes all Au atoms for the vibrational
calculation while retaining Hg as mobile. The clean Au system therefore
has zero vibrational modes in this control, whereas Hg/Au retains only
the three Hg Cartesian vibrational modes. The Hg/Au control frequencies
are

```
22.2127, 22.3355, and 58.3297 cm^-1,
```

with no significant imaginary modes.

The resulting vibrational ZPE difference is
$`\Delta \text{ZPE} = 0.00637762`$ eV.

At 298.15 K and 1e-5 bar, the frozen-slab control gives:

```
Delta U_vib       =  0.07091897 eV
H_Hg,thermal      =  0.06423143 eV
Delta H_ads       = -0.48825772 eV
Delta S_vib       =  7.52526635e-4 eV/K
Delta S_ads       = -2.05302780e-3 eV/K
Delta F_vib       = -0.14706923 eV
mu_Hg,thermal     = -0.77224463 eV
Delta G_ads       =  0.12385252 eV
```

The corresponding crossover is $`T_{\text{cross}} = 237.957615`$ K, with
a reported residual of approximately $`1.41 \times 10^{-8}`$ eV.

At 1 bar the same frozen-slab vibrational control gives

```
Delta G_ads(298.15 K) = -0.17194413 eV
T_cross               =  461.759520 K
```

The frozen-slab control therefore shifts the 1e-5-bar crossover from
235.025 K to 237.958 K, while shifting the 1-bar crossover from 450.921
K to 461.760 K. This control changes only the vibrational coordinate
space; the electronic adsorption energy remains
$`\Delta E_{\text{ads}} = -0.5013228844`$ eV.

### 20.4 Frozen-slab temperature scan

**At 1e-5 bar:**

<table style="width:93%;">
<colgroup>
<col style="width: 13%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads(eV)</th>
<th>Delta S_ads(eV/K)</th>
<th>Delta G_ads(eV)</th>
</tr>
</thead>
<tbody>
<tr>
<td>50</td>
<td>-0.49790851</td>
<td>-2.11795637e-3</td>
<td>-0.39201069</td>
</tr>
<tr>
<td>100</td>
<td>-0.49636703</td>
<td>-2.09724285e-3</td>
<td>-0.28664274</td>
</tr>
<tr>
<td>200</td>
<td>-0.49237968</td>
<td>-2.06978456e-3</td>
<td>-0.07842277</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48825772</td>
<td>-2.05302780e-3</td>
<td><blockquote>
<p>0.12385252</p>
</blockquote></td>
</tr>
<tr>
<td>400</td>
<td>-0.48392500</td>
<td>-2.04052892e-3</td>
<td><blockquote>
<p>0.33228657</p>
</blockquote></td>
</tr>
<tr>
<td>500</td>
<td>-0.47964895</td>
<td>-2.03098780e-3</td>
<td><blockquote>
<p>0.53584495</p>
</blockquote></td>
</tr>
<tr>
<td>600</td>
<td>-0.47536204</td>
<td>-2.02317206e-3</td>
<td><blockquote>
<p>0.73854119</p>
</blockquote></td>
</tr>
</tbody>
</table>

**At 1 bar:**

<table style="width:93%;">
<colgroup>
<col style="width: 13%" />
<col style="width: 25%" />
<col style="width: 30%" />
<col style="width: 23%" />
</colgroup>
<thead>
<tr>
<th>T (K)</th>
<th>Delta H_ads(eV)</th>
<th>Delta S_ads(eV/K)</th>
<th>Delta G_ads(eV)</th>
</tr>
</thead>
<tbody>
<tr>
<td>50</td>
<td>-0.49790851</td>
<td>-1.12584956e-3</td>
<td>-0.44161603</td>
</tr>
<tr>
<td>100</td>
<td>-0.49636703</td>
<td>-1.10513604e-3</td>
<td>-0.38585342</td>
</tr>
<tr>
<td>200</td>
<td>-0.49237968</td>
<td>-1.07767775e-3</td>
<td>-0.27684413</td>
</tr>
<tr>
<td>298.15</td>
<td>-0.48825772</td>
<td>-1.06092098e-3</td>
<td>-0.17194413</td>
</tr>
<tr>
<td>400</td>
<td>-0.48392500</td>
<td>-1.04842210e-3</td>
<td>-0.06455616</td>
</tr>
<tr>
<td>500</td>
<td>-0.47964895</td>
<td>-1.03888098e-3</td>
<td><blockquote>
<p>0.03979154</p>
</blockquote></td>
</tr>
<tr>
<td>600</td>
<td>-0.47536204</td>
<td>-1.03106524e-3</td>
<td><blockquote>
<p>0.14327710</p>
</blockquote></td>
</tr>
</tbody>
</table>

The associated ASE vibrational thermochemistry output for the
frozen-slab control reports, at 298.15 K,

```
ZPE   = 0.0063776208 eV
F_vib = -0.1470692256 eV
U_vib = 0.0772965905 eV
S_vib = 7.52526635e-4 eV/K
```

These quantities refer to the Hg/Au control spectrum because the clean
Au vibrational coordinate space is empty when every Au atom is frozen.

### 20.5 Backend and log provenance

The supplied production logs explicitly identify

```
Vibrational backend: ASE Vibrations
```

and, for the thermochemical runs,

```
Thermodynamics backend: ASE HarmonicThermo + IdealGasThermo
```

Thus the numerical pressure series and the detailed values in sections
20.1--20.4 are ASE-backend production results. They should be
distinguished from the project-custom equations/backend discussed
elsewhere in this document.

The agreement of the electronic energies across all uploaded runs
provides a useful reproducibility check. Likewise, the identical
standard-vibrational spectra across the pressure series show that
changing Hg pressure affects the gas thermodynamic reference but does
not alter the calculated slab Hessians in these runs.

## 21. Current package status

The source package described by this document contains the modular
geometry, energy, phonon, thermochemistry, backend-dispatch, scan, and
Delta G=0 functionality described above.

The exact uploaded package used as the basis for this document does not
contain a separate `[restart]` section for reading previously saved
relaxed structures and vibrational data. Therefore this document does
not claim that file-based restart/post-processing is part of that exact
package. If a restart extension is added, its saved-structure and
saved-vibrational-data conventions should be documented separately and
validated against the metadata of the calculation that produced those
files.

## 22. Consolidated-source completeness statement

This document consolidates the applied theory, computational
methodology, validation information, numerical benchmarks, and
reproducibility instructions contained in the supplied project
documentation.

In particular, it preserves:

- the corrected full-system constrained finite-difference phonon
  construction;
- the Gamma-point-only scientific scope;
- the significant-imaginary-mode stopping policy;
- the disabled acoustic sum-rule correction for the constrained slab;
- the exact default slab/model/pressure/temperature parameters;
- the complete harmonic vibrational statistical mechanics;
- the complete monatomic ideal-gas Hg statistical mechanics;
- the pressure and reference-pressure dependence;
- the H-TS and F-mu free-energy formulations and their consistency
  check;
- the explicit temperature-scan values;
- the two documented vibrational-constraint comparisons;
- the historical incorrect-phonon numerical results and their cause;
- the backend validation record;
- the exact recommended finite-difference convergence points;
- the expected mode and force-evaluation counts;
- the adsorption-site validation recommendation;
- the publication-oriented output figure names; and
- the reproducibility and publication-scope limitations.

The document intentionally distinguishes representative numerical
results from converged physical predictions. Gamma-point harmonic
thermochemistry, finite-size slabs, frozen-layer constraints, and the
chosen MACE model define the computational model whose $`\Delta G`$ and
`T_cross` are reported.
