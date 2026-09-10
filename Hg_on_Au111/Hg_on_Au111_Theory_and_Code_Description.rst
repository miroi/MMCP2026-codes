Hg adsorption on Au(111): Applied Theory and Code Description
===============================================================

Purpose
-------

This document describes the physical model, equations, numerical procedures,
configuration controls, file outputs, and program architecture implemented in
the modular Hg-on-Au(111) workflow supplied with this project.  It is intended
as a stand-alone Supplementary-Information-style description of the calculation
rather than as a line-by-line programming manual.

Equation numbering and cross-references
----------------------------------------

All displayed mathematical equations are assigned a stable reference label and
a sequential equation number, shown as ``(Eq. N)``.  Within an RST/Sphinx
document, an equation can be referenced with ``:ref:`Eq. N <eq-NNN>```;
for example, ``:ref:`Eq. 1 <eq-001>```.  The explicit labels are intended to
remain stable even if surrounding prose is edited, while the displayed
equation number provides a convenient publication-style reference.

The workflow combines:

* construction of a periodic Au(111) slab with ASE;
* placement of one Hg atom at a selected adsorption site;
* constrained geometry optimization using a MACE interatomic potential;
* Gamma-point finite-displacement vibrational analysis;
* two selectable vibrational backends, ``custom`` and ``ase``;
* two selectable thermochemistry backends, ``custom`` and ``ase``;
* monatomic ideal-gas thermodynamics for the gas-phase Hg reference;
* temperature-dependent adsorption thermodynamics;
* a temperature scan; and
* a robust bisection search for the temperature at which
  ``Delta G_ads(T) = 0``.

The present implementation is a constrained-coordinate, Gamma-point harmonic
model.  It is not a Brillouin-zone-converged phonon calculation.


0. Complete computational specification
----------------------------------------

The following parameter set is the explicit computational specification
documented in the supplied supplementary-information files.  It should be
treated as the reference setup when reproducing the representative numerical
results unless a different configuration is stated explicitly.

* surface: Au(111);
* adsorbate: one Hg atom;
* surface construction: ASE ``fcc111``;
* surface cell: 4 x 4;
* slab thickness: 6 Au layers;
* vacuum thickness: 18 Angstrom;
* lattice constant: 4.08 Angstrom;
* adsorption site: fcc;
* initial Hg height: 3.0 Angstrom;
* electronic-structure/force model: MACE;
* MACE model: ``mace-mp-0b3-medium.model``;
* ASE version: 3.29.0;
* MACE version: 0.3.16;
* phonon method: finite-displacement Hessian;
* default finite-displacement amplitude: 0.005 Angstrom;
* thermodynamic frequency cutoff: 1 cm^-1;
* representative Hg pressure: 1e-5 bar;
* temperature scan: 50--600 K;
* representative scan temperatures: 50, 100, 200, 298.15, 400, 500,
  and 600 K.

The project configuration is expressed in bar for pressure, Angstrom for
geometric displacements, eV for energies, cm^-1 for vibrational frequencies,
and K for temperature.  Conversion to SI pressure is performed internally
where required.

The representative numerical values in this document are therefore
model-specific.  They are not intended to establish convergence of the
physical problem without the convergence studies described later.

1. Physical system and structural model
----------------------------------------

1.1 Au(111) slab
~~~~~~~~~~~~~~~~

The clean surface is generated with ASE's ``fcc111`` construction.  The
standard configuration supplied with the project uses

* surface cell: ``4 x 4``;
* slab thickness: 6 Au layers;
* lattice constant: 4.08 Angstrom;
* vacuum thickness: 18.0 Angstrom.

The resulting slab contains 96 Au atoms, i.e. 16 Au atoms per layer for six
layers.  Periodic boundary conditions are inherited from the ASE surface
construction and are retained during optimization and force evaluation.

The bottom Au layers can be constrained during geometry relaxation through
``n_frozen_layers``.  In the standard setup this value is 2.

1.2 Hg adsorption geometry
~~~~~~~~~~~~~~~~~~~~~~~~~~

One Hg atom is added to the optimized clean slab.  The adsorption site is
selected with ``adsorption_site``; the supplied configuration uses the fcc
site.  ``hg_height`` specifies the initial vertical height used when placing
Hg.  The standard initial height is 3.0 Angstrom.

The Hg/Au system therefore contains 97 atoms before vibrational coordinate
selection: 96 Au atoms plus one Hg atom.

1.3 Separation of relaxation and vibrational constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A central feature of the implementation is that the constraints used for
geometry relaxation are not automatically the constraints used to define the
vibrational coordinate space.

``n_frozen_layers`` controls the geometry optimization.  A separate parameter,
``n_frozen_layers_vibrations``, controls which bottom Au layers are excluded
from the vibrational coordinate space.  This second setting does not alter the
optimized structure.

The vibrational freezing can therefore be increased after optimization.  It
may be equal to the total number of Au layers.  Hg is not frozen by this Au
layer setting.

For example, with a six-layer slab and
``n_frozen_layers_vibrations = 6``:

* all 96 Au atoms are frozen in the vibrational coordinate space for clean Au;
* no Au vibrational coordinates remain for clean Au;
* the single Hg atom remains mobile for Hg/Au;
* the Hg/Au vibrational Hessian consequently contains 3 coordinates and 3
  modes.

The frozen atoms are not deleted from the physical system.  They remain in
all force evaluations.  Only the selected mobile coordinates are displaced and
included in the Hessian.

2. Electronic-energy model and adsorption energy
-------------------------------------------------

The atomic energies are evaluated with a MACE calculator.  The default model
is ``mace-mp-0b3-medium.model``.  The calculator is configured for CUDA when a
CUDA-capable PyTorch installation is available; otherwise CPU execution is
used.

The adsorption energy is defined as

.. _eq-008:

  Delta E_ads = E_HgAu - E_Au - E_Hg

**(Eq. 8)**

where

* ``E_HgAu`` is the potential energy of the optimized Hg/Au slab;
* ``E_Au`` is the potential energy of the optimized clean Au slab; and
* ``E_Hg`` is the potential energy of an isolated Hg atom evaluated with the
  same MACE calculator.

Negative ``Delta E_ads`` therefore denotes energetically favorable adsorption
at the electronic-energy level.

For the representative calculation documented in the project outputs,

  E_Au    = -304.8265897028 eV
  E_HgAu  = -305.4524593272 eV
  E_Hg    =   -0.1245467400 eV

and hence

  Delta E_ads = -0.5013228844 eV.

Using 96.48533212 kJ mol^-1 per eV gives approximately

  Delta E_ads = -48.3703 kJ mol^-1.

The electronic adsorption energy is kept separate from vibrational and gas
thermodynamic corrections throughout the implementation.

3. Geometry optimization
------------------------

The clean Au slab and Hg/Au system are optimized independently.  The
optimizer is configurable; the supplied configuration uses FIRE.

The force convergence thresholds are independently configurable:

* ``fmax_au`` for the clean slab;
* ``fmax_hg`` for the Hg/Au system.

The geometry constraints are applied to the specified bottom Au layers.  The
reported maximum force is evaluated on the mobile atoms.

The optimization stage produces optimized atomic structures and optional
optimizer trajectories.  The final potential energies are then used in the
adsorption-energy expression above.

The Hg--surface distance is reported before and after optimization, together
with its change.

4. Harmonic vibrational theory
------------------------------

4.1 Constrained Hessian
~~~~~~~~~~~~~~~~~~~~~~~~

Let the mobile vibrational Cartesian coordinates be collected into a vector
q of dimension 3N_mobile.  The harmonic potential around the optimized
structure is

.. _eq-001:

  V(q) = V_0 + 1/2 q^T H q

**(Eq. 1)**

where H is the mass-unweighted Cartesian Hessian in the selected coordinate
space.

The key constrained-system construction is:

* all atoms remain present in the force calculation;
* only mobile atoms are displaced;
* forces on the complete system are calculated after every displacement; and
* only the Hessian block associated with mobile coordinates is diagonalized.

This preserves the interaction of mobile atoms with the vibrationally frozen
atoms.

4.2 Finite differences
~~~~~~~~~~~~~~~~~~~~~~~

The custom backend constructs the Hessian from finite differences of forces.
For a central two-point scheme, the derivative is represented schematically by

.. _eq-002:

  H_ij = - [F_i(q_j + d) - F_i(q_j - d)] / (2d)

**(Eq. 2)**

where d is the configured displacement in Angstrom.  The minus sign appears
because forces are the negative gradient of the potential.

The implementation also supports the ASE finite-displacement route.  For the
ASE backend, ASE's ``Vibrations`` machinery is used with the explicitly
selected mobile atom indices.  The full ``Atoms`` object and its calculator
remain available to the force calculations.

The displacement is controlled by

  phonon_displacement = 0.005

in Angstrom.  For ASE, ``ase_vibration_nfree`` controls the finite-difference
scheme exposed by ASE.


For a four-point central finite-difference force derivative, when the selected
implementation uses ``nfree = 4``, the corresponding first derivative can be
formed from the forces at -2d, -d, +d, and +2d.  With the sign convention
F = -dV/dq, the Hessian is obtained from the negative force derivative.  The
two-point central formula remains

.. _eq-003:

  H_ij
    = -[F_i(q_j+d)-F_i(q_j-d)]/(2d).

**(Eq. 3)**

For a structure with N_mobile mobile atoms, the central two-point Hessian
requires

  2 * (3 N_mobile)

complete force evaluations, because each of the 3 N_mobile Cartesian
coordinates is displaced in both directions.

For the default 4 x 4 x 6 slab:

* clean Au: 64 mobile atoms -> 192 Hessian coordinates -> 384 force
  evaluations;
* Hg/Au: 65 mobile atoms -> 195 Hessian coordinates -> 390 force
  evaluations.

These counts refer to force evaluations on the full physical systems.  The
frozen atoms are retained in every such calculation.

4.3 Mass weighting and normal modes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

After construction of the Cartesian Hessian, the Hessian is mass weighted.
For Cartesian coordinates i and j associated with atoms a and b,

.. _eq-004:

  D_ij = H_ij / sqrt(m_a m_b)

**(Eq. 4)**

with the appropriate atomic masses in atomic mass units.

Diagonalization gives eigenvalues lambda_k and normal-mode eigenvectors.  The
vibrational angular frequencies are related to the eigenvalues by

.. _eq-005:

  omega_k = sqrt(lambda_k)

**(Eq. 5)**

for positive lambda_k.  Frequencies are converted to wavenumbers in cm^-1.

The code retains signed frequencies for diagnostics, so imaginary modes can be
identified explicitly.

4.4 Imaginary-mode handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The configured ``frequency_cutoff_cm`` determines which positive frequencies
are retained for thermochemistry.  The supplied default is 1.0 cm^-1.

Significant imaginary modes are not silently converted into thermodynamic
modes.  With

  fail_on_imaginary = true

the calculation stops when a significant imaginary instability is detected.
This is preferable to masking an unstable optimized structure by simply
removing negative frequencies.

4.5 Acoustic modes and constrained slabs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a freely translating isolated system, translational acoustic modes can be
expected near zero frequency.  Here, however, bottom Au atoms are constrained
and the vibrational problem is a constrained-coordinate problem.  Consequently
a global translation of all mobile atoms is not necessarily a zero mode.

The implementation therefore does not apply an acoustic sum-rule correction by
default.  This is particularly important when the bottom layers are frozen.

5. Hg mode projection
---------------------

For Hg/Au, the code identifies which normal modes are dominated by Hg motion.
The Hg atom index in the mobile-coordinate space is determined and the mode
projection onto the Hg Cartesian subspace is evaluated.

A mode is classified as Hg-dominated when its Hg projection is greater than or
equal to

  hg_mode_projection_threshold

with the supplied default value 0.20.

This diagnostic is not a replacement for the complete vibrational partition
function.  Thermochemistry uses the selected positive vibrational frequencies
of the complete mobile-coordinate system.  Hg-dominated frequencies are
reported separately to help interpret the low-frequency adsorbate motion.


6.0 Thermodynamic state-function construction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The calculation distinguishes the electronic potential-energy contribution
from thermal and vibrational contributions.  The adsorbed state is represented
by the optimized Hg/Au slab and its selected mobile-coordinate harmonic modes;
the reference state is the optimized clean Au slab plus one gas-phase Hg atom.

For the adsorption reaction

.. _eq-006:

  Hg(g,T,P) + Au(111) -> Hg*/Au(111),

**(Eq. 6)**

the working thermodynamic convention is

.. _eq-007:

  Delta G_ads(T,P)
    = Delta E_ads + Delta F_vib(T) - mu_Hg,thermal(T,P).

**(Eq. 7)**

The same quantity is evaluated independently through

.. _eq-009:

  Delta G_ads
    = Delta H_ads - T Delta S_ads.

**(Eq. 9)**

The implementation checks that these two expressions agree numerically.

The electronic term is

.. _eq-010:

  Delta E_ads
    = E_HgAu - E_Au - E_Hg.

**(Eq. 10)**

The solid-state vibrational free-energy difference is

.. _eq-011:

  Delta F_vib
    = F_vib,HgAu - F_vib,Au.

**(Eq. 11)**

The gas chemical potential is defined relative to the isolated-Hg electronic
energy already included in Delta E_ads; hence the gas term used here is the
thermal ideal-gas contribution rather than a second electronic atomic-energy
term.

6.0.1 Harmonic oscillator energy spectrum
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a normal mode of frequency nu,

.. _eq-012:

  epsilon = h nu = h c tilde{nu},

**(Eq. 12)**

where tilde{nu} is the wavenumber in cm^-1.

The dimensionless temperature parameter is

.. _eq-013:

  x = epsilon/(k_B T) = h nu/(k_B T).

**(Eq. 13)**

The quantum harmonic-oscillator energy levels are

  E_n = h nu (n + 1/2),    n = 0,1,2,...

The single-mode vibrational partition function is

.. _eq-014:

  q_vib = sum_n exp[-beta h nu (n+1/2)]
        = exp(-x/2)/(1-exp(-x)),

**(Eq. 14)**

with

  beta = 1/(k_B T).

The corresponding Helmholtz free energy is

.. _eq-015:

  F_vib,mode
    = -k_B T ln(q_vib)
    = h nu/2 + k_B T ln(1-exp(-x)).

**(Eq. 15)**

The mean energy is

.. _eq-016:

  U_vib,mode
    = - d ln(q_vib)/d beta
    = h nu/2 + h nu/[exp(x)-1].

**(Eq. 16)**

Thus the zero-point term and the finite-temperature thermal term are
naturally separated as

.. _eq-017:

  ZPE_mode = h nu/2

**(Eq. 17)**

and

.. _eq-018:

  U_vib,thermal,mode = h nu/[exp(x)-1].

**(Eq. 18)**

6.0.2 Vibrational entropy from the partition function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a mode,

.. _eq-019:

  S_vib,mode = (U_vib,mode - F_vib,mode)/T.

**(Eq. 19)**

Substitution gives

.. _eq-020:

  S_vib,mode
    = k_B [
        x/(exp(x)-1)
        - ln(1-exp(-x))
      ].

**(Eq. 20)**

For independent harmonic modes, the total quantities are sums over all
retained positive modes.

6.0.3 Heat capacity
^^^^^^^^^^^^^^^^^^

Differentiation of the thermal vibrational energy gives

.. _eq-021:

  C_V,vib,mode
    = d U_vib,thermal/dT
    = k_B x^2 exp(x)/(exp(x)-1)^2.

**(Eq. 21)**

The total harmonic vibrational heat capacity is the sum over modes.

The implementation does not require C_V explicitly for the adsorption free
energy, but this relation is useful for interpreting the temperature
dependence of the vibrational internal energy and enthalpy.

6.0.4 Thermodynamic identities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For each harmonic state,

.. _eq-022:

  F_vib = U_vib,total - T S_vib,

**(Eq. 22)**

and therefore

.. _eq-023:

  Delta F_vib
    = Delta U_vib,total - T Delta S_vib.

**(Eq. 23)**

When ZPE is included,

.. _eq-024:

  Delta U_vib,total
    = Delta ZPE + Delta U_vib,thermal,

**(Eq. 24)**

so that

.. _eq-025:

  Delta F_vib
    = Delta ZPE
      + Delta U_vib,thermal
      - T Delta S_vib.

**(Eq. 25)**

This identity is the algebraic bridge between the two adsorption-free-energy
routes implemented in the code.

6. Vibrational statistical mechanics
-------------------------------------

For each retained positive frequency nu in cm^-1, the quantum of vibrational
energy is

.. _eq-026:

  epsilon = h c nu.

**(Eq. 26)**

In electron-volt units the implementation uses

.. _eq-027:

  epsilon(eV) = 1.2398419843320026e-4 * nu(cm^-1).

**(Eq. 27)**

6.1 Zero-point energy
~~~~~~~~~~~~~~~~~~~~~

For a set of harmonic modes,

  ZPE = 1/2 sum_k epsilon_k.

The adsorption ZPE correction is

.. _eq-028:

  Delta ZPE = ZPE_HgAu - ZPE_Au.

**(Eq. 28)**

The code retains the raw ZPE difference and applies the ``include_zpe`` switch
explicitly.

6.2 Thermal vibrational internal energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The thermal excitation energy of one harmonic mode is

.. _eq-029:

  U_k,thermal = epsilon_k / [ exp(epsilon_k/(k_B T)) - 1 ].

**(Eq. 29)**

Therefore

.. _eq-030:

  U_vib,thermal(T) = sum_k U_k,thermal.

**(Eq. 30)**

The adsorption vibrational thermal-energy correction is

.. _eq-031:

  Delta U_vib,thermal = U_HgAu,thermal - U_Au,thermal.

**(Eq. 31)**

If ZPE is included, the total vibrational internal energy is

.. _eq-032:

  U_vib,total = ZPE + U_vib,thermal.

**(Eq. 32)**

6.3 Vibrational Helmholtz free energy
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For a harmonic oscillator, including zero-point energy,

  F_k = 1/2 epsilon_k + k_B T ln[1 - exp(-epsilon_k/(k_B T))].

Thus

.. _eq-033:

  F_vib(T) = sum_k F_k

**(Eq. 33)**

and

.. _eq-034:

  Delta F_vib = F_HgAu,vib - F_Au,vib.

**(Eq. 34)**

When ZPE is disabled, the zero-point term is removed consistently from the
reported vibrational free-energy correction.

6.4 Vibrational entropy
~~~~~~~~~~~~~~~~~~~~~~~

The vibrational entropy can be obtained from

.. _eq-035:

  S_vib = (U_vib,total - F_vib) / T.

**(Eq. 35)**

The adsorption vibrational entropy change is

.. _eq-036:

  Delta S_vib = S_HgAu,vib - S_Au,vib.

**(Eq. 36)**

The implementation evaluates this quantity from the same harmonic mode sets
used for the vibrational energy and free energy.


6.5 Complete harmonic-thermodynamic formula set
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section collects the complete set of thermodynamic equations used by the
workflow in a form suitable for reproducing the implementation independently
of the Python code.  The surface and adsorbate are treated as harmonic
vibrational systems, while gas-phase Hg is treated as a monatomic ideal gas.
All thermodynamic quantities below are per adsorbed Hg atom unless explicitly
stated otherwise.

6.5.1 Energy units and frequency conversion
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code uses eV for energies and eV/K for entropies.  Vibrational frequencies
are stored as wavenumbers in cm^-1.  The photon-like quantum associated with a
vibrational wavenumber nu_bar is

.. _eq-037:

  epsilon = h c nu_bar.

**(Eq. 37)**

In the project's eV units,

.. _eq-038:

  epsilon [eV] = (h c)_eV_cm * nu_bar [cm^-1],

**(Eq. 38)**

where ``HC_EV_CM`` is the corresponding conversion constant.  The dimensionless
harmonic-oscillator variable is

.. _eq-039:

  x = epsilon / (k_B T).

**(Eq. 39)**

The same x is used in the energy, entropy, and free-energy expressions.

6.5.2 One harmonic vibrational mode
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For one mode of quantum energy epsilon = h nu, the harmonic-oscillator energy
levels are

  E_n = (n + 1/2) epsilon,    n = 0, 1, 2, ... .

The canonical partition function is

  q_vib = sum_{n=0}^infinity exp[-beta (n+1/2) epsilon]
        = exp(-beta epsilon/2) / [1 - exp(-beta epsilon)],

with

  beta = 1/(k_B T).

For a set of independent modes, the total vibrational partition function is

the product

.. _eq-040:

  Q_vib = product_k q_k.

**(Eq. 40)**

Equivalently,

.. _eq-041:

  ln Q_vib = sum_k ln q_k.

**(Eq. 41)**

6.5.3 Vibrational zero-point energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The zero-point energy of one harmonic mode is

  ZPE_k = 1/2 epsilon_k.

For all retained positive modes,

  ZPE = 1/2 sum_k epsilon_k
      = 1/2 sum_k h c nu_bar_k.

For adsorption, the vibrational zero-point correction is the difference
between the adsorbed system and the clean surface:

.. _eq-042:

  Delta ZPE = ZPE_HgAu - ZPE_Au.

**(Eq. 42)**

If ``include_zpe = false``, the implementation sets this adsorption correction
to zero in the thermochemical combination.  The raw ZPE is still available as
a diagnostic quantity, and the ASE free-energy output is adjusted so that the
switch has the same meaning in both backends.

6.5.4 Thermal vibrational internal energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For one harmonic mode, excluding its zero-point contribution,

.. _eq-043:

  U_vib,thermal,k
      = epsilon_k / [exp(epsilon_k/(k_B T)) - 1]
      = epsilon_k / [exp(x_k) - 1].

**(Eq. 43)**

Thus

.. _eq-044:

  U_vib,thermal(T) = sum_k epsilon_k/[exp(x_k)-1].

**(Eq. 44)**

The total harmonic vibrational internal energy, when ZPE is included, is

.. _eq-045:

  U_vib,total(T) = ZPE + U_vib,thermal(T).

**(Eq. 45)**

For an adsorption process,

.. _eq-046:

  Delta U_vib,thermal
      = U_vib,thermal,HgAu - U_vib,thermal,Au,

**(Eq. 46)**

and

.. _eq-047:

  Delta U_vib,total
      = Delta ZPE + Delta U_vib,thermal

**(Eq. 47)**

when ZPE is enabled.

The implementation evaluates the Bose-Einstein denominator numerically as
``expm1(x)``.  For very large x, the thermal contribution is numerically set to
zero once it is below floating-point relevance.

6.5.5 Vibrational Helmholtz free energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For one harmonic mode,

.. _eq-048:

  F_vib,k = -k_B T ln(q_k).

**(Eq. 48)**

Using the harmonic partition function gives

.. _eq-049:

  F_vib,k
      = 1/2 epsilon_k
        + k_B T ln[1 - exp(-epsilon_k/(k_B T))].

**(Eq. 49)**

Therefore,

.. _eq-050:

  F_vib(T) = sum_k { 1/2 epsilon_k
                   + k_B T ln[1-exp(-x_k)] }.

**(Eq. 50)**

The free energy can also be written as

.. _eq-051:

  F_vib = ZPE + k_B T sum_k ln[1-exp(-x_k)].

**(Eq. 51)**

The adsorption vibrational free-energy correction is

.. _eq-052:

  Delta F_vib(T) = F_vib,HgAu(T) - F_vib,Au(T).

**(Eq. 52)**

If ZPE is excluded by the project option, the corresponding ZPE term is
removed consistently from this expression.  With ZPE included, the identity

.. _eq-053:

  F_vib = U_vib,total - T S_vib

**(Eq. 53)**

holds for the harmonic oscillator.

6.5.6 Vibrational entropy
^^^^^^^^^^^^^^^^^^^^^^^^^

For one harmonic mode, the entropy is

.. _eq-054:

  S_vib,k
    = k_B [ x_k/(exp(x_k)-1)
            - ln(1-exp(-x_k)) ].

**(Eq. 54)**

Hence

.. _eq-055:

  S_vib(T) = k_B sum_k
             [ x_k/(exp(x_k)-1)
               - ln(1-exp(-x_k)) ].

**(Eq. 55)**

The equivalent thermodynamic identity used as a check is

.. _eq-056:

  S_vib = [U_vib,total - F_vib]/T.

**(Eq. 56)**

Because the zero-point term cancels between U and F, the same entropy results
whether U is written as ZPE plus thermal energy or in the corresponding
partition-function form.

The adsorption vibrational entropy change is

.. _eq-057:

  Delta S_vib(T) = S_vib,HgAu(T) - S_vib,Au(T).

**(Eq. 57)**

At fixed positive frequency, the low-temperature limit is S_vib -> 0.  For a
very soft mode, however, the harmonic entropy becomes large and therefore the
numerical treatment of low frequencies is particularly important.

6.5.7 Vibrational heat capacity (derived quantity)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Although heat capacity is not required as an independent input to the final
Delta G calculation, it follows from the same harmonic model and is useful for
interpreting the temperature dependence:

.. _eq-058:

  C_V,vib = (dU_vib,total/dT)_V
          = k_B sum_k x_k^2 exp(x_k)/[exp(x_k)-1]^2.

**(Eq. 58)**

The zero-point term has zero temperature derivative and therefore does not
contribute to C_V.

6.5.8 Thermodynamic identities used by the implementation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For the harmonic surface model,

.. _eq-059:

  F_vib = U_vib,total - T S_vib,

**(Eq. 59)**

and therefore

.. _eq-060:

  Delta F_vib = Delta U_vib,total - T Delta S_vib.

**(Eq. 60)**

When ZPE is included,

.. _eq-061:

  Delta U_vib,total = Delta ZPE + Delta U_vib,thermal.

**(Eq. 61)**

Consequently,

.. _eq-062:

  Delta F_vib
    = Delta ZPE + Delta U_vib,thermal - T Delta S_vib.

**(Eq. 62)**

This identity is the algebraic bridge between the explicit ``Delta H - T
Delta S`` route and the compact ``Delta E + Delta F - mu`` route used in the
code.

6.5.9 Clean-surface and adsorbed-system definitions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The two solid-state vibrational states are treated separately:

  Au : frequencies {nu_k^Au}

  HgAu : frequencies {nu_k^HgAu}.

For any vibrational function X,

.. _eq-063:

  Delta X_vib = X_vib,HgAu - X_vib,Au.

**(Eq. 63)**

Thus, explicitly,

.. _eq-064:

  Delta ZPE
    = 1/2 sum_k epsilon_k^HgAu
      - 1/2 sum_j epsilon_j^Au,

**(Eq. 64)**

.. _eq-065:

  Delta U_vib,thermal(T)
    = sum_k epsilon_k^HgAu/[exp(epsilon_k^HgAu/(k_B T))-1]
      - sum_j epsilon_j^Au/[exp(epsilon_j^Au/(k_B T))-1],

**(Eq. 65)**

.. _eq-066:

  Delta F_vib(T)
    = sum_k [epsilon_k^HgAu/2
             + k_B T ln(1-exp(-epsilon_k^HgAu/(k_B T)))]
      - sum_j [epsilon_j^Au/2
             + k_B T ln(1-exp(-epsilon_j^Au/(k_B T)))],

**(Eq. 66)**

.. _eq-067:

  Delta S_vib(T)
    = k_B sum_k [x_k^HgAu/(exp(x_k^HgAu)-1)
                  - ln(1-exp(-x_k^HgAu))]
      - k_B sum_j [x_j^Au/(exp(x_j^Au)-1)
                    - ln(1-exp(-x_j^Au))].

**(Eq. 67)**

Only frequencies passing the positive-frequency selection are included in
these thermodynamic sums.


7.0 Complete ideal-gas Hg statistical mechanics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gas reference is one monatomic Hg atom.  There are no molecular rotational
or vibrational degrees of freedom.  The thermal gas model therefore contains
translational motion only.

7.0.1 Translational partition function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For N identical monatomic particles in volume V, the translational partition
function is

.. _eq-068:

  Q_trans = q_trans^N / N!,

**(Eq. 68)**

where the one-particle partition function is

.. _eq-069:

  q_trans
    = V / Lambda^3
    = V (2 pi m k_B T / h^2)^(3/2),

**(Eq. 69)**

and the thermal de Broglie wavelength is

.. _eq-070:

  Lambda = h / sqrt(2 pi m k_B T).

**(Eq. 70)**

The Hg atomic mass used by the implementation is 200.59 u, converted to kg
per atom where required.

7.0.2 Ideal-gas entropy
^^^^^^^^^^^^^^^^^^^^^^^

Using Stirling's approximation for N! gives the Sackur--Tetrode form

.. _eq-071:

  S_trans
    = N k_B [
        ln(V/(N Lambda^3)) + 5/2
      ].

**(Eq. 71)**

Per Hg atom,

.. _eq-072:

  S_Hg
    = k_B [
        ln(V/(N Lambda^3)) + 5/2
      ].

**(Eq. 72)**

Using the ideal-gas equation

.. _eq-073:

  P V = N k_B T

**(Eq. 73)**

gives

.. _eq-074:

  S_Hg(T,P)
    = k_B [
        ln(k_B T/(P Lambda^3)) + 5/2
      ].

**(Eq. 74)**

Relative to a reference pressure P0,

.. _eq-075:

  S_Hg(T,P)
    = S_Hg(T,P0) - k_B ln(P/P0).

**(Eq. 75)**

Thus lowering the gas pressure increases the gas entropy and makes adsorption
less favorable through the chemical-potential term.

7.0.3 Ideal-gas internal energy and enthalpy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

A monatomic ideal gas has three translational quadratic degrees of freedom.
Therefore, per atom,

  U_gas = 3/2 k_B T.

The ideal-gas equation gives

.. _eq-076:

  P V = k_B T

**(Eq. 76)**

per atom, and hence

.. _eq-077:

  H_gas = U_gas + P V
        = 5/2 k_B T.

**(Eq. 77)**

The implementation calls this quantity ``H_gas_thermal``.

7.0.4 Gas chemical potential
^^^^^^^^^^^^^^^^^^^^^^^^^^^

The thermal chemical-potential contribution is evaluated as

.. _eq-078:

  mu_Hg,thermal(T,P)
    = H_gas,thermal(T) - T S_Hg(T,P).

**(Eq. 78)**

Equivalently,

.. _eq-079:

  mu_Hg,thermal
    = -k_B T ln(q_trans/N)

**(Eq. 79)**

for the corresponding classical ideal-gas reference, with the same pressure,
temperature, and standard-state convention.

The pressure dependence can be written as

.. _eq-080:

  mu_Hg,thermal(T,P)
    = mu_Hg,thermal(T,P0)
      + k_B T ln(P/P0).

**(Eq. 80)**

Consequently, increasing Hg pressure makes the gas chemical potential less
negative and adsorption thermodynamically more favorable.

7.0.5 Unit conversion
^^^^^^^^^^^^^^^^^^^^^

The project uses

.. _eq-081:

  1 bar = 10^5 Pa

**(Eq. 81)**

and

.. _eq-082:

  1 eV = 1.602176634 x 10^-19 J.

**(Eq. 82)**

For vibrational wavenumbers,

.. _eq-083:

  h c = 1.2398419843320026 x 10^-4 eV cm,

**(Eq. 83)**

so

.. _eq-084:

  epsilon(eV)
    = 1.2398419843320026 x 10^-4
      * tilde{nu}(cm^-1).

**(Eq. 84)**

7. Gas-phase Hg thermodynamics
------------------------------

The reference state for Hg is a monatomic ideal gas.  The gas reference is not
represented by a harmonic molecular vibration calculation because an isolated
atom has no molecular vibrational or rotational degrees of freedom.

7.1 Ideal-gas PV contribution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For one mole of ideal gas,

.. _eq-085:

  PV = RT.

**(Eq. 85)**

Per particle,

.. _eq-086:

  PV = k_B T.

**(Eq. 86)**

Therefore the thermal enthalpy contribution of a monatomic ideal gas is

  H_gas,thermal = 5/2 k_B T.

This term is subtracted from the adsorbed-state enthalpy correction because the
gas reference contains translational PV work whereas the adsorbed atom does
not.

7.2 Translational entropy
~~~~~~~~~~~~~~~~~~~~~~~~~

The translational partition function is

.. _eq-087:

  q_trans = (2 pi m k_B T / h^2)^(3/2) V.

**(Eq. 87)**

The Sackur--Tetrode expression used by the custom backend can be written as

.. _eq-088:

  S_trans = k_B [ ln(q_trans/N) + 5/2 ].

**(Eq. 88)**

Using the ideal-gas relation

.. _eq-089:

  V = N k_B T / P

**(Eq. 89)**

introduces the pressure dependence.  Relative to a reference pressure P0,

.. _eq-090:

  S(T,P) = S(T,P0) - k_B ln(P/P0).

**(Eq. 90)**

The gas chemical potential contribution is represented thermodynamically as

.. _eq-091:

  mu_Hg,thermal(T,P) = H_gas,thermal - T S_gas(T,P).

**(Eq. 91)**

The pressure used by the project is entered in bar and converted internally to
SI pressure where required.


7.3 Complete ideal-gas Hg statistical mechanics
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The gas reference is a single monatomic Hg atom.  Internal molecular
rotational and vibrational contributions are absent because Hg is monatomic.
The model therefore contains translational motion only.

7.3.1 Translational partition function
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For one classical particle in a volume V, the translational partition function
is

.. _eq-092:

  q_trans = V / Lambda^3,

**(Eq. 92)**

where the thermal de Broglie wavelength is

.. _eq-093:

  Lambda = h / sqrt(2 pi m k_B T).

**(Eq. 93)**

For N indistinguishable particles,

.. _eq-094:

  Q_trans = q_trans^N / N!.

**(Eq. 94)**

Using Stirling's approximation for large N and the thermodynamic limit gives
the standard ideal-gas expressions used by the Sackur--Tetrode formulation.

For a monatomic ideal gas,

  U_trans = 3/2 N k_B T,

and

  H_trans = U_trans + P V = 5/2 N k_B T,

because PV = N k_B T.

Per Hg atom, the thermal enthalpy used by the code is therefore

  H_gas,thermal(T) = 5/2 k_B T.

The isolated-Hg electronic energy is already included separately through
Delta E_ads; consequently this gas correction contains only the thermal
ideal-gas contribution.

7.3.2 Sackur--Tetrode entropy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For a monatomic ideal gas, the entropy per particle can be written

.. _eq-095:

  S/N = k_B [ ln( V/N Lambda^3 ) + 5/2 ].

**(Eq. 95)**

Using the ideal-gas equation

.. _eq-096:

  P V = N k_B T

**(Eq. 96)**

gives

.. _eq-097:

  V/N = k_B T/P,

**(Eq. 97)**

and hence

.. _eq-098:

  S_gas(T,P)
    = k_B [ ln( k_B T/(P Lambda^3) ) + 5/2 ].

**(Eq. 98)**

The pressure dependence at fixed T is therefore

.. _eq-099:

  S_gas(T,P)
    = S_gas(T,P0) - k_B ln(P/P0),

**(Eq. 99)**

where P0 is the configured reference pressure.

This is the explicit pressure correction implemented by the custom backend.
The pressure is supplied in bar in the input file and converted to pascal for
this SI statistical-mechanical expression.

7.3.3 Chemical potential
^^^^^^^^^^^^^^^^^^^^^^^^

For an ideal gas, the chemical potential can be written as

.. _eq-100:

  mu = -k_B T ln(q_trans/N)

**(Eq. 100)**

up to the standard indistinguishability formulation, or equivalently through
thermodynamic identities.  The implementation uses the enthalpy-entropy form
for the thermal contribution:

.. _eq-101:

  mu_Hg,thermal(T,P)
      = H_gas,thermal(T) - T S_gas(T,P).

**(Eq. 101)**

Substitution of the monatomic ideal-gas enthalpy gives

.. _eq-102:

  mu_Hg,thermal(T,P)
      = 5/2 k_B T - T S_gas(T,P).

**(Eq. 102)**

Because

.. _eq-103:

  S_gas(T,P) = S_gas(T,P0) - k_B ln(P/P0),

**(Eq. 103)**

one obtains the explicit pressure dependence

.. _eq-104:

  mu_Hg,thermal(T,P)
      = mu_Hg,thermal(T,P0) + k_B T ln(P/P0).

**(Eq. 104)**

Thus increasing the gas pressure increases the gas chemical potential.  Since
the adsorption free energy contains ``-mu_Hg,thermal``, increasing Hg pressure
makes adsorption thermodynamically more favorable in this convention.

7.3.4 Standard pressure versus actual pressure
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The code distinguishes the actual Hg pressure P from a reference pressure
P0.  The standard-pressure gas entropy is

.. _eq-105:

  S_gas,standard(T) = S_gas(T,P0),

**(Eq. 105)**

and the actual-pressure entropy is

.. _eq-106:

  S_gas(T,P) = S_gas,standard(T) - k_B ln(P/P0).

**(Eq. 106)**

The thermodynamic result therefore depends on both the configured pressure
and the configured reference pressure.  The latter fixes the zero of the gas
chemical-potential convention; the former represents the physical gas
condition being evaluated.

7.3.5 Pressure derivative of the adsorption free energy
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

At fixed temperature, the ideal-gas contribution implies

.. _eq-107:

  d mu_Hg,thermal / d ln P = k_B T.

**(Eq. 107)**

Therefore, for the adsorption free energy

.. _eq-108:

  Delta G_ads = Delta E_ads + Delta F_vib - mu_Hg,thermal,

**(Eq. 108)**

one has

.. _eq-109:

  d Delta G_ads / d ln P = -k_B T.

**(Eq. 109)**

This relation is a useful analytical check on pressure-dependent scans.

8. Adsorption enthalpy, entropy and Gibbs free energy
------------------------------------------------------

8.1 Reaction definition and thermodynamic state functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The modeled adsorption reaction is

.. _eq-110:

  Hg(g) + Au(surface) -> Hg* / Au(surface),

**(Eq. 110)**

where the asterisk denotes the adsorbed Hg atom in the selected optimized
adsorption state.  The electronic-energy contribution is

.. _eq-111:

  Delta E_ads = E_HgAu - E_Au - E_Hg.

**(Eq. 111)**

The thermodynamic adsorption free energy is constructed by adding the
vibrational correction for the solid states and subtracting the gas-phase Hg
chemical-potential contribution.

The central state-function definitions are

.. _eq-112:

  H = U + P V,

**(Eq. 112)**

.. _eq-113:

  F = U - T S,

**(Eq. 113)**

.. _eq-114:

  G = H - T S = F + P V,

**(Eq. 114)**

and for a gas-phase species the relevant reservoir quantity is its chemical
potential ``mu`` (Gibbs free energy per particle in the ideal-gas limit).

For the adsorption reaction, the code uses

.. _eq-115:

  Delta H_ads = H_HgAu - H_Au - H_Hg(g),

**(Eq. 115)**

.. _eq-116:

  Delta S_ads = S_HgAu - S_Au - S_Hg(g),

**(Eq. 116)**

.. _eq-117:

  Delta G_ads = Delta H_ads - T Delta S_ads.

**(Eq. 117)**

Because the electronic energy is separated from the thermal gas contribution,
these become the explicit equations given below.

8.2 Explicit enthalpy equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The surface and adsorbate electronic energies are supplied by MACE.  The
vibrational correction to their internal energies is harmonic.  The gas phase
contributes its thermal ideal-gas enthalpy.  Thus

.. _eq-118:

  Delta H_ads(T)
    = Delta E_ads
      + Delta ZPE
      + Delta U_vib,thermal(T)
      - H_gas,thermal(T),

**(Eq. 118)**

with

  H_gas,thermal(T) = 5/2 k_B T.

If ZPE is disabled,

.. _eq-119:

  Delta H_ads(T)
    = Delta E_ads
      + Delta U_vib,thermal(T)
      - H_gas,thermal(T).

**(Eq. 119)**

8.3 Explicit entropy equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The solid-state entropy change is

.. _eq-120:

  Delta S_vib(T) = S_vib,HgAu(T) - S_vib,Au(T).

**(Eq. 120)**

The reaction entropy is

.. _eq-121:

  Delta S_ads(T,P)
    = Delta S_vib(T) - S_gas(T,P).

**(Eq. 121)**

With the ideal-gas expression,

.. _eq-122:

  Delta S_ads(T,P)
    = Delta S_vib(T)
      - k_B [ ln(k_B T/(P Lambda^3)) + 5/2 ].

**(Eq. 122)**

Equivalently, relative to P0,

.. _eq-123:

  Delta S_ads(T,P)
    = Delta S_ads(T,P0) + k_B ln(P/P0).

**(Eq. 123)**

8.4 Explicit Gibbs-energy equation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Combining the preceding expressions gives

.. _eq-124:

  Delta G_ads(T,P)
    = Delta E_ads
      + Delta ZPE
      + Delta U_vib,thermal(T)
      - H_gas,thermal(T)
      - T Delta S_vib(T)
      + T S_gas(T,P).

**(Eq. 124)**

Using

.. _eq-125:

  Delta F_vib
    = Delta ZPE + Delta U_vib,thermal - T Delta S_vib

**(Eq. 125)**

and

.. _eq-126:

  mu_Hg,thermal = H_gas,thermal - T S_gas,

**(Eq. 126)**

this reduces exactly to

.. _eq-127:

  Delta G_ads(T,P)
    = Delta E_ads + Delta F_vib(T) - mu_Hg,thermal(T,P).

**(Eq. 127)**

This is the compact expression used for the internal cross-check.

8.5 Separation into electronic, vibrational and gas terms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is useful to display the free energy as

.. _eq-128:

  Delta G_ads
    = Delta E_ads
      + [Delta ZPE + Delta U_vib,thermal - T Delta S_vib]
      - [H_gas,thermal - T S_gas].

**(Eq. 128)**

The three physically distinct contributions are therefore:

* electronic adsorption energy: ``Delta E_ads``;
* surface/adsorbate vibrational free-energy correction: ``Delta F_vib``; and
* gas-phase Hg thermal chemical potential: ``mu_Hg,thermal``.

No gas-phase Hg electronic energy is added a second time: ``E_Hg`` is already
part of ``Delta E_ads``.

8.6 Pressure dependence of Delta G_ads
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

At fixed T, the pressure dependence follows directly from the ideal-gas
chemical potential:

.. _eq-129:

  mu_Hg(T,P) = mu_Hg(T,P0) + k_B T ln(P/P0).

**(Eq. 129)**

Therefore

.. _eq-130:

  Delta G_ads(T,P)
    = Delta G_ads(T,P0) - k_B T ln(P/P0).

**(Eq. 130)**

This is the analytical pressure dependence of the present model.

8.7 Alternative free-energy bookkeeping
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The code's preferred compact route is

.. _eq-131:

  Delta G_ads = Delta E_ads + Delta F_vib - mu_Hg,thermal.

**(Eq. 131)**

The expanded route is

.. _eq-132:

  Delta G_ads = Delta H_ads - T Delta S_ads.

**(Eq. 132)**

The implementation evaluates both and raises an error if they disagree beyond
the configured numerical tolerance.  Agreement is expected because

.. _eq-133:

  Delta F_vib = Delta ZPE + Delta U_vib,thermal - T Delta S_vib

**(Eq. 133)**

and

.. _eq-134:

  mu_Hg,thermal = H_gas,thermal - T S_gas.

**(Eq. 134)**

8.8 ZPE switch semantics
~~~~~~~~~~~~~~~~~~~~~~~~

When ``include_zpe = true``:

.. _eq-135:

  Delta H_ads = Delta E_ads + Delta ZPE + Delta U_vib,thermal - H_gas,

**(Eq. 135)**

.. _eq-136:

  Delta F_vib = Delta ZPE + Delta U_vib,thermal - T Delta S_vib.

**(Eq. 136)**

When ``include_zpe = false``:

.. _eq-137:

  Delta H_ads = Delta E_ads + Delta U_vib,thermal - H_gas,

**(Eq. 137)**

and the vibrational free-energy correction is correspondingly evaluated with
the ZPE term removed.  Entropy is unchanged because ZPE is temperature
independent.

The implemented adsorption enthalpy is

.. _eq-138:

  Delta H_ads(T)
  = Delta E_ads
    + Delta ZPE
    + Delta U_vib,thermal(T)
    - H_gas,thermal(T).

**(Eq. 138)**

The adsorption entropy is

.. _eq-139:

  Delta S_ads(T,P)
  = Delta S_vib(T) - S_gas(T,P).

**(Eq. 139)**

The adsorption Gibbs free energy is then

.. _eq-140:

  Delta G_ads(T,P)
  = Delta H_ads(T,P) - T Delta S_ads(T,P).

**(Eq. 140)**

An algebraically equivalent and useful implementation check is

.. _eq-141:

  Delta G_ads
  = Delta E_ads + Delta F_vib - mu_Hg,thermal.

**(Eq. 141)**

The program evaluates both routes and checks that they agree within numerical
tolerance.  This provides an internal thermodynamic-consistency test.

The sign convention is therefore:

* Delta G_ads < 0: adsorption is thermodynamically favorable under the stated
  T and P conditions;
* Delta G_ads > 0: desorption is thermodynamically favored relative to the
  chosen gas reference;
* Delta G_ads = 0: adsorption and desorption are at the thermodynamic
  crossover in this model.

9. Temperature scan
-------------------

The temperature scan does not repeat the electronic-structure calculation and
does not recompute the Hessian.  Once the optimized structures, adsorption
energy, and vibrational frequencies are available, the scan repeatedly
re-evaluates the analytical thermodynamic expressions at each configured
temperature.

For each T the program evaluates and records at least:

* Delta E_ads;
* Delta ZPE;
* Delta U_vib,thermal;
* Delta U_vib,total;
* Delta F_vib;
* S_gas;
* Delta S_vib;
* Delta S_ads;
* H_gas,thermal;
* Delta H_ads; and
* Delta G_ads.

The scan is written to

  hg_au111_thermodynamics_scan.csv

when the output prefix is the supplied default.  In ``mode = scan`` the values
are also printed in a compact terminal table.

For the representative scan in the project documentation, the values are:

+----------+-----------------+---------------------+----------------+
| T (K)    | Delta H_ads     | Delta S_ads         | Delta G_ads    |
+==========+=================+=====================+================+
| 50.00    | -0.49620164 eV  | -2.12731440e-3 eV/K | -0.38983592 eV |
+----------+-----------------+---------------------+----------------+
| 100.00   | -0.49539024 eV  | -2.11729969e-3 eV/K | -0.28366027 eV |
+----------+-----------------+---------------------+----------------+
| 200.00   | -0.49187176 eV  | -2.09333473e-3 eV/K | -0.07320481 eV |
+----------+-----------------+---------------------+----------------+
| 298.15   | -0.48791444 eV  | -2.07726519e-3 eV/K |  0.13142218 eV |
+----------+-----------------+---------------------+----------------+
| 400.00   | -0.48366842 eV  | -2.06502003e-3 eV/K |  0.34233959 eV |
+----------+-----------------+---------------------+----------------+
| 500.00   | -0.47944343 eV  | -2.05559378e-3 eV/K |  0.54835346 eV |
+----------+-----------------+---------------------+----------------+
| 600.00   | -0.47519066 eV  | -2.04784062e-3 eV/K |  0.75351372 eV |
+----------+-----------------+---------------------+----------------+

10. Determination of T_cross
-----------------------------

The code can determine the temperature at which

  Delta G_ads(T) = 0.

The configured ``temperatures`` list defines the search interval.  The code
sorts the positive temperatures, evaluates Delta G on that grid, and searches
for the first adjacent sign change.

If a grid point is numerically zero, it is accepted directly.  Otherwise, if
a sign change is present, a deterministic bisection method is applied.

The bisection uses a temperature tolerance of 1e-4 K and a Gibbs-energy
tolerance of 1e-10 eV, with a maximum of 100 iterations.  Each bisection step
only evaluates the thermodynamic equations; no MACE force calculation and no
phonon calculation is performed.

For the representative calculation,

  T_cross approximately 235.025029 K.

At the reported solution the residual was approximately

  Delta G_ads(T_cross) = -3.18e-8 eV,

with a final numerical bracket of approximately

  [235.02498236, 235.02507596] K.

The root should be interpreted as a result of the stated harmonic,
Gamma-point, slab, gas-reference, pressure, and electronic-structure model.
It is not an experimentally universal desorption temperature.

10.1 Mathematical definition and numerical root criterion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The crossover temperature is the positive solution of

  f(T) = Delta G_ads(T,P) = 0.

The workflow does not solve this equation symbolically.  It evaluates the
thermodynamic function numerically using the already computed electronic
energy and frequency sets.

Given a bracket [T_low, T_high] satisfying

.. _eq-142:

  f(T_low) f(T_high) < 0,

**(Eq. 142)**

the midpoint is

  T_mid = 1/2 (T_low + T_high).

If

  f(T_low) f(T_mid) <= 0,

the new bracket is [T_low, T_mid]; otherwise it is [T_mid, T_high].  The
process is repeated until either the temperature bracket is narrower than the
configured temperature tolerance or the absolute Gibbs-energy residual is
below the configured Gibbs-energy tolerance.

For a bracket containing a single continuous crossing, bisection is guaranteed
to converge.  The reported ``T_cross`` is therefore a numerical root of the
specific model function Delta G_ads(T,P), not an independently optimized
structure at T_cross.

10.2 Analytical pressure shift of the crossover condition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Because

.. _eq-143:

  Delta G_ads(T,P)
    = Delta G_ads(T,P0) - k_B T ln(P/P0),

**(Eq. 143)**

the crossover condition can also be written

.. _eq-144:

  Delta G_ads(T_cross,P0)
    = k_B T_cross ln(P/P0).

**(Eq. 144)**

This relation does not by itself give a closed-form T_cross because the
vibrational terms are temperature dependent, but it provides a useful check on
pressure-dependent crossover calculations.

11. Selectable computational backends
--------------------------------------

11.1 Vibrational backend
~~~~~~~~~~~~~~~~~~~~~~~~

The configuration keyword is

  [phonons]
  method = custom

or

  [phonons]
  method = ase

``custom`` invokes the project's direct finite-difference Hessian
implementation.  ``ase`` invokes ASE's ``Vibrations`` implementation.  Both
operate on the same optimized structure and the same selected mobile
coordinates; they differ in the Hessian/mode implementation.

The terminal output explicitly reports the requested and actual vibrational
backend so that accidental backend substitution is visible.

11.2 Thermodynamics backend
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The configuration keyword is

  [thermodynamics]
  method = custom

or

  [thermodynamics]
  method = ase

The ``custom`` branch evaluates the project's explicit harmonic and ideal-gas
equations described above.

The ``ase`` branch uses ASE's ``HarmonicThermo`` for the vibrational systems and
ASE's ``IdealGasThermo`` for monatomic Hg(g).  The project retains the same
adsorption-energy convention and explicitly handles the ``include_zpe`` switch
so that the two backends are thermodynamically comparable.

The thermodynamics backend is reported once at the beginning of the
thermochemistry section.  It is deliberately not printed inside the repeated
temperature-scan or bisection evaluations.

12. Program architecture
-------------------------

The main execution path is

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

12.1 ``run_hg_au111.py``
~~~~~~~~~~~~~~~~~~~~~~~~~

This is the command-line entry point.  It loads the configuration, applies
command-line overrides, creates the MACE calculator, executes the workflow,
and prints final results.

Typical execution is

  python run_hg_au111.py -c hg_au111_config.ini

12.2 ``config.py``
~~~~~~~~~~~~~~~~~~

``Config`` stores defaults, reads the INI file, converts strings to numerical
values, and validates incompatible settings.

Important controls include:

* slab dimensions and lattice constant;
* geometry-relaxation frozen layers;
* vibrationally frozen layers;
* adsorption site and Hg height;
* optimization thresholds and optimizer;
* phonon backend and finite-difference settings;
* thermodynamics backend;
* temperature and pressure;
* scan temperatures;
* Delta G = 0 search activation; and
* MACE model path.

The validation layer ensures, among other things, that the vibrationally frozen
layer count is not smaller than the geometry-relaxation frozen layer count and
does not exceed the number of Au layers.

12.3 ``surface.py``
~~~~~~~~~~~~~~~~~~~

Builds the Au(111) slab, identifies Au layers, applies bottom-layer geometry
constraints, determines adsorption-site coordinates, adds Hg, and reports
surface details.

12.4 ``optimization.py``
~~~~~~~~~~~~~~~~~~~~~~~~~

Performs constrained structural optimization using ASE optimizers and the MACE
calculator.  The relaxation constraints remain separate from the later
vibrational coordinate selection.

12.5 ``energies.py``
~~~~~~~~~~~~~~~~~~~~

Contains adsorption-energy bookkeeping and isolated-Hg energy evaluation and
prints the electronic-energy summary.

12.6 ``phonons.py``
~~~~~~~~~~~~~~~~~~~

Contains the constrained finite-displacement vibrational machinery.  It:

* identifies the Au layers;
* constructs the vibrational frozen mask;
* retains all atoms in force evaluations;
* identifies mobile vibrational coordinates;
* dispatches to ``custom`` or ``ase``;
* builds/reads the mode information used by the workflow;
* identifies significant imaginary modes;
* applies the positive-frequency cutoff;
* performs Hg-mode projection diagnostics; and
* writes frequency and diagnostic files.

The important conceptual distinction is that ``phonon_atoms`` represent the
coordinate space, not the physical system passed to the calculator.

12.7 ``thermodynamics.py``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Implements the harmonic oscillator formulas, ideal-gas Hg thermodynamics, the
custom thermochemistry backend, and the ASE thermochemistry backend.  It also
contains internal consistency checks between the H-TS and F-mu routes.

12.8 ``workflow.py``
~~~~~~~~~~~~~~~~~~~~

Coordinates the complete calculation.  It builds and optimizes both systems,
calculates the adsorption energy, obtains phonons, evaluates thermochemistry,
performs the scan, optionally searches for T_cross, and writes the final JSON
results.

The temperature scan and T_cross search intentionally reuse the already
calculated vibrational frequencies.

12.9 ``calculator.py``
~~~~~~~~~~~~~~~~~~~~~~

Creates the MACE calculator and provides conservative CUDA-memory cleanup
functions around energy and force evaluations.

12.10 ``constants.py``
~~~~~~~~~~~~~~~~~~~~~~

Contains the physical constants and unit-conversion factors needed for the
vibrational frequency and thermodynamic calculations.

13. Input configuration
-----------------------

A representative configuration is structured as follows:

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
  output_prefix = hg_au111
  verbose = true
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

The pressure is specified in bar in the project configuration.  Internal
conversions to pascal are made where required by ASE or the statistical-
mechanical formulas.

14. Output files and provenance
--------------------------------

The workflow writes files using the configured ``output_prefix``.  Depending
on the selected options and backend, the principal outputs include:

* optimizer trajectory files when ``save_trajectories = true``;
* ``Au_phonons_frequencies_cm-1.dat``;
* ``Au_phonons_signed_frequencies_cm-1.dat``;
* ``Au_phonons_gamma_eigenvalues.dat``;
* ``Au_phonons_mode_diagnostics.dat``;
* corresponding ``HgAu_phonons_*`` files;
* phonon metadata JSON files;
* ``hg_au111_thermodynamics_scan.csv``; and
* ``hg_au111_results.json``.

The frequency file without ``signed`` contains the positive modes retained for
thermochemistry.  The signed-frequency file is intended for stability and
imaginary-mode diagnostics.

The final JSON combines the electronic adsorption energy, thermodynamic
corrections, backend-independent summary quantities, vibrational mode counts,
Hg-dominated frequencies, and the optional T_cross result.


14.5 Corrected constrained-phonon validation
--------------------------------------------

The critical correction in the present workflow is that frozen Au atoms are
retained in every finite-difference force calculation.  They are not removed
from the physical system.  Only their coordinates are excluded from the
vibrational coordinate vector.

For the default 96-atom clean slab with 32 geometry-frozen atoms, force
evaluations still contain all 96 Au atoms, while the mobile-coordinate Hessian
is 192 x 192.

For Hg/Au, force evaluations contain all 97 atoms and the mobile-coordinate
Hessian is 195 x 195.

The expected metadata for the corrected custom finite-difference calculation
contains the conceptual flags

  frozen_atoms_retained_in_force_calculations: true

and

  method: central_finite_difference_full_system_mobile_hessian

These provide a direct provenance check that the old, incorrect
``remove-frozen-atoms`` construction has not been used.

14.6 Historical phonon error and numerical consequences
--------------------------------------------------------

An earlier version of the workflow removed the frozen Au atoms before
calculating displaced forces.  That procedure changed the physical mechanical
response and produced a different vibrational spectrum.

The previously reported values from that incorrect calculation included:

* minimum frequency: about 10.7 cm^-1;
* ZPE correction: about 0.00891 eV;
* Delta G_ads(298 K): about -0.16293 eV;
* Delta G_ads = 0 crossover: about 448.887 K.

These values must not be mixed with results from the corrected constrained
finite-difference construction.

The corrected workflow must be rerun from the optimized structures.  A change
relative to the historical values is expected and is not, by itself, evidence
of a regression: it reflects the correction of the force-evaluation system.

14.7 Backend validation
-----------------------

The code provides two thermochemistry backends:

* ``custom``: explicit project equations for harmonic vibrations and ideal-gas
  Hg;
* ``ase``: ASE ``HarmonicThermo`` plus ``IdealGasThermo``.

The supplied validation record reports that:

* all Python modules compile successfully with ``python -m py_compile``;
* custom backend dispatch was exercised;
* the ASE backend was exercised against a local test double implementing the
  documented ASE thermochemistry equations;
* the resulting thermodynamic quantities agree with the custom backend to
  numerical precision;
* a full ASE/MACE calculation should nevertheless be run in the intended
  ``mace_env`` environment for end-to-end validation.

The implementation also checks internally that

  Delta G_ads = Delta H_ads - T Delta S_ads

and

  Delta G_ads = Delta E_ads + Delta F_vib - mu_Hg,thermal

agree within numerical tolerance.

14.8 Explicit convergence protocol
-----------------------------------

Before treating Delta G or T_cross as publication-level quantitative results,
repeat the phonon calculation with at least

  phonon_displacement = 0.003 Angstrom
  phonon_displacement = 0.005 Angstrom
  phonon_displacement = 0.008 Angstrom

and examine the stability of the frequencies, ZPE, Delta F_vib, Delta G, and
T_cross.

Also test:

* a larger surface cell than 4 x 4;
* increased slab thickness beyond six layers;
* the number of vibrationally frozen layers;
* adsorption sites fcc, hcp, bridge, and top;
* finite-temperature structural effects;
* explicit coverage dependence;
* and, where quantitative periodic vibrational free energies are required,
  Brillouin-zone/q-point sampling.

14.9 Expected mode counts and stability checks
----------------------------------------------

For the default model:

  clean Au: 64 mobile atoms -> 192 modes

  Hg/Au: 65 mobile atoms -> 195 modes.

The positive-frequency thermochemistry list is obtained only after the signed
spectrum has been examined.

With

  fail_on_imaginary = true,

a significant imaginary mode terminates the thermochemical stage rather than
being silently discarded.  A low-frequency or imaginary result should instead
trigger inspection of the optimized geometry, displacement amplitude, slab
model, frozen-layer definition, and convergence.

Because the bottom Au layers are fixed, the constrained system does not have
the same translational zero modes as a freely translating isolated system.
Accordingly, the acoustic sum-rule correction is disabled.

15. Computational scope and limitations
----------------------------------------

15.1 Gamma-point approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The vibrational calculation is a Gamma-point calculation for the finite
surface supercell.  It does not integrate phonon free energies over the
Brillouin zone.  For a periodic surface, quantitative vibrational free-energy
convergence may require explicit q-point sampling or an appropriate periodic
phonon treatment.

15.2 Harmonic approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The thermochemistry is harmonic.  Anharmonicity, temperature-dependent
structural relaxation, phonon-phonon interactions, and thermal expansion are
not included.

15.3 Frozen-layer approximation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The vibrational free energy corresponds to the selected constrained coordinate
space.  Changing ``n_frozen_layers_vibrations`` changes the vibrational model
without changing the optimized geometry.  This parameter therefore represents
a deliberate physical/modeling choice and should be tested for convergence.

15.4 Surface and slab convergence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The 4x4x6 model should not automatically be regarded as converged.  Before
using quantitative Delta G or T_cross values for publication, test at least
surface-cell size, slab thickness, number of frozen layers, finite-difference
displacement, and the treatment of low-frequency modes.

15.5 Adsorption-site comparison
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The supplied calculation uses an fcc starting site.  Assigning fcc as the
thermodynamic equilibrium site requires comparison with competing sites such
as top, bridge, and hcp, with consistent geometry and vibrational treatment.

15.6 Gas-reference pressure
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The calculated Delta G_ads and T_cross are pressure dependent because the gas
entropy and chemical potential depend on Hg pressure.  Report the pressure
explicitly with any thermodynamic result.


15.7 Comparison of vibrational constraint models
------------------------------------------------

Two vibrational treatments were documented for the same electronic adsorption
energy.  The purpose of this comparison is to show the effect of the selected
vibrational coordinate space, not to change the relaxed geometry.

Standard treatment:

  Delta ZPE             = 0.00881388 eV
  Delta H_ads           = -0.48791444 eV
  Delta S_ads           = -2.07727 x 10^-3 eV/K
  Delta G_ads(298 K)    =  0.13142 eV
  Delta G_ads = 0 T     = about 235 K

Frozen-slab vibrational treatment:

  Delta ZPE             = 0.00637762 eV
  Delta H_ads           = -0.48825772 eV
  Delta S_ads           = -2.05303 x 10^-3 eV/K
  Delta G_ads(298 K)    =  0.12385 eV
  Delta G_ads = 0 T     = about 238 K

The electronic adsorption energy is unchanged between these treatments.  The
difference comes from the modified vibrational spectrum and therefore from
Delta F_vib, Delta U_vib, and Delta S_vib.

16. Representative numerical result
------------------------------------

For the documented custom-backend calculation with two geometry-frozen Au
layers and two vibrationally frozen Au layers, the representative values at
298.15 K and 1e-5 bar were approximately:

  Delta E_ads       = -0.50132288 eV
  Delta ZPE        =  0.00881388 eV
  Delta U_vib      =  0.06882601 eV
  H_gas,thermal    =  0.06423145 eV
  Delta H_ads      = -0.48791444 eV
  Delta S_ads      = -2.077265e-3 eV/K
  Delta F_vib      = -0.13949947 eV
  mu_Hg,thermal    = -0.77224453 eV
  Delta G_ads      =  0.13142218 eV
  T_cross          =  235.025029 K

These numbers are model-specific and should be reproduced only with the same
structure, calculator/model, vibrational coordinate definition, frequency
cutoff, pressure, and thermochemistry conventions.


16.1 Explicit representative temperature scan
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For the representative custom-backend calculation at 1e-5 bar, the documented
scan is:

  T (K)     Delta H_ads (eV)     Delta S_ads (eV/K)     Delta G_ads (eV)
  50.00     -0.49620164          -2.12731440e-3         -0.38983592
  100.00    -0.49539024          -2.11729969e-3         -0.28366027
  200.00    -0.49187176          -2.09333473e-3         -0.07320481
  298.15    -0.48791444          -2.07726519e-3          0.13142218
  400.00    -0.48366842          -2.06502003e-3          0.34233959
  500.00    -0.47944343          -2.05559378e-3          0.54835346
  600.00    -0.47519066          -2.04784062e-3          0.75351372

The corresponding thermodynamic scan is written to

  hg_au111_thermodynamics_scan.csv

Publication-oriented plots documented by the supplied supplementary material
are:

  DeltaH_ads_vs_T.png
  DeltaS_ads_vs_T.png
  DeltaG_ads_vs_T.png

These figures visualize the temperature dependence of adsorption enthalpy,
entropy, and Gibbs free energy.

16.2 Physical interpretation of the temperature dependence
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The electronic interaction is favorable because

  Delta E_ads < 0.

Nevertheless, adsorption removes the large translational entropy associated
with a gas-phase Hg atom.  Thus

  Delta S_ads < 0,

and the term

  -T Delta S_ads

becomes increasingly positive as temperature rises.

The calculated Delta G_ads therefore becomes less favorable with increasing
temperature in the documented pressure range.  The crossover temperature is
the point where the favorable electronic and vibrational contributions are
balanced by the gas-phase thermodynamic penalty.

17. Recommended reproducibility protocol
------------------------------------------

For a publication-quality calculation:

1. Record the exact MACE model file and software versions.
2. Record the complete INI configuration.
3. Save the optimized Au and Hg/Au structures.
4. Save the signed and positive vibrational frequencies and metadata.
5. Verify that significant imaginary modes are absent.
6. Check convergence with respect to finite-difference displacement.
7. Check convergence with respect to slab thickness and surface cell.
8. Check convergence with respect to the number of vibrationally frozen Au
   layers.
9. Repeat with both ``custom`` and ``ase`` backends as an implementation
   cross-check.
10. Report temperature, Hg pressure, gas reference pressure, and ZPE
    convention together with Delta G and T_cross.

The resulting Delta G(T) and T_cross should be regarded as properties of the
specified computational model, not as standalone experimental observables.


17.3 Results from the supplied ASE/MACE production logs
-------------------------------------------------------

The uploaded production logs provide a consistent set of electronic energies
and ASE vibrational/thermochemical results for several Hg pressures.  The
electronic energies are pressure-independent:

  E_Au     = -304.8265897028 eV
  E_HgAu   = -305.4524593272 eV
  E_Hg     = -0.1245467400 eV

and therefore

  Delta E_ads = -0.5013228844 eV.

For the standard vibrational treatment (all 192 clean-Au mobile modes and all
195 Hg/Au mobile modes retained), the spectra reported in the logs contain no
significant imaginary modes.  The clean Au frequencies span 4.309--140.294
cm^-1 and the Hg/Au frequencies span 4.281--140.090 cm^-1.  Two Hg-dominated
Hg/Au modes are identified at approximately 19.312 and 19.517 cm^-1 using the
reported projection threshold of 0.200.

At 298.15 K and 1e-5 bar, the standard treatment gives:

  Delta ZPE          =  0.00881388 eV
  Delta U_vib       =  0.06882601 eV
  H_Hg,thermal      =  0.06423145 eV
  Delta H_ads       = -0.48791444 eV
  S_Hg              =  2.80555417e-3 eV/K
  Delta S_vib       =  7.28288986e-4 eV/K
  Delta S_ads       = -2.07726519e-3 eV/K
  Delta F_vib       = -0.13949947 eV
  mu_Hg,thermal     = -0.77224453 eV
  Delta G_ads       =  0.13142218 eV.

The corresponding zero of the Gibbs free energy is

  T_cross = 235.025029 K,

with the bisection residual reported as

  Delta G_ads(T_cross) = -3.16e-8 eV.

The same pressure-independent electronic adsorption energy is combined with
the pressure-dependent gas chemical potential.  The standard-vibrational ASE
logs give the following pressure dependence at 298.15 K:

  P(Hg) [bar]       Delta G_ads [eV]       T_cross [K]
  ------------------------------------------------------
  1                  -0.16437430            450.921106
  1e-1               -0.10521497            380.457339
  1e-2               -0.04605564            329.270472
  1e-3                0.01310369            290.354684
  1e-4                0.07226302            259.745055
  1e-5                0.13142218            235.025029

This pressure series demonstrates directly the expected logarithmic ideal-gas
pressure dependence.  At fixed temperature, decreasing Hg pressure makes
Delta G_ads less favorable.  Over this pressure interval, the crossover
temperature correspondingly decreases from about 451 K at 1 bar to 235 K at
1e-5 bar.

17.4 Standard-vibrational temperature scan at multiple pressures
-----------------------------------------------------------------

For the standard vibrational treatment, Delta H_ads is essentially
pressure-independent because the pressure enters through the gas entropy and
chemical potential rather than through the electronic or solid vibrational
terms.  The supplied logs give:

  P = 1 bar
  T(K)       Delta H_ads(eV)   Delta S_ads(eV/K)   Delta G_ads(eV)
  50         -0.49620164       -1.13520814e-3       -0.43944123
  100        -0.49539024       -1.12519345e-3       -0.38287090
  200        -0.49187176       -1.10122850e-3       -0.27162606
  298.15     -0.48791444       -1.08515897e-3       -0.16437430
  400        -0.48366843       -1.07291381e-3       -0.05450290
  500        -0.47944344       -1.06348757e-3        0.05230035
  600        -0.47519067       -1.05573441e-3        0.15824998

  P = 1e-1 bar
  T(K)       Delta H_ads(eV)   Delta S_ads(eV/K)   Delta G_ads(eV)
  50         -0.49620164       -1.33362951e-3       -0.42952016
  100        -0.49539024       -1.32361482e-3       -0.36302876
  200        -0.49187176       -1.29964987e-3       -0.23194179
  298.15     -0.48791444       -1.28358033e-3       -0.10521497
  400        -0.48366843       -1.27133518e-3        0.02486564
  500        -0.47944344       -1.26190893e-3        0.15151103
  600        -0.47519067       -1.25415577e-3        0.27730280

  P = 1e-2 bar
  T(K)       Delta H_ads(eV)   Delta S_ads(eV/K)   Delta G_ads(eV)
  50         -0.49620164       -1.53205087e-3       -0.41959909
  100        -0.49539024       -1.52203618e-3       -0.34318662
  200        -0.49187176       -1.49807123e-3       -0.19225751
  298.15     -0.48791444       -1.48200169e-3       -0.04605564
  400        -0.48366843       -1.46975654e-3        0.10423419
  500        -0.47944344       -1.46033030e-3        0.25072171
  600        -0.47519067       -1.45257714e-3        0.39635562

  P = 1e-3 bar
  T(K)       Delta H_ads(eV)   Delta S_ads(eV/K)   Delta G_ads(eV)
  50         -0.49620164       -1.73047223e-3       -0.40967802
  100        -0.49539024       -1.72045755e-3       -0.32334449
  200        -0.49187176       -1.69649259e-3       -0.15257324
  298.15     -0.48791444       -1.68042306e-3        0.01310369
  400        -0.48366843       -1.66817790e-3        0.18360273
  500        -0.47944344       -1.65875166e-3        0.34993239
  600        -0.47519067       -1.65099850e-3        0.51540844

  P = 1e-4 bar
  T(K)       Delta H_ads(eV)   Delta S_ads(eV/K)   Delta G_ads(eV)
  50         -0.49620164       -1.92889360e-3       -0.39975696
  100        -0.49539024       -1.91887891e-3       -0.30350235
  200        -0.49187176       -1.89491396e-3       -0.11288897
  298.15     -0.48791444       -1.87884442e-3        0.07226302
  400        -0.48366843       -1.86659927e-3        0.26297128
  500        -0.47944344       -1.85717302e-3        0.44914307
  600        -0.47519067       -1.84941986e-3        0.63446125

17.5 Frozen-slab vibrational control from the supplied logs
------------------------------------------------------------

A separate control calculation freezes all Au atoms for the vibrational
calculation while retaining Hg as mobile.  The clean Au system therefore has
zero vibrational modes in this control, whereas Hg/Au retains only the three
Hg Cartesian vibrational modes.  The Hg/Au control frequencies are

  22.2127, 22.3355, and 58.3297 cm^-1,

with no significant imaginary modes.

The resulting vibrational ZPE difference is

  Delta ZPE = 0.00637762 eV.

At 298.15 K and 1e-5 bar, the frozen-slab control gives:

  Delta U_vib       =  0.07091897 eV
  H_Hg,thermal      =  0.06423143 eV
  Delta H_ads       = -0.48825772 eV
  Delta S_vib       =  7.52526635e-4 eV/K
  Delta S_ads       = -2.05302780e-3 eV/K
  Delta F_vib       = -0.14706923 eV
  mu_Hg,thermal     = -0.77224463 eV
  Delta G_ads       =  0.12385252 eV.

The corresponding crossover is

  T_cross = 237.957615 K,

with a reported residual of approximately 1.41e-8 eV.

At 1 bar the same frozen-slab vibrational control gives

  Delta G_ads(298.15 K) = -0.17194413 eV
  T_cross               = 461.759520 K.

The frozen-slab control therefore shifts the 1e-5-bar crossover from
235.025 K to 237.958 K, while shifting the 1-bar crossover from 450.921 K to
461.760 K.  This control changes only the vibrational coordinate space; the
electronic adsorption energy remains

  Delta E_ads = -0.5013228844 eV.

17.6 Frozen-slab temperature scan
---------------------------------

At 1e-5 bar, the frozen-slab vibrational control gives:

  T(K)       Delta H_ads(eV)   Delta S_ads(eV/K)   Delta G_ads(eV)
  50         -0.49790851       -2.11795637e-3       -0.39201069
  100        -0.49636703       -2.09724285e-3       -0.28664274
  200        -0.49237968       -2.06978456e-3       -0.07842277
  298.15     -0.48825772       -2.05302780e-3        0.12385252
  400        -0.48392500       -2.04052892e-3        0.33228657
  500        -0.47964895       -2.03098780e-3        0.53584495
  600        -0.47536204       -2.02317206e-3        0.73854119

At 1 bar, the corresponding frozen-slab scan is:

  T(K)       Delta H_ads(eV)   Delta S_ads(eV/K)   Delta G_ads(eV)
  50         -0.49790851       -1.12584956e-3       -0.44161603
  100        -0.49636703       -1.10513604e-3       -0.38585342
  200        -0.49237968       -1.07767775e-3       -0.27684413
  298.15     -0.48825772       -1.06092098e-3       -0.17194413
  400        -0.48392500       -1.04842210e-3       -0.06455616
  500        -0.47964895       -1.03888098e-3        0.03979154
  600        -0.47536204       -1.03106524e-3        0.14327710

The associated ASE vibrational thermochemistry output for the frozen-slab
control reports, at 298.15 K,

  ZPE = 0.0063776208 eV
  F_vib = -0.1470692256 eV
  U_vib = 0.0772965905 eV
  S_vib = 7.52526635e-4 eV/K.

These quantities refer to the Hg/Au control spectrum because the clean Au
vibrational coordinate space is empty when every Au atom is frozen.

17.7 Backend and log provenance
--------------------------------

The supplied production logs explicitly identify

  Vibrational backend: ASE Vibrations

and, for the thermochemical runs,

  Thermodynamics backend: ASE HarmonicThermo + IdealGasThermo.

Thus the numerical pressure series and the detailed values in sections 17.3--
17.6 are ASE-backend production results.  They should be distinguished from
the project-custom equations/backend discussed elsewhere in this document.

The agreement of the electronic energies across all uploaded runs provides a
useful reproducibility check.  Likewise, the identical standard-vibrational
spectra across the pressure series show that changing Hg pressure affects the
gas thermodynamic reference but does not alter the calculated slab Hessians
in these runs.


19. Current package status
--------------------------

The source package described by this document contains the modular geometry,
energy, phonon, thermochemistry, backend-dispatch, scan, and Delta G=0
functionality described above.

The exact uploaded package used as the basis for this document does not contain
a separate ``[restart]`` section for reading previously saved relaxed
structures and vibrational data.  Therefore this document does not claim that
file-based restart/post-processing is part of that exact package.  If a restart
extension is added, its saved-structure and saved-vibrational-data conventions
should be documented separately and validated against the metadata of the
calculation that produced those files.

End of document
---------------


20. Consolidated-source completeness statement
----------------------------------------------

This document consolidates the applied theory, computational methodology,
validation information, numerical benchmarks, and reproducibility instructions
contained in the supplied project documentation.

In particular, it preserves:

* the corrected full-system constrained finite-difference phonon construction;
* the Gamma-point-only scientific scope;
* the significant-imaginary-mode stopping policy;
* the disabled acoustic sum-rule correction for the constrained slab;
* the exact default slab/model/pressure/temperature parameters;
* the complete harmonic vibrational statistical mechanics;
* the complete monatomic ideal-gas Hg statistical mechanics;
* the pressure and reference-pressure dependence;
* the H-TS and F-mu free-energy formulations and their consistency check;
* the explicit temperature-scan values;
* the two documented vibrational-constraint comparisons;
* the historical incorrect-phonon numerical results and their cause;
* the backend validation record;
* the exact recommended finite-difference convergence points;
* the expected mode and force-evaluation counts;
* the adsorption-site validation recommendation;
* the publication-oriented output figure names; and
* the reproducibility and publication-scope limitations.

The document intentionally distinguishes representative numerical results from
converged physical predictions.  Gamma-point harmonic thermochemistry,
finite-size slabs, frozen-layer constraints, and the chosen MACE model define
the computational model whose Delta G and T_cross are reported.