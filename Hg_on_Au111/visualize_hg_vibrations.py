#!/usr/bin/env python3
"""
Visualize Hg-dominated Gamma-point vibrations from an existing ASE Vibrations cache.

This script does NOT run new force calculations.

It:
1. Reads the optimized full Hg/Au(111) structure from an ASE trajectory.
2. Derives the mobile-atom coordinate space directly from the existing ASE cache.
3. Reopens the existing ASE Vibrations cache from --phonon-dir.
4. Reproduces the production dynamical-matrix construction and mode analysis.
5. Identifies Hg-dominated positive-frequency modes.
6. Converts mass-weighted eigenvectors to Cartesian displacement vectors.
7. Embeds the displacements into the full structure, keeping frozen atoms fixed.
8. Writes ASE .traj and XYZ animation files suitable for `ase gui`.

The production phonons.py uses:
    Vibrations(full_atoms, indices=mobile_indices, delta=0.005, nfree=2,
               name=prefix)
followed by:
    vib.run()
    vib.read(method="standard", direction="central")
The cache is indexed by the original/full atom indices.

The same force-constant -> mass-weighting -> eigensolver procedure is
reproduced here so that the visualization corresponds to the production
phonon frequencies and eigenvectors.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from ase import Atoms
from ase.calculators.singlepoint import SinglePointCalculator
from ase.io import read, write
from ase.vibrations import Vibrations


# Conversion used by ASE for sqrt(eV / Angstrom^2 / amu) -> cm^-1.
ASE_TO_CM = 521.47083


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Visualize Hg-dominated Gamma-point modes from an ASE Vibrations cache."
    )

    parser.add_argument(
        "--structure",
        default="HgAu_opt.traj",
        help="Optimized full Hg/Au structure trajectory (default: HgAu_opt.traj).",
    )
    parser.add_argument(
        "--phonon-dir",
        default="HgAu_phonons",
        help="ASE Vibrations cache directory/name (default: HgAu_phonons).",
    )
    parser.add_argument(
        "--output-dir",
        default="Hg_vibrations",
        help="Directory for generated trajectories and diagnostics.",
    )
    parser.add_argument(
        "--frozen",
        type=int,
        default=None,
        help=(
            "Number of frozen atoms at the beginning of the full structure. "
            "By default this is derived from the ASE cache filenames."
        ),
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=0.20,
        help="Minimum Hg projection for automatic mode selection (default: 0.20).",
    )
    parser.add_argument(
        "--cutoff",
        type=float,
        default=1.0,
        help="Positive-frequency cutoff in cm^-1 (default: 1.0).",
    )
    parser.add_argument(
        "--amplitude",
        type=float,
        default=0.35,
        help="Maximum Cartesian displacement in Angstrom for visualization (default: 0.35).",
    )
    parser.add_argument(
        "--frames",
        type=int,
        default=21,
        help="Number of frames per animation (default: 21).",
    )
    parser.add_argument(
        "--phase-cycles",
        type=float,
        default=1.0,
        help="Number of complete oscillation cycles in each trajectory (default: 1.0).",
    )
    parser.add_argument(
        "--n-modes",
        type=int,
        default=3,
        help="Number of Hg-dominated modes to write automatically (default: 3).",
    )
    parser.add_argument(
        "--all-hg",
        action="store_true",
        help="Write all stable modes with Hg projection >= threshold.",
    )
    parser.add_argument(
        "--mode",
        type=int,
        action="append",
        default=None,
        help=(
            "Explicit zero-based mode index to visualize. Repeat the option for "
            "multiple modes, e.g. --mode 0 --mode 1."
        ),
    )

    return parser.parse_args()


def load_structure(path: str) -> Atoms:
    structure = read(path, index=-1)

    if not isinstance(structure, Atoms):
        raise TypeError(f"Could not read an ASE Atoms object from {path!r}.")

    if len(structure) == 0:
        raise ValueError(f"Structure {path!r} contains no atoms.")

    return structure


def discover_cache_indices(phonon_dir: str) -> tuple[np.ndarray, set[int]]:
    """Derive the vibrational mobile indices directly from an ASE cache.

    ASE Vibrations caches displacement keys using the *original/full-system*
    atom indices, e.g. ``64x+``, ``64x-``, ..., ``96z-``.  This function
    reconstructs those original indices instead of assuming that the first
    mobile atom is index zero.

    Every mobile atom must have all six central-difference displacement files
    (x+/x-/y+/y-/z+/z-).  ``cache.eq.json`` is the equilibrium entry and is
    not itself used to infer the mobile atom list.
    """
    cache_dir = Path(phonon_dir)
    if not cache_dir.is_dir():
        raise FileNotFoundError(
            f"Phonon cache directory not found: {cache_dir}"
        )

    import re

    pattern = re.compile(r"^cache\.(\d+)([xyz])([+-])\.json$")
    required = {f"{axis}{sign}" for axis in "xyz" for sign in "+-"}
    found: dict[int, set[str]] = {}

    for path in cache_dir.iterdir():
        match = pattern.match(path.name)
        if match is None:
            continue
        atom_index = int(match.group(1))
        component = match.group(2) + match.group(3)
        found.setdefault(atom_index, set()).add(component)

    if not found:
        raise RuntimeError(
            f"No ASE Vibrations displacement files matching "
            f"cache.<atom><x|y|z><+|->.json were found in {cache_dir}."
        )

    incomplete = {
        index: sorted(required - components)
        for index, components in found.items()
        if components != required
    }
    if incomplete:
        details = "; ".join(
            f"{index}: missing {','.join(missing)}"
            for index, missing in sorted(incomplete.items())
        )
        raise RuntimeError(
            "The ASE cache is incomplete; cannot safely infer the mobile "
            f"atom set. {details}"
        )

    mobile_indices = np.array(sorted(found), dtype=int)
    return mobile_indices, set(int(i) for i in mobile_indices)


def make_mobile_atoms(
    full: Atoms, mobile_indices: np.ndarray
) -> tuple[Atoms, np.ndarray, np.ndarray]:
    """Construct the mobile coordinate system from original atom indices."""
    mobile_indices = np.asarray(mobile_indices, dtype=int)
    if mobile_indices.ndim != 1 or mobile_indices.size == 0:
        raise ValueError("The cache-derived mobile index list is empty or invalid.")
    if np.any(mobile_indices < 0) or np.any(mobile_indices >= len(full)):
        raise ValueError(
            "The phonon cache contains atom indices outside the supplied "
            f"structure (0..{len(full)-1}): {mobile_indices.tolist()}"
        )

    mobile = full[mobile_indices].copy()
    mobile.calc = SinglePointCalculator(
        mobile,
        energy=0.0,
        forces=np.zeros((len(mobile), 3), dtype=float),
    )

    frozen_mask = np.ones(len(full), dtype=bool)
    frozen_mask[mobile_indices] = False
    frozen_indices = np.where(frozen_mask)[0]
    return mobile, mobile_indices, frozen_indices


def load_phonon_cache(
    full: Atoms,
    mobile_indices: np.ndarray,
    phonon_dir: str,
    displacement: float = 0.005,
    nfree: int = 2,
) -> Vibrations:
    """Reopen the production ASE Vibrations cache without recalculating forces."""
    vib = Vibrations(
        full,
        indices=mobile_indices.tolist(),
        delta=float(displacement),
        nfree=int(nfree),
        name=phonon_dir,
    )
    vib.read(method="standard", direction="central")
    return vib

def gamma_dynamical_matrix(vib: Vibrations) -> np.ndarray:
    """Reproduce the production mass-weighted Gamma dynamical matrix exactly."""
    H = np.asarray(vib.H, dtype=float)
    n = len(vib.indices)
    expected = 3 * n
    if H.shape != (expected, expected):
        raise ValueError(
            f"Vibrational Hessian has shape {H.shape}, "
            f"expected {(expected, expected)}"
        )

    masses = np.asarray(vib.atoms.get_masses(), dtype=float)[vib.indices]
    if np.any(masses <= 0.0):
        raise ValueError("All vibrational atoms must have positive masses.")

    m_inv_sqrt = np.repeat(masses ** -0.5, 3)
    D = H * np.outer(m_inv_sqrt, m_inv_sqrt)
    return 0.5 * (D + D.T)

def gamma_modes(D: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Diagonalize the mass-weighted Gamma dynamical matrix.

    This reproduces the production phonons.py convention:
        sign(lambda) * sqrt(abs(lambda)) * ASE_TO_CM
    for signed frequencies in cm^-1.
    """
    eigenvalues, eigenvectors = np.linalg.eigh(D)
    signed_cm = (
        np.sign(eigenvalues)
        * np.sqrt(np.abs(eigenvalues))
        * ASE_TO_CM
    )
    return signed_cm, eigenvectors


def hg_projections(
    mobile: Atoms,
    eigenvectors: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Return Hg atom indices in the mobile system and Hg projection per mode.

    Eigenvectors here are mass-weighted normal-mode eigenvectors. The projection
    is the squared eigenvector amplitude on the Hg Cartesian degrees of freedom,
    matching the diagnostic used by production phonons.py.
    """
    hg_indices = np.array(
        [i for i, atom in enumerate(mobile) if atom.symbol == "Hg"],
        dtype=int,
    )

    if len(hg_indices) == 0:
        raise ValueError("No Hg atom was found in the mobile phonon system.")

    hg_dofs = np.concatenate(
        [3 * hg_indices + axis for axis in range(3)]
    )

    # eigenvectors has shape (3N, 3N), columns are normal modes.
    projection = np.sum(eigenvectors[hg_dofs, :] ** 2, axis=0)

    return hg_indices, projection


def cartesian_mode_vector(
    mobile: Atoms,
    mass_weighted_vector: np.ndarray,
    amplitude: float,
) -> np.ndarray:
    """
    Convert a mass-weighted eigenvector to a Cartesian displacement vector.

    The resulting vector is normalized so that its largest per-atom displacement
    has magnitude `amplitude` Angstrom.
    """
    n = len(mobile)

    masses = np.asarray(mobile.get_masses(), dtype=float)
    sqrt_masses = np.sqrt(masses)

    q = mass_weighted_vector.reshape(n, 3)

    # q is mass-weighted. Convert to Cartesian displacement.
    displacement = q / sqrt_masses[:, None]

    per_atom = np.linalg.norm(displacement, axis=1)
    max_atom = float(np.max(per_atom))

    if not np.isfinite(max_atom) or max_atom <= 0.0:
        raise ValueError("Mode has zero or invalid Cartesian displacement.")

    displacement *= amplitude / max_atom

    return displacement


def build_animation(
    full: Atoms,
    mobile_indices: np.ndarray,
    displacement: np.ndarray,
    frames: int,
    phase_cycles: float,
) -> list[Atoms]:
    """
    Embed a sinusoidal Cartesian mode into the full structure.

    Frozen atoms are left exactly at their equilibrium coordinates.
    """
    if frames < 2:
        raise ValueError("--frames must be at least 2.")

    if phase_cycles <= 0.0:
        raise ValueError("--phase-cycles must be positive.")

    result: list[Atoms] = []

    phases = np.linspace(
        0.0,
        2.0 * np.pi * phase_cycles,
        frames,
        endpoint=True,
    )

    for phase in phases:
        frame = full.copy()
        frame.positions[mobile_indices] += displacement * np.sin(phase)
        result.append(frame)

    return result


def safe_frequency_tag(frequency_cm: float) -> str:
    """
    Create a filename-safe frequency string with three decimal places.
    """
    return f"{frequency_cm:07.3f}"


def write_mode(
    full: Atoms,
    mobile_indices: np.ndarray,
    displacement: np.ndarray,
    mode_index: int,
    frequency_cm: float,
    projection: float,
    output_dir: Path,
    frames: int,
    phase_cycles: float,
) -> dict:
    animation = build_animation(
        full=full,
        mobile_indices=mobile_indices,
        displacement=displacement,
        frames=frames,
        phase_cycles=phase_cycles,
    )

    freq_tag = safe_frequency_tag(frequency_cm)

    traj_path = output_dir / f"Hg_mode_{freq_tag}_cm-1.traj"
    xyz_path = output_dir / f"Hg_mode_{freq_tag}_cm-1.xyz"

    write(traj_path, animation, format="traj")
    write(xyz_path, animation, format="xyz")

    return {
        "mode_index": int(mode_index),
        "frequency_cm-1": float(frequency_cm),
        "hg_projection": float(projection),
        "trajectory": str(traj_path),
        "xyz": str(xyz_path),
    }


def main() -> None:
    args = parse_args()

    if args.amplitude <= 0.0:
        raise ValueError("--amplitude must be positive.")

    if args.frames < 2:
        raise ValueError("--frames must be at least 2.")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    full = load_structure(args.structure)

    # The production cache is authoritative for the vibrational coordinate
    # space.  In particular, do not infer the mobile set from a hard-coded
    # frozen-atom count: production ASE Vibrations stores keys using the
    # original/full atom indices.
    cache_mobile_indices, _ = discover_cache_indices(args.phonon_dir)

    if args.frozen is not None:
        # Retain the CLI option as an explicit validation aid, but never use
        # it to construct the coordinate space when the cache is available.
        expected_mobile = np.arange(args.frozen, len(full), dtype=int)
        if not np.array_equal(expected_mobile, cache_mobile_indices):
            raise ValueError(
                f"--frozen={args.frozen} does not match the mobile indices "
                "encoded by the production cache. "
                f"Cache mobile indices: {cache_mobile_indices.tolist()}"
            )

    mobile, mobile_indices, frozen_indices = make_mobile_atoms(
        full,
        mobile_indices=cache_mobile_indices,
    )

    phon = load_phonon_cache(
        full,
        mobile_indices,
        args.phonon_dir,
        displacement=0.005,
        nfree=2,
    )

    # Verify ASE accepted exactly the cache-derived original indices.
    vib_indices = np.asarray(phon.indices, dtype=int)
    if not np.array_equal(vib_indices, mobile_indices):
        raise RuntimeError(
            "ASE Vibrations returned indices different from those inferred "
            "from the production cache: "
            f"{vib_indices.tolist()} vs {mobile_indices.tolist()}"
        )

    D = gamma_dynamical_matrix(phon)
    frequencies_cm, eigenvectors = gamma_modes(D)
    hg_mobile_indices, hg_projection = hg_projections(
        mobile,
        eigenvectors,
    )

    positive = frequencies_cm >= args.cutoff
    selected_indices = np.where(positive)[0]

    if len(selected_indices) == 0:
        raise RuntimeError(
            f"No positive phonon modes remain at the cutoff "
            f"{args.cutoff:.3f} cm^-1."
        )

    # Significant imaginary modes are diagnostic errors for this stable
    # visualization workflow. Tiny negative/near-zero numerical modes are
    # allowed below the positive cutoff.
    significant_imaginary = frequencies_cm < -args.cutoff

    if np.any(significant_imaginary):
        bad = np.where(significant_imaginary)[0]
        preview = ", ".join(f"{frequencies_cm[i]:.3f}" for i in bad[:10])
        raise RuntimeError(
            "Significant imaginary modes were found: "
            f"{preview} cm^-1. Check the optimized structure and phonon cache."
        )

    if args.mode is not None:
        mode_indices = list(dict.fromkeys(args.mode))

        for mode_index in mode_indices:
            if mode_index < 0 or mode_index >= len(frequencies_cm):
                raise ValueError(
                    f"Mode index {mode_index} is outside the valid range "
                    f"0..{len(frequencies_cm)-1}."
                )
            if frequencies_cm[mode_index] < args.cutoff:
                raise ValueError(
                    f"Mode {mode_index} has frequency "
                    f"{frequencies_cm[mode_index]:.6f} cm^-1, below the "
                    f"{args.cutoff:.3f} cm^-1 cutoff."
                )

    else:
        hg_candidates = [
            int(i)
            for i in selected_indices
            if hg_projection[i] >= args.threshold
        ]

        # Highest Hg projection first; frequency is used as a deterministic
        # tie-breaker.
        hg_candidates.sort(
            key=lambda i: (-float(hg_projection[i]), float(frequencies_cm[i]))
        )

        if args.all_hg:
            mode_indices = hg_candidates
        else:
            mode_indices = hg_candidates[: args.n_modes]

    if not mode_indices:
        raise RuntimeError(
            "No modes satisfy the requested selection. "
            f"Try lowering --threshold (currently {args.threshold}) "
            "or use explicit --mode."
        )

    equilibrium_path = output_dir / "HgAu_equilibrium.traj"
    write(equilibrium_path, full, format="traj")

    diagnostics = {
        "structure": str(Path(args.structure).resolve()),
        "phonon_dir": str(Path(args.phonon_dir).resolve()),
        "n_full_atoms": int(len(full)),
        "n_mobile_atoms": int(len(mobile)),
        "frozen_atoms": int(len(frozen_indices)),
        "frozen_indices": frozen_indices.tolist(),
        "mobile_indices": mobile_indices.tolist(),
        "index_source": "ASE Vibrations cache filenames",
        "hg_mobile_indices": hg_mobile_indices.tolist(),
        "cutoff_cm-1": float(args.cutoff),
        "hg_projection_threshold": float(args.threshold),
        "visualization_amplitude_A": float(args.amplitude),
        "frames": int(args.frames),
        "phase_cycles": float(args.phase_cycles),
        "n_positive_modes": int(len(selected_indices)),
        "frequency_range_positive_cm-1": [
            float(np.min(frequencies_cm[selected_indices])),
            float(np.max(frequencies_cm[selected_indices])),
        ],
        "selected_modes": [],
    }

    print()
    print("Hg vibration visualization")
    print("===========================")
    print(f"Structure       : {args.structure}")
    print(f"Phonon cache    : {args.phonon_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Full atoms      : {len(full)}")
    print(f"Frozen atoms    : {len(frozen_indices)}")
    print(f"Frozen indices  : {frozen_indices.tolist()}")
    print(f"Mobile indices  : {mobile_indices.tolist()}")
    print(f"Mobile atoms    : {len(mobile)}")
    print(f"Hg mobile index : {hg_mobile_indices.tolist()}")
    print(f"Positive modes  : {len(selected_indices)}")
    print(
        f"Positive range  : "
        f"{np.min(frequencies_cm[selected_indices]):.3f}–"
        f"{np.max(frequencies_cm[selected_indices]):.3f} cm^-1"
    )
    print()

    print("Selected modes:")
    for mode_index in mode_indices:
        print(
            f"  mode {mode_index:3d}: "
            f"{frequencies_cm[mode_index]:10.3f} cm^-1, "
            f"Hg projection = {hg_projection[mode_index]:.4f}"
        )

    print()

    for mode_index in mode_indices:
        displacement = cartesian_mode_vector(
            mobile=mobile,
            mass_weighted_vector=eigenvectors[:, mode_index],
            amplitude=args.amplitude,
        )

        info = write_mode(
            full=full,
            mobile_indices=mobile_indices,
            displacement=displacement,
            mode_index=mode_index,
            frequency_cm=frequencies_cm[mode_index],
            projection=hg_projection[mode_index],
            output_dir=output_dir,
            frames=args.frames,
            phase_cycles=args.phase_cycles,
        )

        diagnostics["selected_modes"].append(info)

        print(f"Wrote: {info['trajectory']}")
        print(f"       {info['xyz']}")

    diagnostics_path = output_dir / "visualization_diagnostics.json"
    diagnostics_path.write_text(
        json.dumps(diagnostics, indent=2),
        encoding="utf-8",
    )

    print()
    print(f"Equilibrium: {equilibrium_path}")
    print(f"Diagnostics: {diagnostics_path}")
    print()
    print("Open an animation with ASE GUI, for example:")
    if diagnostics["selected_modes"]:
        print(
            "  ase gui "
            + diagnostics["selected_modes"][0]["trajectory"]
        )
    print()
    print("No new force calculations were performed.")


if __name__ == "__main__":
    main()
