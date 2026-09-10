#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Finite-displacement Gamma-point vibrational analysis for Hg/Au(111).

Geometry-relaxation constraints and vibrational constraints are deliberately
separate.  During a vibrational calculation the full optimized structure is
passed to the calculator, so forces retain interactions with atoms that are
frozen in the vibrational coordinate space.  Only the selected mobile atoms
are displaced and therefore only their Hessian block is diagonalized.

The key control is n_frozen_layers_vibrations:
    - it counts bottom Au layers frozen for vibrations only;
    - it may equal the total number of Au layers;
    - Hg is never frozen by this Au-layer setting;
    - the original geometry-relaxation constraints are never modified.
"""

from __future__ import annotations

import gc
import json
import os
import shutil

import numpy as np
from ase.constraints import FixAtoms
from ase.vibrations import Vibrations

from constants import ASE_TO_CM


def clear_memory():
    gc.collect()
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
            try:
                torch.cuda.ipc_collect()
            except Exception:
                pass
    except Exception:
        pass


def get_frozen_mask(atoms):
    """Return atoms frozen by the geometry-relaxation FixAtoms constraints."""
    frozen = np.zeros(len(atoms), dtype=bool)
    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            if hasattr(constraint, "get_indices"):
                idx = constraint.get_indices()
            elif hasattr(constraint, "index"):
                idx = constraint.index
            elif hasattr(constraint, "indices"):
                idx = constraint.indices
            else:
                continue
            frozen[np.asarray(idx, dtype=int)] = True
    return frozen


def _au_layers(atoms):
    """Return Au atom indices and bottom-to-top layer z values.

    Relaxation can make atoms within one Au layer acquire slightly different
    z coordinates. Therefore, exact/rounded z values must not be used to
    identify layers after optimization. Au(111) layers are separated by
    roughly 2.35 A here, whereas intra-layer relaxation is much smaller.
    Cluster sorted Au z coordinates using a conservative 0.5 A gap threshold
    and use the mean z of each cluster as the layer position.
    """
    symbols = np.asarray(atoms.get_chemical_symbols())
    au_indices = np.where(symbols == "Au")[0]
    if au_indices.size == 0:
        return au_indices, np.empty(0, dtype=float)

    z = np.asarray(atoms.positions[au_indices, 2], dtype=float)
    order = np.argsort(z)
    sorted_indices = au_indices[order]
    z_sorted = z[order]

    # Keep the actual atom membership of each layer.  Using the layer mean
    # together with a tight isclose tolerance would fail after relaxation,
    # because atoms within one layer can have different z coordinates.
    layer_gap_threshold = 0.5
    split_points = np.where(np.diff(z_sorted) > layer_gap_threshold)[0] + 1
    groups = np.split(np.arange(len(sorted_indices)), split_points)

    layer_indices = [sorted_indices[group] for group in groups]
    layer_z = np.array(
        [np.mean(z_sorted[group]) for group in groups], dtype=float
    )
    return au_indices, layer_z, layer_indices


def get_vibrational_frozen_mask(atoms, n_frozen_layers_vibrations=None):
    """
    Build the vibrational frozen mask without changing atoms.constraints.

    Existing geometry-frozen atoms remain frozen.  In addition, the requested
    number of bottom Au layers are frozen for vibrational coordinates only.
    Non-Au atoms, including Hg, are not frozen by the layer setting.
    """
    geometry_frozen = get_frozen_mask(atoms)
    frozen = geometry_frozen.copy()

    if n_frozen_layers_vibrations is None:
        return frozen

    n_layers = int(n_frozen_layers_vibrations)
    if n_layers < 0:
        raise ValueError("n_frozen_layers_vibrations must be >= 0")

    au_indices, layers, layer_indices = _au_layers(atoms)
    if n_layers > len(layers):
        raise ValueError(
            f"n_frozen_layers_vibrations={n_layers} exceeds the "
            f"{len(layers)} Au layers present in the structure"
        )

    if n_layers == 0:
        return frozen

    # Freeze by layer membership, not by comparing relaxed z coordinates to
    # a layer-average z.  This guarantees that all atoms belonging to a
    # relaxed Au layer are selected, even when their z values differ.
    for layer_atoms in layer_indices[:n_layers]:
        frozen[np.asarray(layer_atoms, dtype=int)] = True
    return frozen


def _make_vibrational_selection(atoms, n_frozen_layers_vibrations=None):
    """Return the full-system vibrational mask and selected mobile indices."""
    frozen = get_vibrational_frozen_mask(
        atoms, n_frozen_layers_vibrations=n_frozen_layers_vibrations
    )
    mobile_indices = np.where(~frozen)[0]
    return mobile_indices, frozen


def _gamma_dynamical_matrix_from_hessian(atoms, hessian):
    """Mass-weight a selected-atom Hessian in eV/A^2."""
    H = np.asarray(hessian, dtype=float)
    n = len(atoms)
    expected = 3 * n
    if H.shape != (expected, expected):
        raise ValueError(
            f"Vibrational Hessian has shape {H.shape}, "
            f"expected {(expected, expected)}"
        )

    masses = np.asarray(atoms.get_masses(), dtype=float)
    if np.any(masses <= 0):
        raise ValueError("All vibrational atoms must have positive masses")

    m_inv_sqrt = np.repeat(masses ** -0.5, 3)
    D = H * np.outer(m_inv_sqrt, m_inv_sqrt)
    return 0.5 * (D + D.T)


def _gamma_modes(phon_atoms, hessian):
    D = _gamma_dynamical_matrix_from_hessian(phon_atoms, hessian)
    eigenvalues, eigenvectors = np.linalg.eigh(D)
    signed_cm = (
        np.sign(eigenvalues)
        * np.sqrt(np.abs(eigenvalues))
        * ASE_TO_CM
    )
    return eigenvalues, signed_cm, eigenvectors


def _hg_projections(phon_atoms, eigenvectors):
    symbols = phon_atoms.get_chemical_symbols()
    hg_atoms = np.array(
        [i for i, symbol in enumerate(symbols) if symbol == "Hg"],
        dtype=int,
    )

    n_modes = eigenvectors.shape[1]
    projection = np.zeros(n_modes, dtype=float)
    if len(hg_atoms) == 0:
        return hg_atoms, projection

    hg_dof = np.concatenate(
        [np.arange(3 * i, 3 * i + 3) for i in hg_atoms]
    )
    for mode in range(n_modes):
        vector = eigenvectors[:, mode]
        norm = np.sum(np.abs(vector) ** 2)
        if norm > 0:
            projection[mode] = float(
                np.sum(np.abs(vector[hg_dof]) ** 2) / norm
            )
    return hg_atoms, projection


def _empty_result(atoms, mobile_indices, frozen, cutoff, displacement,
                  n_frozen_layers_vibrations, hg_projection_threshold):
    """Return a valid zero-mode result for a fully frozen clean slab."""
    phon_atoms = atoms[mobile_indices].copy()
    phon_atoms.set_constraint([])
    empty = np.empty(0, dtype=float)
    metadata = {
        "n_original_atoms": int(len(atoms)),
        "n_frozen_atoms": int(np.sum(frozen)),
        "n_mobile_atoms": int(len(mobile_indices)),
        "mobile_indices": mobile_indices.tolist(),
        "hg_atoms_in_phonon_system": [],
        "n_raw_modes": 0,
        "n_kept_modes": 0,
        "n_positive_modes": 0,
        "n_significant_imaginary_modes": 0,
        "frequency_cutoff_cm-1": float(cutoff),
        "phonon_displacement_A": float(displacement),
        "hg_projection_threshold": float(hg_projection_threshold),
        "n_frozen_layers_vibrations": (
            None if n_frozen_layers_vibrations is None
            else int(n_frozen_layers_vibrations)
        ),
        "force_calculations_include_frozen_atoms": True,
        "phonon_coordinate_space": "mobile atoms only",
        "phonon_backend": "ASE Vibrations",
    }
    return {
        "frequencies_cm": empty,
        "signed_frequencies_cm": empty,
        "eigenvalues": empty,
        "eigenvectors": np.empty((0, 0), dtype=float),
        "hg_projection": empty,
        "hg_mode_mask": np.zeros(0, dtype=bool),
        "phonon": None,
        "phonon_atoms": phon_atoms,
        "mobile_indices": mobile_indices,
        "frozen_mask": frozen,
        "imaginary_frequencies_cm": empty,
        "positive_frequencies_cm": empty,
        "frequency_cutoff_cm": float(cutoff),
        "metadata": metadata,
    }


def _custom_hessian(atoms, calc, mobile_indices, displacement, nfree, symmetrize=True):
    """Build the mobile-coordinate Hessian directly from finite-difference forces."""
    n_mobile = len(mobile_indices)
    ndof = 3 * n_mobile
    hessian = np.zeros((ndof, ndof), dtype=float)
    work = atoms.copy()
    work.calc = calc

    for j, atom_index in enumerate(mobile_indices):
        for axis in range(3):
            col = 3 * j + axis
            coordinate = 3 * atom_index + axis

            if nfree == 2:
                plus = work.copy()
                minus = work.copy()
                plus.calc = calc
                minus.calc = calc
                plus.positions[atom_index, axis] += displacement
                minus.positions[atom_index, axis] -= displacement
                f_plus = np.asarray(plus.get_forces(), dtype=float)[mobile_indices]
                f_minus = np.asarray(minus.get_forces(), dtype=float)[mobile_indices]
                hessian[:, col] = -(
                    (f_plus - f_minus).reshape(-1) / (2.0 * displacement)
                )
            else:
                pp = work.copy()
                p = work.copy()
                m = work.copy()
                mm = work.copy()
                for obj in (pp, p, m, mm):
                    obj.calc = calc
                pp.positions[atom_index, axis] += 2.0 * displacement
                p.positions[atom_index, axis] += displacement
                m.positions[atom_index, axis] -= displacement
                mm.positions[atom_index, axis] -= 2.0 * displacement
                f_pp = np.asarray(pp.get_forces(), dtype=float)[mobile_indices]
                f_p = np.asarray(p.get_forces(), dtype=float)[mobile_indices]
                f_m = np.asarray(m.get_forces(), dtype=float)[mobile_indices]
                f_mm = np.asarray(mm.get_forces(), dtype=float)[mobile_indices]
                # Fourth-order central derivative of force, followed by the
                # minus sign converting dF/dx into the force-constant Hessian.
                hessian[:, col] = -(
                    (f_mm - 8.0 * f_m + 8.0 * f_p - f_pp).reshape(-1)
                    / (12.0 * displacement)
                )

    # Numerical finite differences need not be exactly symmetric.
    if symmetrize:
        hessian = 0.5 * (hessian + hessian.T)
    return hessian


def _calculate_phonons_custom(
    atoms, calc, mobile_indices, frozen, displacement, prefix,
    frequency_cutoff_cm, n_frozen_layers_vibrations,
    hg_projection_threshold, nfree, symmetrize, fail_on_imaginary, verbose
):
    """Custom finite-difference Hessian/mode calculation."""
    print("Vibrational backend: custom project Hessian")
    print("Custom Hessian: direct finite-difference forces")
    clear_memory()
    hessian = _custom_hessian(
        atoms, calc, mobile_indices, displacement, int(nfree),
        symmetrize=bool(symmetrize)
    )
    clear_memory()

    phon_atoms = atoms[mobile_indices].copy()
    phon_atoms.set_constraint([])
    eigenvalues, signed_cm, eigenvectors = _gamma_modes(phon_atoms, hessian)

    cutoff = float(frequency_cutoff_cm)
    keep = signed_cm >= cutoff if cutoff > 0 else signed_cm > 0.0
    frequencies_cm = signed_cm[keep]
    imaginary = signed_cm[signed_cm < -cutoff] if cutoff > 0 else signed_cm[signed_cm < 0]

    if len(imaginary) > 0:
        message = (
            f"{len(imaginary)} significant imaginary phonon mode(s) detected. "
            "The structure is not a stable harmonic minimum."
        )
        if fail_on_imaginary:
            raise RuntimeError(message)
        if verbose:
            print(f"WARNING: {message}")
    positive = signed_cm[signed_cm > cutoff] if cutoff > 0 else signed_cm[signed_cm > 0]

    print(f"Raw Gamma modes: {len(signed_cm)}")
    print(f"Positive modes kept (nu >= {frequency_cutoff_cm:g} cm^-1): {len(frequencies_cm)}")
    if len(imaginary):
        print(f"Significant imaginary modes: {len(imaginary)}")
        print("Most negative modes (cm^-1):", np.sort(imaginary)[:10])
    else:
        print("No significant imaginary modes.")
    if len(positive):
        print(f"Positive frequency range: {positive.min():.3f}–{positive.max():.3f} cm^-1")
        print("Lowest positive modes (cm^-1):", np.array2string(np.sort(positive)[:10], precision=4))

    hg_atoms, hg_projection = _hg_projections(phon_atoms, eigenvectors)
    if len(hg_atoms):
        print("Hg atom index in mobile phonon coordinates:", hg_atoms.tolist())
        hg_mode_mask = keep & (signed_cm > 0.0) & (
            hg_projection >= float(hg_projection_threshold)
        )
        hg_modes = signed_cm[hg_mode_mask]
        print(f"Hg-dominated modes (projection >= {hg_projection_threshold:.3f}):", len(hg_modes))
        if len(hg_modes):
            print("Hg-dominated frequencies (cm^-1):", np.array2string(np.sort(hg_modes), precision=3))
    else:
        hg_mode_mask = np.zeros(len(signed_cm), dtype=bool)

    np.savetxt(prefix + "_frequencies_cm-1.dat", frequencies_cm, fmt="%.10f")
    np.savetxt(prefix + "_signed_frequencies_cm-1.dat", signed_cm, fmt="%.10f")
    np.savetxt(prefix + "_gamma_eigenvalues.dat", eigenvalues, fmt="%.16e")
    diagnostics = np.column_stack((
        np.arange(len(signed_cm)), signed_cm, hg_projection,
        keep.astype(int), hg_mode_mask.astype(int)
    ))
    np.savetxt(prefix + "_mode_diagnostics.dat", diagnostics,
               header="mode frequency_cm-1 hg_projection kept_positive hg_dominated",
               fmt=["%d", "%.10f", "%.10f", "%d", "%d"])

    metadata = {
        "n_original_atoms": int(len(atoms)),
        "n_frozen_atoms": int(np.sum(frozen)),
        "n_mobile_atoms": int(len(mobile_indices)),
        "mobile_indices": mobile_indices.tolist(),
        "hg_atoms_in_phonon_system": hg_atoms.tolist(),
        "n_raw_modes": int(len(signed_cm)),
        "n_kept_modes": int(len(frequencies_cm)),
        "n_positive_modes": int(len(positive)),
        "n_significant_imaginary_modes": int(len(imaginary)),
        "frequency_cutoff_cm-1": cutoff,
        "phonon_displacement_A": float(displacement),
        "hg_projection_threshold": float(hg_projection_threshold),
        "n_frozen_layers_vibrations": None if n_frozen_layers_vibrations is None else int(n_frozen_layers_vibrations),
        "nfree": int(nfree),
        "symmetrize": bool(symmetrize),
        "fail_on_imaginary": bool(fail_on_imaginary),
        "force_calculations_include_frozen_atoms": True,
        "phonon_coordinate_space": "mobile atoms only",
        "phonon_backend": "custom project Hessian",
    }
    with open(prefix + "_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    clear_memory()
    return {
        "frequencies_cm": frequencies_cm,
        "signed_frequencies_cm": signed_cm,
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,
        "hg_projection": hg_projection,
        "hg_mode_mask": hg_mode_mask,
        "phonon": None,
        "phonon_atoms": phon_atoms,
        "mobile_indices": mobile_indices,
        "frozen_mask": frozen,
        "imaginary_frequencies_cm": imaginary,
        "positive_frequencies_cm": positive,
        "frequency_cutoff_cm": cutoff,
        "metadata": metadata,
    }


def calculate_phonons(
    atoms,
    calc,
    displacement=0.005,
    prefix="phonons",
    verbose=True,
    frequency_cutoff_cm=1.0,
    symmetrize=True,
    n_frozen_layers_vibrations=None,
    hg_projection_threshold=0.20,
    nfree=2,
    method="custom",
    fail_on_imaginary=True,
):
    """Run a Gamma-point vibrational calculation using the selected backend."""

    method = str(method).strip().lower()
    if method not in {"custom", "ase"}:
        raise ValueError("phonons method must be either custom or ase")
    if displacement <= 0:
        raise ValueError("displacement must be > 0")
    if frequency_cutoff_cm < 0:
        raise ValueError("frequency_cutoff_cm must be >= 0")
    if not (0.0 < hg_projection_threshold <= 1.0):
        raise ValueError("hg_projection_threshold must be in (0, 1]")
    if int(nfree) not in (2, 4):
        raise ValueError("nfree must be 2 or 4")
    fail_on_imaginary = bool(fail_on_imaginary)

    print("\nPHONON CALCULATION")
    print("-" * 70)

    clear_memory()
    mobile_indices, frozen = _make_vibrational_selection(
        atoms, n_frozen_layers_vibrations=n_frozen_layers_vibrations
    )

    symbols = atoms.get_chemical_symbols()
    au_indices, layers, layer_indices = _au_layers(atoms)
    frozen_au = np.sum(frozen[au_indices]) if len(au_indices) else 0
    mobile_au = np.sum(~frozen[au_indices]) if len(au_indices) else 0
    mobile_hg = np.sum((~frozen) & (np.asarray(symbols) == "Hg"))

    print(f"Requested vibrational backend: {method}")
    print(f"Original atoms: {len(atoms)}")
    print(f"Frozen atoms: {int(np.sum(frozen))}")
    print(f"Mobile atoms: {len(mobile_indices)}")
    print(f"Vibrational frozen Au atoms: {int(frozen_au)}")
    print(f"Vibrationally mobile Au atoms: {int(mobile_au)}")
    print(f"Vibrationally mobile Hg atoms: {int(mobile_hg)}")
    if n_frozen_layers_vibrations is not None:
        print("Vibrational frozen Au layers:", int(n_frozen_layers_vibrations), f"of {len(layers)}")
    print("Mobile composition:", " ".join(
        f"{s}:{sum(1 for i in mobile_indices if symbols[i] == s)}"
        for s in sorted(set(symbols[i] for i in mobile_indices))
    ) if len(mobile_indices) else "none")
    print("Force calculations include frozen atoms: YES")
    print("Phonon coordinate space: mobile atoms only")
    print("Hg is included in the phonon system:", "YES" if mobile_hg else "NO")
    print("Phonon configuration:")
    print(f"  phonon_displacement     : {displacement:.6f} Å")
    print(f"  frequency_cutoff_cm    : {frequency_cutoff_cm:.6f} cm^-1")
    print(f"  n_frozen_layers_vibrations: {n_frozen_layers_vibrations}")
    print(f"  hg_projection_threshold: {hg_projection_threshold:.3f}")
    print(f"  symmetrize              : {bool(symmetrize)}")
    print(f"  fail_on_imaginary       : {fail_on_imaginary}")

    if len(mobile_indices) == 0:
        if any(symbol == "Hg" for symbol in symbols):
            raise ValueError("The vibrational selection froze all atoms of a system containing Hg. Hg must remain mobile.")
        print("All atoms are frozen for vibrations; no vibrational modes are present.")
        metadata = {
            "n_original_atoms": int(len(atoms)),
            "n_frozen_atoms": int(np.sum(frozen)),
            "n_mobile_atoms": 0,
            "mobile_indices": [],
            "hg_atoms_in_phonon_system": [],
            "n_raw_modes": 0,
            "n_kept_modes": 0,
            "n_positive_modes": 0,
            "n_significant_imaginary_modes": 0,
            "frequency_cutoff_cm-1": float(frequency_cutoff_cm),
            "phonon_displacement_A": float(displacement),
            "hg_projection_threshold": float(hg_projection_threshold),
            "n_frozen_layers_vibrations": None if n_frozen_layers_vibrations is None else int(n_frozen_layers_vibrations),
            "force_calculations_include_frozen_atoms": True,
            "phonon_coordinate_space": "mobile atoms only",
            "phonon_backend": method,
        }
        with open(prefix + "_metadata.json", "w", encoding="utf-8") as handle:
            json.dump(metadata, handle, indent=2)
        for suffix in ("_frequencies_cm-1.dat","_signed_frequencies_cm-1.dat","_gamma_eigenvalues.dat","_mode_diagnostics.dat"):
            open(prefix + suffix, "w", encoding="utf-8").close()
        return {
            "frequencies_cm": np.empty(0), "signed_frequencies_cm": np.empty(0),
            "eigenvalues": np.empty(0), "eigenvectors": np.empty((0,0)),
            "hg_projection": np.empty(0), "hg_mode_mask": np.zeros(0, dtype=bool),
            "phonon": None, "phonon_atoms": atoms[mobile_indices].copy(),
            "mobile_indices": mobile_indices, "frozen_mask": frozen,
            "imaginary_frequencies_cm": np.empty(0), "positive_frequencies_cm": np.empty(0),
            "frequency_cutoff_cm": float(frequency_cutoff_cm), "metadata": metadata,
        }

    if method == "custom":
        return _calculate_phonons_custom(
            atoms, calc, mobile_indices, frozen, displacement, prefix,
            frequency_cutoff_cm, n_frozen_layers_vibrations,
            hg_projection_threshold, nfree, symmetrize,
            fail_on_imaginary, verbose
        )

    print("Vibrational backend: ASE Vibrations")
    if os.path.exists(prefix):
        shutil.rmtree(prefix)
    vib = Vibrations(
        atoms, indices=mobile_indices.tolist(), delta=displacement,
        nfree=int(nfree), name=prefix,
    )
    try:
        vib.run()
        vib.read(method="standard", direction="central")
        hessian = np.asarray(vib.H, dtype=float)
    finally:
        clear_memory()

    phon_atoms = atoms[mobile_indices].copy()
    phon_atoms.set_constraint([])
    if symmetrize:
        hessian = 0.5 * (hessian + hessian.T)
    eigenvalues, signed_cm, eigenvectors = _gamma_modes(phon_atoms, hessian)
    cutoff = float(frequency_cutoff_cm)
    keep = signed_cm >= cutoff if cutoff > 0 else signed_cm > 0.0
    frequencies_cm = signed_cm[keep]
    imaginary = signed_cm[signed_cm < -cutoff] if cutoff > 0 else signed_cm[signed_cm < 0]
    if len(imaginary) > 0:
        message = (
            f"{len(imaginary)} significant imaginary phonon mode(s) detected. "
            "The structure is not a stable harmonic minimum."
        )
        if fail_on_imaginary:
            raise RuntimeError(message)
        if verbose:
            print(f"WARNING: {message}")
    positive = signed_cm[signed_cm > cutoff] if cutoff > 0 else signed_cm[signed_cm > 0]

    print(f"Raw Gamma modes: {len(signed_cm)}")
    print(f"Positive modes kept (nu >= {frequency_cutoff_cm:g} cm^-1): {len(frequencies_cm)}")
    print(f"WARNING: {len(imaginary)} significant imaginary modes." if len(imaginary) else "No significant imaginary modes.")
    if len(positive):
        print(f"Positive frequency range: {positive.min():.3f}–{positive.max():.3f} cm^-1")
        print("Lowest positive modes (cm^-1):", np.array2string(np.sort(positive)[:10], precision=4))

    hg_atoms, hg_projection = _hg_projections(phon_atoms, eigenvectors)
    if len(hg_atoms):
        print("Hg atom index in mobile phonon coordinates:", hg_atoms.tolist())
        hg_mode_mask = keep & (signed_cm > 0.0) & (hg_projection >= float(hg_projection_threshold))
        hg_modes = signed_cm[hg_mode_mask]
        print(f"Hg-dominated modes (projection >= {hg_projection_threshold:.3f}):", len(hg_modes))
        if len(hg_modes):
            print("Hg-dominated frequencies (cm^-1):", np.array2string(np.sort(hg_modes), precision=3))
    else:
        hg_mode_mask = np.zeros(len(signed_cm), dtype=bool)

    np.savetxt(prefix + "_frequencies_cm-1.dat", frequencies_cm, fmt="%.10f")
    np.savetxt(prefix + "_signed_frequencies_cm-1.dat", signed_cm, fmt="%.10f")
    np.savetxt(prefix + "_gamma_eigenvalues.dat", eigenvalues, fmt="%.16e")
    diagnostics = np.column_stack((np.arange(len(signed_cm)), signed_cm, hg_projection, keep.astype(int), hg_mode_mask.astype(int)))
    np.savetxt(prefix + "_mode_diagnostics.dat", diagnostics,
               header="mode frequency_cm-1 hg_projection kept_positive hg_dominated",
               fmt=["%d","%.10f","%.10f","%d","%d"])
    metadata = {
        "n_original_atoms": int(len(atoms)), "n_frozen_atoms": int(np.sum(frozen)),
        "n_mobile_atoms": int(len(mobile_indices)), "mobile_indices": mobile_indices.tolist(),
        "hg_atoms_in_phonon_system": hg_atoms.tolist(), "n_raw_modes": int(len(signed_cm)),
        "n_kept_modes": int(len(frequencies_cm)), "n_positive_modes": int(len(positive)),
        "n_significant_imaginary_modes": int(len(imaginary)),
        "frequency_cutoff_cm-1": cutoff, "phonon_displacement_A": float(displacement),
        "hg_projection_threshold": float(hg_projection_threshold),
        "n_frozen_layers_vibrations": None if n_frozen_layers_vibrations is None else int(n_frozen_layers_vibrations),
        "nfree": int(nfree),
        "symmetrize": bool(symmetrize),
        "fail_on_imaginary": bool(fail_on_imaginary),
        "force_calculations_include_frozen_atoms": True,
        "phonon_coordinate_space": "mobile atoms only", "phonon_backend": "ASE Vibrations",
    }
    with open(prefix + "_metadata.json", "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
    clear_memory()
    return {
        "frequencies_cm": frequencies_cm, "signed_frequencies_cm": signed_cm,
        "eigenvalues": eigenvalues, "eigenvectors": eigenvectors,
        "hg_projection": hg_projection, "hg_mode_mask": hg_mode_mask,
        "phonon": vib, "phonon_atoms": phon_atoms, "mobile_indices": mobile_indices,
        "frozen_mask": frozen, "imaginary_frequencies_cm": imaginary,
        "positive_frequencies_cm": positive, "frequency_cutoff_cm": cutoff,
        "metadata": metadata,
    }


def get_hg_modes(phonon_result, threshold=0.20, positive_only=True):
    """Return Hg-dominated modes from a calculate_phonons result."""
    freq = np.asarray(phonon_result["signed_frequencies_cm"], dtype=float)
    proj = np.asarray(phonon_result["hg_projection"], dtype=float)
    mask = proj >= float(threshold)
    if positive_only:
        mask &= freq > 0.0
    return freq[mask]


def get_hg_projected_frequencies(phonon_result, threshold=0.20):
    return get_hg_modes(
        phonon_result,
        threshold=threshold,
        positive_only=True,
    )
