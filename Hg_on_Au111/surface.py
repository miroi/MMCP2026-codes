#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Au(111) construction, constraints, and Hg adsorption geometry."""

from __future__ import annotations

import numpy as np
from ase import Atoms
from ase.build import add_adsorbate, fcc111
from ase.constraints import FixAtoms


def _au_layer_z(atoms):
    au_indices = [
        i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Au"
    ]
    if not au_indices:
        raise ValueError("No Au atoms found")
    z = np.asarray(atoms.positions[au_indices, 2], dtype=float)
    return au_indices, np.unique(np.round(z, 6))


def freeze_bottom_layers(atoms, n_layers):
    """Freeze only Au atoms in the bottom n_layers."""
    au_indices, layers = _au_layer_z(atoms)

    if n_layers < 0 or n_layers >= len(layers):
        raise ValueError(
            f"Invalid n_layers={n_layers}; slab has {len(layers)} Au layers"
        )

    frozen_z = layers[:n_layers]
    mask = np.zeros(len(atoms), dtype=bool)

    for i in au_indices:
        if np.any(np.isclose(atoms.positions[i, 2], frozen_z, atol=1e-5)):
            mask[i] = True

    atoms.set_constraint(FixAtoms(mask=mask))
    return atoms


def get_frozen_mask(atoms):
    """Return a reliable boolean mask from all FixAtoms constraints."""
    mask = np.zeros(len(atoms), dtype=bool)

    for constraint in atoms.constraints:
        if isinstance(constraint, FixAtoms):
            if hasattr(constraint, "get_indices"):
                idx = np.asarray(constraint.get_indices(), dtype=int)
            elif hasattr(constraint, "index"):
                idx = np.asarray(constraint.index, dtype=int)
            elif hasattr(constraint, "indices"):
                idx = np.asarray(constraint.indices, dtype=int)
            else:
                continue
            mask[idx] = True

    return mask


def get_mobile_indices(atoms):
    return np.where(~get_frozen_mask(atoms))[0].tolist()


def build_au111(config):
    slab = fcc111(
        "Au",
        size=config.surface_size,
        a=config.lattice_constant,
        vacuum=config.vacuum_size,
        orthogonal=True,
    )
    slab = freeze_bottom_layers(slab, config.n_frozen_layers)
    return slab


def get_adsorption_site_position(slab, site="fcc", height=3.0):
    temp = slab.copy()
    temp.set_constraint([])
    add_adsorbate(temp, "Hg", height, position=site)
    pos = temp.positions[-1]
    return float(pos[0]), float(pos[1])


def get_hg_surface_distance(atoms):
    """Return the vertical Hg-to-top-Au-surface distance in Angstrom."""
    symbols = atoms.get_chemical_symbols()
    au_indices = [i for i, s in enumerate(symbols) if s == "Au"]
    hg_indices = [i for i, s in enumerate(symbols) if s == "Hg"]

    if not au_indices:
        raise ValueError("No Au atoms found")
    if len(hg_indices) != 1:
        raise ValueError(
            f"Expected exactly one Hg atom, found {len(hg_indices)}"
        )

    top_au_z = float(np.max(atoms.positions[au_indices, 2]))
    hg_z = float(atoms.positions[hg_indices[0], 2])
    return hg_z - top_au_z


def add_hg_atom(slab, config, site_position):
    atoms = slab.copy()
    atoms.set_constraint([])

    au_indices = [
        i for i, s in enumerate(atoms.get_chemical_symbols()) if s == "Au"
    ]
    top_z = np.max(atoms.positions[au_indices, 2])

    x, y = site_position
    hg = Atoms("Hg", positions=[[x, y, top_z + config.hg_height]])
    combined = atoms + hg

    combined = freeze_bottom_layers(
        combined, config.n_frozen_layers
    )

    return combined


def print_surface_details(atoms, n_layers):
    symbols = atoms.get_chemical_symbols()
    au_indices = [i for i, s in enumerate(symbols) if s == "Au"]

    layers = np.unique(np.round(atoms.positions[au_indices, 2], 6))

    print("\nSURFACE DETAILS")
    print("-" * 60)
    print("Total atoms:", len(atoms))
    print("Au atoms:", len(au_indices))
    print("Au layers:", len(layers))
    print("Frozen layers:", n_layers)

    for i, z in enumerate(layers):
        count = np.sum(
            np.isclose(atoms.positions[au_indices, 2], z, atol=1e-5)
        )
        state = "FROZEN" if i < n_layers else "FREE"
        print(
            f"Layer {i+1}: z={z:.6f} Å  {count} Au  {state}"
        )

    if "Hg" in symbols:
        hg_index = symbols.index("Hg")
        top_au = np.max(atoms.positions[au_indices, 2])
        print(
            "\nHg–Au surface distance:",
            f"{atoms.positions[hg_index,2] - top_au:.6f} Å",
        )
