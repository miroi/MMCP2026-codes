#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Electronic adsorption energies for Hg on Au(111)."""

from __future__ import annotations

import numpy as np
from ase import Atoms

from constants import EV_TO_KJ


def isolated_hg_energy(calc):
    hg = Atoms("Hg", positions=[[0.0, 0.0, 0.0]])
    hg.calc = calc
    return float(hg.get_potential_energy())


def adsorption_energy(E_HgAu, E_Au, E_Hg):
    values = [E_HgAu, E_Au, E_Hg]
    if not all(np.isfinite(x) for x in values):
        raise ValueError("Non-finite energy detected")

    return float(E_HgAu - E_Au - E_Hg)


def ev_to_kjmol(energy):
    return float(energy) * EV_TO_KJ


def print_energy_summary(E_HgAu, E_Au, E_Hg, E_ads):
    print("\n" + "=" * 70)
    print("ENERGY SUMMARY")
    print("=" * 70)
    print(f"E(Hg/Au)  : {E_HgAu: .10f} eV")
    print(f"E(Au)     : {E_Au: .10f} eV")
    print(f"E(Hg)     : {E_Hg: .10f} eV")
    print(f"ΔE_ads    : {E_ads: .10f} eV")
    print(f"          : {ev_to_kjmol(E_ads): .4f} kJ/mol")
    print("=" * 70)
