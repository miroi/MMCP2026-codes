#!/usr/bin/env python3
"""
solve_for_Ta.py — corrected compact Eichler–Zvara deposition-temperature solver.

The project equation is

    ΔH_des/(R Ta) = ΔS_des/R
        + ln[s0*nu*sqrt(2*pi*M)*Q*g/(s*Ta**(3/2))]

For this compact form, M is entered explicitly in g/mol, following the
thermochromatographic notation used by the project/reference material.

The original project had a 60x flow conversion error:
    0.002 sccm = 3.333333333e-11 m^3/s
not 2e-9 m^3/s.

This script also reports the published mobile-adsorption entropy as a
diagnostic, using the single-particle mass m = M/N_A.
"""

from __future__ import annotations
import math
import os
import sys
from typing import Dict
from scipy.optimize import brentq

AVOGADRO = 6.02214076e23


def read_parameters(filename: str = "input_params.txt") -> Dict[str, float]:
    if not os.path.exists(filename):
        raise FileNotFoundError(f"Input file not found: {filename}")
    params: Dict[str, float] = {}
    with open(filename, "r", encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            if "=" not in line:
                continue
            name, value = line.split("=", 1)
            try:
                params[name.strip()] = float(value.strip())
            except ValueError as exc:
                raise ValueError(f"Could not parse: {line!r}") from exc
    if not params:
        raise ValueError(f"No numeric parameters found in {filename}")
    return params


def sccm_to_m3_s(q_sccm: float) -> float:
    """1 sccm = 1 cm^3/min = 1e-6/60 m^3/s."""
    return q_sccm * 1.0e-6 / 60.0


def m3_s_to_sccm(q_m3_s: float) -> float:
    return q_m3_s * 60.0e6


def mobile_adsorption_entropy(
    T_K: float,
    molar_mass_g_mol: float,
    nu_s_inv: float,
    standard_length_m: float = 0.01,
) -> float:
    """
    Published mobile-adsorption entropy diagnostic:

    ΔS = R ln[(1/L0)(1/nu_B)sqrt(k_B*T/(2*pi*m))] + R/2

    with L0 = 1 cm and m = mass of one particle.
    """
    k_B = 1.380649e-23
    R = 8.31446261815324
    m_particle_kg = molar_mass_g_mol * 1.0e-3 / AVOGADRO
    arg = (
        (1.0 / standard_length_m)
        * (1.0 / nu_s_inv)
        * math.sqrt(k_B * T_K / (2.0 * math.pi * m_particle_kg))
    )
    return R * math.log(arg) + 0.5 * R


def main() -> None:
    filename = sys.argv[1] if len(sys.argv) > 1 else "input_params.txt"
    p = read_parameters(filename)

    R = p["R_J_mol_K"]
    H_ads = p["delta_H_ads_J_mol"]
    S_ads = p["delta_S_ads_J_mol_K"]

    M_g_mol = p["M_adsorbate_g_mol"]
    nu = p["nu_phonon_s_inv"]
    s0 = p["s0_m2"]
    Q_sccm = p["Q_sccm"]
    g = p["temperature_gradient_K_m"]
    d = p["column_diameter_m"]

    Q = sccm_to_m3_s(Q_sccm)
    s = math.pi * (d / 2.0) ** 2

    if min(M_g_mol, nu, s0, Q, g, s) <= 0:
        raise ValueError("Mass, frequency, area, flow and gradient must be positive.")

    factor = s0 * nu * math.sqrt(2.0 * math.pi * M_g_mol) * Q * g / s

    H_des = -H_ads
    S_des = -S_ads

    def residual(T: float) -> float:
        if T <= 0:
            return float("inf")
        return H_des / (R * T) - S_des / R - math.log(factor / T**1.5)

    Tmin = p.get("T_min_K", 50.0)
    Tmax = p.get("T_max_K", 1000.0)

    # Bracket the root rather than relying on an arbitrary fsolve guess.
    grid = [Tmin + (Tmax - Tmin) * i / 2000.0 for i in range(2001)]
    root = None
    for a, b in zip(grid[:-1], grid[1:]):
        fa, fb = residual(a), residual(b)
        if fa == 0:
            root = a
            break
        if fa * fb < 0:
            root = brentq(residual, a, b, xtol=1e-10, rtol=1e-12)
            break
    if root is None:
        raise RuntimeError(f"No root found in [{Tmin}, {Tmax}] K")

    T_eq = H_ads / S_ads

    print("=" * 72)
    print("EICHLER–ZVARA COMPACT DEPOSITION-TEMPERATURE CALCULATION")
    print("=" * 72)
    print(f"Input file:                 {filename}")
    print()
    print("Thermodynamic parameters")
    print(f"  ΔH_ads                    {H_ads/1000:.6f} kJ/mol")
    print(f"  ΔS_ads                    {S_ads:.6f} J/(mol K)")
    print(f"  ΔH_des                    {H_des/1000:.6f} kJ/mol")
    print(f"  ΔS_des                    {S_des:.6f} J/(mol K)")
    print(f"  R                         {R:.9f} J/(mol K)")
    print()
    print("Experimental/model parameters")
    print(f"  M                         {M_g_mol:.6f} g/mol")
    print(f"  nu_B                      {nu:.6e} s^-1")
    print(f"  s0                        {s0:.6e} m^2")
    print(f"  Q                         {Q_sccm:.6f} sccm")
    print(f"  Q                         {Q:.9e} m^3/s")
    print(f"  g                         {g:.6f} K/m")
    print(f"  d                         {d:.6f} m")
    print(f"  column area s             {s:.9e} m^2")
    print()
    print(f"  compact prefactor         {factor:.9e}")
    print()
    print(f"Solved deposition Ta        {root:.6f} K")
    print(f"                            {root - 273.15:.6f} °C")
    print(f"Equation residual           {residual(root):.3e}")
    print()
    print(f"Constant-H,S ΔG=0 T        {T_eq:.6f} K")
    print(f"                            {T_eq - 273.15:.6f} °C")
    print()
    print("Mobile-adsorption entropy diagnostic")
    print(f"  ΔS_mobile(Ta)             "
          f"{mobile_adsorption_entropy(root, M_g_mol, nu):.6f} J/(mol K)")
    print()
    print("Flow conversion check")
    print(f"  {Q_sccm:.6f} sccm = {Q:.12e} m^3/s")
    print("  The old 2e-9 m^3/s value corresponds to 0.12 sccm.")
    print("=" * 72)


if __name__ == "__main__":
    main()
