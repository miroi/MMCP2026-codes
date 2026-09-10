#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Complete Hg adsorption workflow on Au(111), including G_ads(T)=0 search."""

from __future__ import annotations

import gc
import json
import math

import numpy as np

from calculator import clear_cuda_cache
from energies import (
    adsorption_energy,
    isolated_hg_energy,
    print_energy_summary,
)
from optimization import optimize_structure
from phonons import calculate_phonons
from surface import (
    add_hg_atom,
    build_au111,
    get_adsorption_site_position,
    get_hg_surface_distance,
    print_surface_details,
)
from thermodynamics import (
    calculate_adsorption_thermodynamics,
    print_ase_hgau_thermo_batch,
)


def _save_json(filename, data):
    """Save JSON-safe scalar/array results."""
    clean = {}
    for key, value in data.items():
        if isinstance(value, (np.floating, np.integer)):
            clean[key] = value.item()
        elif isinstance(value, np.ndarray):
            clean[key] = value.tolist()
        else:
            clean[key] = value

    with open(filename, "w", encoding="utf-8") as handle:
        json.dump(clean, handle, indent=2)


def _g_ads(Eads, Au_freq, HgAu_freq, temperature, pressure, include_zpe, reference_pressure=1.0, method="custom"):
    """Evaluate Delta G_ads at one temperature."""
    thermo = calculate_adsorption_thermodynamics(
        Eads,
        HgAu_freq,
        Au_freq,
        float(temperature),
        pressure,
        include_zpe=include_zpe,
        reference_pressure=reference_pressure,
        method=method,
    )
    return float(thermo["G_ads"]), thermo


def _thermo_scan(Eads, Au_freq, HgAu_freq, config):
    """Evaluate the thermodynamics on the configured temperature grid."""
    rows = []

    for T in config.temperatures:
        g_ads, thermo = _g_ads(
            Eads,
            Au_freq,
            HgAu_freq,
            T,
            config.pressure,
            config.include_zpe,
            config.gas_reference_pressure,
            config.thermo_method,
        )
        rows.append({"T_K": float(T), **thermo})

    filename = config.output_prefix + "_thermodynamics_scan.csv"
    keys = [
        "T_K",
        "E_ads",
        "ZPE_ads",
        "delta_U_vib_thermal",
    "delta_U_vib_total",
        "delta_F_vib",
        "S_gas",
        "S_ads_vib",
        "S_surface_vib",
        "S_ads",
        "H_gas_thermal",
        "H_ads",
        "G_ads",
    ]

    with open(filename, "w", encoding="utf-8") as handle:
        handle.write(",".join(keys) + "\n")
        for row in rows:
            handle.write(
                ",".join(f"{row.get(k, float('nan')):.12g}" for k in keys)
                + "\n"
            )

    # Always show the temperature scan in the terminal when mode=scan.
    # The scan reuses the already calculated phonon frequencies; no new
    # electronic-structure or phonon calculations are performed here.
    print("\n" + "=" * 70)
    print("THERMODYNAMICS TEMPERATURE SCAN")
    print("-" * 70)
    print(
        f"{'T (K)':>12}  {'ΔH_ads (eV)':>16}  "
        f"{'ΔS_ads (eV/K)':>18}  {'ΔG_ads (eV)':>16}"
    )
    print("-" * 70)
    for row in rows:
        print(
            f"{row['T_K']:12.2f}  "
            f"{row['H_ads']:16.8f}  "
            f"{row['S_ads']:18.8e}  "
            f"{row['G_ads']:16.8f}"
        )
    print("-" * 70)
    print(
        "Scan temperatures: "
        + ", ".join(f"{float(T):.2f}" for T in config.temperatures)
        + " K"
    )
    print(f"Scan output: {filename}")
    print("=" * 70)

    return rows


def _find_g_ads_zero(
    Eads,
    Au_freq,
    HgAu_freq,
    config,
    temperatures=None,
):
    """
    Find the temperature at which Delta G_ads(T) = 0.

    The search first checks the supplied temperature grid. If no adjacent
    grid points bracket a root, it searches a continuous interval using
    config.temperature_root_min/max_K.  A safeguarded bisection solve is
    then used, so no electronic-structure or phonon calculations are needed.
    """
    if temperatures is None:
        temperatures = np.asarray(config.temperatures, dtype=float)
    else:
        temperatures = np.asarray(temperatures, dtype=float)

    temperatures = temperatures[np.isfinite(temperatures)]
    temperatures = temperatures[temperatures > 0.0]

    if len(temperatures) < 2:
        raise ValueError(
            "At least two positive temperatures are required in the "
            "[calculation] temperatures setting to locate a ΔG_ads = 0 root."
        )

    # The existing [calculation] temperatures define the root-search interval.
    # No additional temperature-range keywords are required.
    tmin = float(np.min(temperatures))
    tmax = float(np.max(temperatures))

    if tmin <= 0.0 or tmax <= tmin:
        raise ValueError(
            "The [calculation] temperatures setting must contain at least "
            "two distinct positive temperatures."
        )

    # Include the configured scan points that fall inside the root interval.
    grid = temperatures[(temperatures >= tmin) & (temperatures <= tmax)]
    grid = np.unique(np.concatenate(([tmin], grid, [tmax])))
    grid.sort()

    values = []
    for T in grid:
        g, _ = _g_ads(
            Eads, Au_freq, HgAu_freq, T,
            config.pressure, config.include_zpe,
            config.gas_reference_pressure, config.thermo_method
        )
        values.append(g)
    values = np.asarray(values, dtype=float)

    # Exact grid-point root.
    exact = np.where(np.isclose(values, 0.0, atol=1e-14, rtol=0.0))[0]
    if len(exact):
        idx = int(exact[0])
        T = float(grid[idx])
        g, thermo = _g_ads(
            Eads, Au_freq, HgAu_freq, T,
            config.pressure, config.include_zpe,
            config.gas_reference_pressure, config.thermo_method
        )
        return {
            "found": True,
            "temperature_K": T,
            "G_ads_eV": g,
            "method": "exact_grid_point",
            "bracket_K": [T, T],
            "iterations": 0,
            "thermodynamics": thermo,
        }

    # Locate the first sign change in ascending temperature.
    bracket = None
    for i in range(len(grid) - 1):
        if values[i] * values[i + 1] < 0.0:
            bracket = (float(grid[i]), float(grid[i + 1]))
            break

    if bracket is None:
        return {
            "found": False,
            "temperature_K": None,
            "G_ads_eV": None,
            "method": "no_sign_change",
            "bracket_K": None,
            "iterations": 0,
            "search_min_K": tmin,
            "search_max_K": tmax,
            "G_ads_at_min_eV": float(values[0]),
            "G_ads_at_max_eV": float(values[-1]),
            "message": (
                "No G_ads=0 sign change was found in the configured "
                "temperature interval."
            ),
        }

    a, b = bracket
    fa = float(values[np.where(grid == a)[0][0]])
    fb = float(values[np.where(grid == b)[0][0]])

    # Bisection: robust and fully deterministic.
    # Fixed numerical tolerances for the one-dimensional bisection.
    temperature_tol = 1.0e-4
    g_tol = 1.0e-10
    max_iter = 100

    iterations = 0
    mid = 0.5 * (a + b)
    fm = None
    thermo_mid = None

    while iterations < max_iter:
        iterations += 1
        mid = 0.5 * (a + b)
        fm, thermo_mid = _g_ads(
            Eads, Au_freq, HgAu_freq, mid,
            config.pressure, config.include_zpe,
            config.gas_reference_pressure, config.thermo_method
        )

        if abs(fm) <= g_tol or (b - a) <= temperature_tol:
            break

        if fa * fm < 0.0:
            b = mid
            fb = fm
        else:
            a = mid
            fa = fm

    # Re-evaluate at the reported midpoint for a self-consistent result.
    g_final, thermo_final = _g_ads(
        Eads, Au_freq, HgAu_freq, mid,
        config.pressure, config.include_zpe,
        config.gas_reference_pressure, config.thermo_method
    )

    return {
        "found": True,
        "temperature_K": float(mid),
        "G_ads_eV": float(g_final),
        "method": "bisection",
        "bracket_K": [float(a), float(b)],
        "iterations": iterations,
        "temperature_tolerance_K": temperature_tol,
        "G_ads_tolerance_eV": g_tol,
        "thermodynamics": thermo_final,
    }


def run_workflow(config, calc):
    print("\n" + "=" * 70)
    print("Hg ON Au(111) MODULAR WORKFLOW")
    print("=" * 70)

    # ---------------------------------------------------------
    # Clean Au slab
    # ---------------------------------------------------------
    Au = build_au111(config)
    Au.calc = calc

    print("\nAu surface created")
    print_surface_details(Au, config.n_frozen_layers)

    Au = optimize_structure(
        Au,
        calc,
        name="Au111_opt",
        fmax=config.fmax_au,
        optimizer=config.optimizer,
        trajectory=config.save_trajectories,
    )

    E_Au = float(Au.get_potential_energy())

    # ---------------------------------------------------------
    # Hg/Au system
    # ---------------------------------------------------------
    site = get_adsorption_site_position(
        Au,
        config.adsorption_site,
        config.hg_height,
    )

    HgAu = add_hg_atom(
        Au,
        config,
        site,
    )
    HgAu.calc = calc

    print("\nHg/Au system created")
    print_surface_details(HgAu, config.n_frozen_layers)
    hg_surface_distance_before = get_hg_surface_distance(HgAu)
    print(
        f"Hg–Au surface distance BEFORE optimization: "
        f"{hg_surface_distance_before:.6f} Å"
    )

    HgAu = optimize_structure(
        HgAu,
        calc,
        name="HgAu_opt",
        fmax=config.fmax_hg,
        optimizer=config.optimizer,
        trajectory=config.save_trajectories,
    )

    E_HgAu = float(HgAu.get_potential_energy())

    hg_surface_distance_after = get_hg_surface_distance(HgAu)
    print(
        f"Hg–Au surface distance AFTER optimization: "
        f"{hg_surface_distance_after:.6f} Å"
    )
    print(
        f"Hg–Au surface distance change: "
        f"{hg_surface_distance_after - hg_surface_distance_before:+.6f} Å"
    )

    # ---------------------------------------------------------
    # Isolated Hg atom
    # ---------------------------------------------------------
    E_Hg = isolated_hg_energy(calc)

    Eads = adsorption_energy(
        E_HgAu,
        E_Au,
        E_Hg,
    )

    print_energy_summary(
        E_HgAu,
        E_Au,
        E_Hg,
        Eads,
    )

    # ---------------------------------------------------------
    # Phonons
    # ---------------------------------------------------------
    clear_cuda_cache()

    Au_ph = calculate_phonons(
        Au,
        calc,
        displacement=config.phonon_displacement,
        prefix="Au_phonons",
        verbose=config.verbose,
        frequency_cutoff_cm=config.phonon_cutoff_cm,
        symmetrize=config.phonon_symmetrize,
        n_frozen_layers_vibrations=config.n_frozen_layers_vibrations,
        hg_projection_threshold=config.hg_mode_projection_threshold,
        nfree=config.ase_vibration_nfree,
        method=config.phonon_method,
        fail_on_imaginary=config.fail_on_imaginary,
    )

    if config.cuda_clear_between_phonons:
        clear_cuda_cache()
        gc.collect()

    HgAu_ph = calculate_phonons(
        HgAu,
        calc,
        displacement=config.phonon_displacement,
        prefix="HgAu_phonons",
        verbose=config.verbose,
        frequency_cutoff_cm=config.phonon_cutoff_cm,
        symmetrize=config.phonon_symmetrize,
        n_frozen_layers_vibrations=config.n_frozen_layers_vibrations,
        hg_projection_threshold=config.hg_mode_projection_threshold,
        nfree=config.ase_vibration_nfree,
        method=config.phonon_method,
        fail_on_imaginary=config.fail_on_imaginary,
    )

    Au_freq = Au_ph["frequencies_cm"]
    HgAu_freq = HgAu_ph["frequencies_cm"]

    # ---------------------------------------------------------
    # Thermodynamics backend
    # ---------------------------------------------------------
    print(
        "Thermodynamics backend:",
        "project custom equations"
        if config.thermo_method == "custom"
        else "ASE HarmonicThermo + IdealGasThermo",
    )

    # ---------------------------------------------------------
    # Thermodynamics at configured temperature
    # ---------------------------------------------------------
    thermo = calculate_adsorption_thermodynamics(
        Eads,
        HgAu_freq,
        Au_freq,
        config.temperature,
        config.pressure,
        include_zpe=config.include_zpe,
        reference_pressure=config.gas_reference_pressure,
        method=config.thermo_method,
    )

    if config.thermo_method == "ase":
        print_ase_hgau_thermo_batch(
            HgAu_freq,
            config.temperatures,
        )

    print("\n" + "=" * 70)
    print("ADSORPTION THERMODYNAMICS")
    print("=" * 70)
    print(f"T = {config.temperature:.2f} K")
    print(f"p = {config.pressure:.6g} bar")
    print(f"ΔE_ads          = {thermo['E_ads']:.8f} eV")
    print(f"ΔZPE            = {thermo['ZPE_ads']:.8f} eV")
    print(f"ΔU_vib          = {thermo['delta_U_vib_thermal']:.8f} eV")
    print(f"H_gas thermal   = {thermo['H_gas_thermal']:.8f} eV")
    print(f"ΔH_ads          = {thermo['H_ads']:.8f} eV")
    print(f"S_gas           = {thermo['S_gas']:.8e} eV/K")
    print(f"ΔS_vib          = {thermo['delta_S_vib']:.8e} eV/K")
    print(f"ΔS_ads          = {thermo['S_ads']:.8e} eV/K")
    print(f"ΔF_vib          = {thermo['delta_F_vib']:.8f} eV")
    print(f"μ_gas thermal   = {thermo['mu_gas_thermal']:.8f} eV")
    print(f"ΔG_ads          = {thermo['G_ads']:.8f} eV")
    print("=" * 70)

    results = {
        "E_Au": E_Au,
        "E_HgAu": E_HgAu,
        "E_Hg": E_Hg,
        "E_ads": Eads,
        "Hg_surface_distance_initial_A": float(hg_surface_distance_before),
        "Hg_surface_distance_final_A": float(hg_surface_distance_after),
        "Hg_surface_distance_change_A": float(
            hg_surface_distance_after - hg_surface_distance_before
        ),
        **thermo,
        "vibrational_n_frozen_layers": int(config.n_frozen_layers_vibrations),
        "Au_n_mobile": int(len(Au_ph["phonon_atoms"])),
        "HgAu_n_mobile": int(len(HgAu_ph["phonon_atoms"])),
        "Au_n_modes": int(len(Au_freq)),
        "HgAu_n_modes": int(len(HgAu_freq)),
        "Hg_mode_frequencies_cm-1": (
            HgAu_ph["signed_frequencies_cm"][
                HgAu_ph["hg_mode_mask"]
            ].tolist()
        ),
    }

    # ---------------------------------------------------------
    # Temperature scan
    # ---------------------------------------------------------
    if config.calculate_temperature_scan or config.mode == "scan":
        scan = _thermo_scan(
            Eads,
            Au_freq,
            HgAu_freq,
            config,
        )
        results["temperature_scan"] = scan

    # ---------------------------------------------------------
    # Optional G_ads = 0 root
    # ---------------------------------------------------------
    if config.calculate_g_ads_zero_temperature:
        root = _find_g_ads_zero(
            Eads,
            Au_freq,
            HgAu_freq,
            config,
        )
        results["G_ads_zero_temperature"] = root

        print("\n" + "=" * 70)
        print("ΔG_ads = 0 TEMPERATURE")
        print("=" * 70)

        if root["found"]:
            print(
                f"T_cross        = {root['temperature_K']:.6f} K"
            )
            print(
                f"ΔG_ads(T_cross)= {root['G_ads_eV']:.6e} eV"
            )
            print(f"Method          = {root['method']}")
            print(
                "Final bracket   = "
                f"[{root['bracket_K'][0]:.8f}, "
                f"{root['bracket_K'][1]:.8f}] K"
            )
            print(f"Iterations      = {root['iterations']}")
        else:
            print("No ΔG_ads = 0 crossing found.")
            print(
                f"Search interval = "
                f"[{root['search_min_K']:.2f}, "
                f"{root['search_max_K']:.2f}] K"
            )
            print(
                f"ΔG_ads at limits = "
                f"{root['G_ads_at_min_eV']:.8f}, "
                f"{root['G_ads_at_max_eV']:.8f} eV"
            )
            print(root["message"])

    print("=" * 70)

    # Save only after all results have been assembled.
    _save_json(
        config.output_prefix + "_results.json",
        results,
    )

    print("\nWorkflow completed.")
    return results
