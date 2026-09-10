#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""MACE calculator setup with conservative CUDA cleanup."""

from __future__ import annotations

import gc
import platform

import ase
import mace
import torch
from mace.calculators import MACECalculator


def clear_cuda_cache() -> None:
    gc.collect()
    if torch.cuda.is_available():
        try:
            torch.cuda.empty_cache()
            torch.cuda.ipc_collect()
        except Exception:
            pass


def print_version() -> None:
    print("\n" + "=" * 70)
    print("VERSION INFORMATION")
    print("=" * 70)
    print("Python :", platform.python_version())
    print("ASE    :", ase.__version__)
    print("MACE   :", mace.__version__)
    print("Torch  :", torch.__version__)
    print("CUDA available:", torch.cuda.is_available())
    if torch.cuda.is_available():
        print("GPU:", torch.cuda.get_device_name(0))
    print("=" * 70)


def get_calculator(config):
    clear_cuda_cache()

    device = "cuda" if torch.cuda.is_available() else "cpu"

    print(f"Using device: {device}")
    print("MACE model:", config.model_path)

    calc = MACECalculator(
        model_paths=[config.model_path],
        device=device,
        default_dtype="float64",
    )

    clear_cuda_cache()
    return calc


def safe_energy(atoms):
    clear_cuda_cache()
    energy = float(atoms.get_potential_energy())
    clear_cuda_cache()
    return energy


def safe_forces(atoms):
    clear_cuda_cache()
    forces = atoms.get_forces()
    clear_cuda_cache()
    return forces
