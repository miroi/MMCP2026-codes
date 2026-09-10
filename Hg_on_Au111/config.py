#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Configuration management for Hg adsorption on Au(111)."""

from __future__ import annotations

import configparser
import os


class Config:
    def __init__(self):
        self.vacuum_size = 18.0
        self.surface_size = (4, 4, 6)
        self.lattice_constant = 4.08
        self.n_frozen_layers = 2
        # Total number of bottom Au layers frozen only for vibrational/thermochemical calculations.
        # This does not affect geometry relaxation and may equal the full slab thickness.
        self.n_frozen_layers_vibrations = 2

        self.adsorption_site = "fcc"
        self.hg_height = 3.0

        self.fmax_au = 0.01
        self.fmax_hg = 0.005
        self.optimizer = "FIRE"

        self.phonon_displacement = 0.005
        self.phonon_method = "custom"
        self.phonon_cutoff_cm = 1.0
        self.phonon_symmetrize = True
        self.ase_vibration_nfree = 2
        self.fail_on_imaginary = True

        self.use_hg_projected_modes = True
        self.hg_mode_projection_threshold = 0.20

        self.temperature = 298.15
        self.temperatures = [200, 250, 298.15, 350, 400, 500, 600]
        self.pressure = 1.0
        self.gas_reference_pressure = 1.0

        self.include_zpe = True
        self.thermo_method = "custom"
        self.calculate_temperature_scan = False
        self.calculate_g_ads_zero_temperature = False

        self.mode = "single"
        self.output_prefix = "hg_au111"
        self.verbose = True
        self.save_trajectories = True
        self.cuda_clear_between_phonons = True

        self.model_path = None
        self.config_file = "hg_au111_config.ini"

    @staticmethod
    def _as_bool(value: str, default: bool) -> bool:
        if value is None:
            return default
        return value.strip().lower() in {"true", "yes", "1", "on"}

    def load(self, filename: str) -> bool:
        self.config_file = filename

        if not os.path.exists(filename):
            print(f"WARNING: {filename} not found. Using defaults.")
            return False

        cfg = configparser.ConfigParser()
        cfg.read(filename)

        if "system" in cfg:
            s = cfg["system"]
            self.vacuum_size = s.getfloat("vacuum_size", self.vacuum_size)
            self.surface_size = tuple(
                int(x.strip())
                for x in s.get("surface_size", "4,4,6").split(",")
            )
            if len(self.surface_size) != 3:
                raise ValueError("surface_size must contain exactly 3 integers")

            self.lattice_constant = s.getfloat(
                "lattice_constant", self.lattice_constant
            )
            self.n_frozen_layers = s.getint(
                "n_frozen_layers", self.n_frozen_layers
            )
            self.n_frozen_layers_vibrations = s.getint(
                "n_frozen_layers_vibrations", self.n_frozen_layers_vibrations
            )
            self.adsorption_site = s.get(
                "adsorption_site", self.adsorption_site
            ).strip().lower()
            self.hg_height = s.getfloat("hg_height", self.hg_height)

            self.fmax_au = s.getfloat("fmax_au", self.fmax_au)
            self.fmax_hg = s.getfloat("fmax_hg", self.fmax_hg)
            self.phonon_displacement = s.getfloat(
                "phonon_displacement", self.phonon_displacement
            )

        if "calculation" in cfg:
            c = cfg["calculation"]
            self.mode = c.get("mode", self.mode).strip().lower()
            self.temperature = c.getfloat(
                "single_temperature", self.temperature
            )
            self.calculate_g_ads_zero_temperature = self._as_bool(
                c.get("calculate_g_ads_zero_temperature"),
                self.calculate_g_ads_zero_temperature,
            )

            temps = c.get("temperatures", "").strip()
            if temps:
                self.temperatures = [
                    float(x.strip()) for x in temps.split(",") if x.strip()
                ]

            self.pressure = c.getfloat("pressure", self.pressure)
            self.output_prefix = c.get(
                "output_prefix", self.output_prefix
            )
            self.verbose = self._as_bool(
                c.get("verbose"), self.verbose
            )
            self.save_trajectories = self._as_bool(
                c.get("save_trajectories"), self.save_trajectories
            )
            self.optimizer = c.get("optimizer", self.optimizer).strip()

        if "phonons" in cfg:
            p = cfg["phonons"]
            self.phonon_method = p.get("method", self.phonon_method).strip().lower()
            self.phonon_cutoff_cm = p.getfloat(
                "frequency_cutoff_cm", self.phonon_cutoff_cm
            )
            self.phonon_symmetrize = self._as_bool(
                p.get("symmetrize"), self.phonon_symmetrize
            )
            self.fail_on_imaginary = self._as_bool(
                p.get("fail_on_imaginary"), self.fail_on_imaginary
            )
            self.use_hg_projected_modes = self._as_bool(
                p.get("use_hg_projected_modes"),
                self.use_hg_projected_modes,
            )
            self.hg_mode_projection_threshold = p.getfloat(
                "hg_mode_projection_threshold",
                self.hg_mode_projection_threshold,
            )

        if "thermodynamics" in cfg:
            t = cfg["thermodynamics"]
            self.thermo_method = t.get("method", self.thermo_method).strip().lower()
            self.include_zpe = self._as_bool(
                t.get("include_zpe"), self.include_zpe
            )
            self.gas_reference_pressure = t.getfloat(
                "gas_reference_pressure",
                self.gas_reference_pressure,
            )
            self.calculate_temperature_scan = self._as_bool(
                t.get("calculate_temperature_scan"),
                self.calculate_temperature_scan,
            )

        if "runtime" in cfg:
            r = cfg["runtime"]
            self.cuda_clear_between_phonons = self._as_bool(
                r.get("cuda_clear_between_phonons"),
                self.cuda_clear_between_phonons,
            )

        if "model" in cfg:
            self.model_path = cfg["model"].get(
                "path", self.model_path
            ).strip()

        return True

    def validate(self) -> bool:
        if len(self.surface_size) != 3:
            raise ValueError("surface_size must have three entries")

        if self.n_frozen_layers < 0:
            raise ValueError("n_frozen_layers must be >= 0")

        n_layers = int(self.surface_size[2])
        if self.n_frozen_layers >= n_layers:
            raise ValueError(
                "n_frozen_layers must be smaller than the number of layers "
                "because the slab must retain mobile atoms during geometry relaxation"
            )
        if self.n_frozen_layers_vibrations < self.n_frozen_layers:
            raise ValueError(
                "n_frozen_layers_vibrations must be >= n_frozen_layers"
            )
        if self.n_frozen_layers_vibrations > n_layers:
            raise ValueError(
                "n_frozen_layers_vibrations must be <= the number of Au layers"
            )

        if self.phonon_method not in {"custom", "ase"}:
            raise ValueError("phonons method must be either custom or ase")

        if self.phonon_displacement <= 0:
            raise ValueError("phonon_displacement must be > 0")

        if self.phonon_cutoff_cm < 0:
            raise ValueError("phonon frequency cutoff must be >= 0")

        if not (0.0 < self.hg_mode_projection_threshold <= 1.0):
            raise ValueError(
                "hg_mode_projection_threshold must be in (0, 1]"
            )

        if self.thermo_method not in {"custom", "ase"}:
            raise ValueError("thermodynamics method must be either 'custom' or 'ase'")

        if self.temperature <= 0:
            raise ValueError("temperature must be > 0")

        if self.pressure <= 0:
            raise ValueError("pressure must be > 0")

        if self.gas_reference_pressure <= 0:
            raise ValueError("gas_reference_pressure must be > 0")

        if not self.model_path:
            raise RuntimeError("No MACE model specified")

        if not os.path.exists(self.model_path):
            raise FileNotFoundError(self.model_path)

        return True
