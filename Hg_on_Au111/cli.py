#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Command-line interface."""

import argparse


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Hg adsorption thermodynamics on Au(111)"
    )

    parser.add_argument(
        "-c", "--config",
        default="hg_au111_config.ini",
        help="INI configuration file",
    )
    parser.add_argument(
        "--model",
        help="Override MACE model path",
    )
    parser.add_argument(
        "--mode",
        choices=["single", "scan"],
        help="Override calculation mode",
    )
    parser.add_argument(
        "-t", "--temperature",
        type=float,
        help="Override temperature in K",
    )
    parser.add_argument(
        "-p", "--pressure",
        type=float,
        help="Override Hg gas pressure in bar",
    )

    return parser.parse_args()
