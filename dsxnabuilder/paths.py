# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

"""
dsxnabuilder.paths
==================

Where the package finds its data. The data files live next to the package in
the repository (``data/``) and the prebuilt modXNA residues in ``web/lib/``,
which the static web app serves as well; the package is therefore used from
a checkout (``pip install -e .``).
"""

import os

ROOT = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
DATA = os.path.join(ROOT, "data")
RESIDUE_DIR = os.path.join(ROOT, "web", "lib", "residues")
FRCMOD = os.path.join(ROOT, "web", "lib", "frcmod.modxna")
