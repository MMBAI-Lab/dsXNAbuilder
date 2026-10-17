# Copyright (C) 2026 Pablo D. Dans
# SPDX-License-Identifier: GPL-3.0-or-later

import os
import sys

import pytest

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

THREEDNA = os.path.join(ROOT, "external", "3dna")


@pytest.fixture(scope="session")
def session():
    from dsxnabuilder.service import Session
    return Session()


@pytest.fixture(scope="session")
def lib(session):
    return session.fibers


@pytest.fixture(scope="session")
def catalog(session):
    return session.catalog


@pytest.fixture(scope="session")
def aliases(session):
    return session.aliases


@pytest.fixture(scope="session")
def raw_models():
    import json
    with open(os.path.join(ROOT, "data", "fiber_models.json"), encoding="utf-8") as fh:
        return json.load(fh)


needs_3dna = pytest.mark.skipif(not os.path.isdir(os.path.join(THREEDNA, "fiber")),
                                reason="3DNA data not fetched (tools/fetch_3dna.py)")
