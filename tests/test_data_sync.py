"""The runtime JSON tables are in sync with their YAML sources."""

import json
import os

import pytest

yaml = pytest.importorskip("yaml")
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))


@pytest.mark.parametrize("name", ["aliases", "base_pairing"])
def test_yaml_json_in_sync(name):
    with open(os.path.join(ROOT, "data", f"{name}.yaml"), encoding="utf-8") as fh:
        y = yaml.safe_load(fh)
    with open(os.path.join(ROOT, "data", f"{name}.json"), encoding="utf-8") as fh:
        j = json.load(fh)
    assert y == j, "run: python tools/build_data.py"
