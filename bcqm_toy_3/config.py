# config.py
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List, Any, Dict

import yaml


@dataclass
class SimulationConfig:
    dt: float
    n_steps: int
    n_ensembles: int
    seed: Optional[int]
    sign_mode: int


@dataclass
class ScanConfig:
    wcoh_values: List[float]
    label: str


@dataclass
class OutputConfig:
    base_dir: Path


@dataclass
class Config:
    simulation: SimulationConfig
    scan: ScanConfig
    output: OutputConfig


def _expect(mapping: Dict[str, Any], key: str) -> Any:
    if key not in mapping:
        raise KeyError(f"Missing required key '{key}' in configuration section.")
    return mapping[key]


def load_config(path: Path) -> Config:
    """
    Load a YAML configuration file into a structured Config object.
    """
    with path.open("r", encoding="utf-8") as f:
        data = yaml.safe_load(f)

    sim_data = _expect(data, "simulation")
    scan_data = _expect(data, "scan")
    out_data = _expect(data, "output")

    sim = SimulationConfig(
        dt=float(_expect(sim_data, "dt")),
        n_steps=int(_expect(sim_data, "n_steps")),
        n_ensembles=int(_expect(sim_data, "n_ensembles")),
        seed=None if sim_data.get("seed") in (None, "null") else int(sim_data["seed"]),
        sign_mode=int(sim_data.get("sign_mode", 1)),
    )

    scan = ScanConfig(
        wcoh_values=[float(x) for x in _expect(scan_data, "wcoh_values")],
        label=str(scan_data.get("label", "scan")),
    )

    output = OutputConfig(
        base_dir=Path(out_data.get("base_dir", "outputs")),
    )

    return Config(simulation=sim, scan=scan, output=output)
