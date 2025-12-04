# bcqm_cluster_toy/config.py

from dataclasses import dataclass
from pathlib import Path
from typing import List
import yaml


@dataclass
class SimulationConfig:
    dt: float
    n_steps: int
    n_ensembles: int
    seed: int


@dataclass
class ClusterScanConfig:
    N_values: List[int]
    label: str


@dataclass
class OutputConfig:
    base_dir: str


@dataclass
class Config:
    simulation: SimulationConfig
    cluster: ClusterScanConfig
    output: OutputConfig


def load_config(path: str | Path) -> Config:
    """Load YAML config into a strongly-typed Config object."""
    path = Path(path)
    with path.open("r") as f:
        data = yaml.safe_load(f)

    sim = data["simulation"]
    cl = data["cluster"]
    out = data["output"]

    sim_cfg = SimulationConfig(
        dt=float(sim["dt"]),
        n_steps=int(sim["n_steps"]),
        n_ensembles=int(sim["n_ensembles"]),
        seed=int(sim["seed"]),
    )

    cluster_cfg = ClusterScanConfig(
        N_values=[int(x) for x in cl["N_values"]],
        label=str(cl["label"]),
    )

    out_cfg = OutputConfig(base_dir=str(out["base_dir"]))

    return Config(simulation=sim_cfg, cluster=cluster_cfg, output=out_cfg)