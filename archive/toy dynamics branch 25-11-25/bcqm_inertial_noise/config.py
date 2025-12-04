"""
Configuration handling for bcqm_inertial_noise.

The preferred format is YAML. A minimal configuration might look like:

output_dir: outputs/wcoh_scan_v1
random_seed: 12345

simulation:
  n_trajectories: 256
  n_steps: 8192
  dt: 1.0
  w_coh_list: [50.0, 100.0]

graph:
  n_vertices: 512
  hop_radius: 3
  initial_condition: "single_site"
  drift_bias: 0.1

analysis:
  window: "hann"
  segment_length: 4096
  overlap: 0.5
  reference_frequency_index: 10
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Any, Dict

try:
    import yaml  # type: ignore
except ImportError:  # pragma: no cover - yaml may not be present in all environments
    yaml = None  # type: ignore


@dataclass
class SimulationConfig:
    n_trajectories: int
    n_steps: int
    dt: float
    w_coh_list: Sequence[float]


@dataclass
class GraphConfig:
    n_vertices: int
    hop_radius: int
    initial_condition: str = "single_site"
    drift_bias: float = 0.0
    edge_deletion_prob: float = 0.0  # for universality/robustness tests


@dataclass
class AnalysisConfig:
    window: str = "hann"
    segment_length: int = 4096
    overlap: float = 0.5
    reference_frequency_index: int = 10


@dataclass
class StorageConfig:
    save_full_trajectories: bool = False
    save_summary_stats: bool = True
    format: str = "npz"  # or "csv"


@dataclass
class Config:
    output_dir: Path
    random_seed: int
    simulation: SimulationConfig
    graph: GraphConfig
    analysis: AnalysisConfig
    storage: StorageConfig


def _validate_positive_int(name: str, value: int) -> None:
    if value <= 0:
        raise ValueError(f"{name} must be positive, got {value!r}")


def _validate_positive_float(name: str, value: float) -> None:
    if value <= 0.0:
        raise ValueError(f"{name} must be positive, got {value!r}")


def load_config(path: str | Path) -> Config:
    """
    Load a configuration file and perform basic validation.

    Parameters
    ----------
    path:
        Path to a YAML configuration file.

    Returns
    -------
    Config
        Parsed and validated configuration dataclass.
    """
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Configuration file not found: {p}")

    if yaml is None:
        raise RuntimeError(
            "PyYAML is required to load configuration files, but it is not installed. "
            "Install it with `pip install pyyaml`."
        )

    with p.open("r", encoding="utf8") as f:
        raw: Dict[str, Any] = yaml.safe_load(f)

    # Top-level fields
    output_dir = Path(raw.get("output_dir", "outputs"))
    random_seed = int(raw.get("random_seed", 0))

    sim_raw = raw.get("simulation", {})
    graph_raw = raw.get("graph", {})
    analysis_raw = raw.get("analysis", {})
    storage_raw = raw.get("storage", {})

    sim = SimulationConfig(
        n_trajectories=int(sim_raw.get("n_trajectories", 16)),
        n_steps=int(sim_raw.get("n_steps", 2048)),
        dt=float(sim_raw.get("dt", 1.0)),
        w_coh_list=list(sim_raw.get("w_coh_list", [100.0])),
    )

    graph = GraphConfig(
        n_vertices=int(graph_raw.get("n_vertices", 256)),
        hop_radius=int(graph_raw.get("hop_radius", 3)),
        initial_condition=str(graph_raw.get("initial_condition", "single_site")),
        drift_bias=float(graph_raw.get("drift_bias", 0.0)),
        edge_deletion_prob=float(graph_raw.get("edge_deletion_prob", 0.0)),
    )

    analysis = AnalysisConfig(
        window=str(analysis_raw.get("window", "hann")),
        segment_length=int(analysis_raw.get("segment_length", 1024)),
        overlap=float(analysis_raw.get("overlap", 0.5)),
        reference_frequency_index=int(analysis_raw.get("reference_frequency_index", 10)),
    )

    storage = StorageConfig(
        save_full_trajectories=bool(storage_raw.get("save_full_trajectories", False)),
        save_summary_stats=bool(storage_raw.get("save_summary_stats", True)),
        format=str(storage_raw.get("format", "npz")),
    )

    # Basic validation
    _validate_positive_int("n_trajectories", sim.n_trajectories)
    _validate_positive_int("n_steps", sim.n_steps)
    _validate_positive_float("dt", sim.dt)
    for w in sim.w_coh_list:
        _validate_positive_float("W_coh", float(w))

    _validate_positive_int("n_vertices", graph.n_vertices)
    if not (1 <= graph.hop_radius <= graph.n_vertices - 1):
        raise ValueError(
            f"hop_radius must be between 1 and n_vertices-1, got {graph.hop_radius}"
        )
    if not (0.0 <= graph.edge_deletion_prob <= 1.0):
        raise ValueError(
            f"edge_deletion_prob must be between 0 and 1, got {graph.edge_deletion_prob}"
        )

    _validate_positive_int("segment_length", analysis.segment_length)
    if not (0.0 < analysis.overlap < 1.0):
        raise ValueError(
            f"analysis.overlap must be in (0,1), got {analysis.overlap}"
        )

    output_dir.mkdir(parents=True, exist_ok=True)

    return Config(
        output_dir=output_dir,
        random_seed=random_seed,
        simulation=sim,
        graph=graph,
        analysis=analysis,
        storage=storage,
    )
