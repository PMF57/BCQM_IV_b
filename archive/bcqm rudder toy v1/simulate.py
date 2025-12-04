"""
Simulation routines for bcqm_inertial_noise.

This module provides helpers to run individual trajectories and ensembles
of probe/cluster paths, and to convert positions into accelerations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np

from .config import Config
from .model import (
    GraphState,
    ProbeState,
    Params,
    initialise_graph,
    initialise_probe,
    make_params,
    current_position,
    step_probe,
    step_probe_bcqm,
)


@dataclass
class TrajectoryResult:
    positions: np.ndarray  # shape (n_steps,)
    accelerations: np.ndarray  # shape (n_steps-2,)


def compute_acceleration(positions: np.ndarray, dt: float) -> np.ndarray:
    """
    Compute a discrete acceleration signal from positions via a second
    central finite difference.

    For positions x_t at times t = n*dt, we use
    a_t = (x_{t+dt} - 2 x_t + x_{t-dt}) / dt^2.

    The returned array has length n_steps - 2.
    """
    positions = np.asarray(positions, dtype=float)
    if positions.ndim != 1:
        raise ValueError("positions must be a 1D array")
    if positions.size < 3:
        raise ValueError("need at least 3 time steps to compute acceleration")

    x_minus = positions[:-2]
    x_0 = positions[1:-1]
    x_plus = positions[2:]
    acc = (x_plus - 2.0 * x_0 + x_minus) / (dt ** 2)
    return acc


def simulate_trajectory(cfg: Config, w_coh: float, rng: np.random.Generator) -> TrajectoryResult:
    """
    Simulate a single probe trajectory for a given coherence horizon.

    Pseudocode:
    - construct graph state;
    - construct simulation parameters (including w_coh);
    - initialise probe state;
    - loop over n_steps, recording positions and stepping the probe;
    - compute acceleration from the positions.

    Returns positions and accelerations as 1D arrays.
    """
    graph: GraphState = initialise_graph(cfg)
    params: Params = make_params(cfg, w_coh)
    probe: ProbeState = initialise_probe(cfg, graph, rng)

    n_steps = cfg.simulation.n_steps
    positions = np.empty(n_steps, dtype=float)

    for n in range(n_steps):
        positions[n] = current_position(probe)
        probe = step_probe_bcqm(probe, graph, params, rng)      #or step_probe

    accelerations = compute_acceleration(positions, params.dt)
    return TrajectoryResult(positions=positions, accelerations=accelerations)


def simulate_ensemble(
    cfg: Config,
    w_coh: float,
    rng: np.random.Generator,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Simulate an ensemble of trajectories for a given coherence horizon.

    Parameters
    ----------
    cfg:
        Top-level configuration.
    w_coh:
        Coherence horizon for this run.
    rng:
        Base random number generator. Individual trajectories use
        derived seeds for reproducibility.

    Returns
    -------
    positions:
        Array of shape (n_trajectories, n_steps)
    accelerations:
        Array of shape (n_trajectories, n_steps-2)
    """
    n_traj = cfg.simulation.n_trajectories
    n_steps = cfg.simulation.n_steps

    positions = np.empty((n_traj, n_steps), dtype=float)
    accs = np.empty((n_traj, n_steps - 2), dtype=float)

    # Different base seed for each ensemble, derived from the passed RNG
    base_seed = int(rng.integers(0, 2**31 - 1))
    
    for k in range(n_traj):
        traj_rng = np.random.default_rng(base_seed + k)
        result = simulate_trajectory(cfg, w_coh, traj_rng)
        positions[k, :] = result.positions
        accs[k, :] = result.accelerations

    return positions, accs
