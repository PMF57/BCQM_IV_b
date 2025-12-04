# bcqm_cluster_toy/cluster_model.py

from __future__ import annotations

from typing import Tuple
import numpy as np


def simulate_cluster_trajectories(
    N: int,
    n_steps: int,
    n_ensembles: int,
    dt: float,
    seed: int,
    gamma: float = 0.01,
    sigma: float = 1.0,
) -> np.ndarray:
    """
    Simulate acceleration trajectories for a cluster of N probes.

    Returns
    -------
    a_com : np.ndarray
        Array of shape (n_ensembles, n_steps) containing the
        centre-of-mass acceleration time series for each ensemble.
    """
    rng = np.random.default_rng(seed + N)  # N-dependent offset for reproducibility

    # (n_ensembles, N, n_steps) would be largest; we instead loop over ensembles.
    a_com = np.empty((n_ensembles, n_steps), dtype=float)

    for e in range(n_ensembles):
        # a_i for each probe i=1..N
        a = np.zeros((N, n_steps), dtype=float)
        for n in range(n_steps - 1):
            noise = rng.normal(loc=0.0, scale=1.0, size=N)
            a[:, n + 1] = (1.0 - gamma) * a[:, n] + sigma * noise
        # centre-of-mass acceleration for this ensemble
        a_com[e, :] = a.mean(axis=0)

    return a_com