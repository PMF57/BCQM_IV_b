# model.py
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


@dataclass
class KernelParameters:
    """
    Parameters controlling the toy acceleration kernel.

    This is deliberately simple: an Ornstein-Uhlenbeck-like update
    for the acceleration with a fixed relaxation time and noise strength.
    """
    relaxation_time: float = 10.0
    noise_strength: float = 1.0


def simulate_acceleration_series(
    n_steps: int,
    dt: float,
    sign_mode: int,
    rng: np.random.Generator,
    params: Optional[KernelParameters] = None,
) -> np.ndarray:
    """
    Simulate a single acceleration time series.

    sign_mode = +1 or -1 parameterises the "sign choice" in the toy dynamics
    (the change we previously tested by flipping + to - in the update rule).
    """
    if params is None:
        params = KernelParameters()

    series = np.empty(n_steps, dtype=float)
    a = 0.0

    # Simple AR(1) / OU-like dynamics for acceleration
    gamma = dt / params.relaxation_time
    decay = max(0.0, 1.0 - gamma)

    for i in range(n_steps):
        xi = rng.normal()
        # The place where the sign variant enters:
        a = decay * a + float(sign_mode) * params.noise_strength * xi
        series[i] = a

    return series
