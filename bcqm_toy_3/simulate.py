# simulate.py
from __future__ import annotations

from pathlib import Path
from typing import Tuple, List

import numpy as np

from .config import SimulationConfig
from .model import simulate_acceleration_series, KernelParameters
from .spectra import estimate_Sa, estimate_amplitude_and_omega_c


def run_ensemble_for_wcoh(
    wcoh: float,
    sim_cfg: SimulationConfig,
    rng: np.random.Generator,
    out_dir: Path,
) -> Tuple[float, float]:
    """
    Run an ensemble of trajectories for a given W_coh, average the PSDs,
    save the mean spectrum, and return (A, omega_c).
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    Sa_list: List[np.ndarray] = []
    freqs: np.ndarray | None = None

    kernel_params = KernelParameters()  # can be made configurable later

    for _ in range(sim_cfg.n_ensembles):
        acc = simulate_acceleration_series(
            n_steps=sim_cfg.n_steps,
            dt=sim_cfg.dt,
            sign_mode=sim_cfg.sign_mode,
            rng=rng,
            params=kernel_params,
        )
        f, Sa = estimate_Sa(acc, sim_cfg.dt)
        freqs = f
        Sa_list.append(Sa)

    Sa_mean = np.mean(np.stack(Sa_list, axis=0), axis=0)

    # Save spectrum
    np.savez(
        out_dir / f"Wcoh_{wcoh:g}.npz",
        freqs=freqs,
        Sa=Sa_mean,
    )

    # Extract amplitude and characteristic frequency
    A, omega_c = estimate_amplitude_and_omega_c(freqs, Sa_mean)
    return A, omega_c
