# bcqm_cluster_toy/cluster_simulate.py

from __future__ import annotations

from pathlib import Path
import csv

from .config import load_config
from .cluster_model import simulate_cluster_trajectories
from .cluster_spectra import ensemble_com_spectrum


def run_cluster_scan(config_path: str | Path) -> None:
    """
    Run an N-scan for the cluster toy and save COM spectra + amplitudes.

    For each N in config.cluster.N_values, this writes:
      - cluster_N{N}.npz : omega, Sa, N, dt
    and a CSV:
      - amplitude_scaling_COM.csv : N, A_com, omega_c_com
    into base_dir / label.
    """
    cfg = load_config(config_path)
    sim = cfg.simulation
    cl = cfg.cluster
    out = cfg.output

    base_dir = Path(out.base_dir).expanduser().resolve() / cl.label
    base_dir.mkdir(parents=True, exist_ok=True)

    rows = []

    for N in cl.N_values:
        print(f"Running cluster simulation for N = {N} ...")
        a_com = simulate_cluster_trajectories(
            N=N,
            n_steps=sim.n_steps,
            n_ensembles=sim.n_ensembles,
            dt=sim.dt,
            seed=sim.seed,
        )

        omega, Sa_avg, A_com, omega_c_com = ensemble_com_spectrum(a_com, sim.dt)

        npz_path = base_dir / f"cluster_N{N}.npz"
        import numpy as np  # local import to keep top of file light

        np.savez(
            npz_path,
            omega=omega,
            Sa=Sa_avg,
            N=N,
            dt=sim.dt,
        )

        rows.append(
            {
                "N": N,
                "A_com": A_com,
                "omega_c_com": omega_c_com,
            }
        )

        print(f"  -> A_com = {A_com:.6g}, omega_c_com = {omega_c_com:.6g}")

    # Write amplitude scaling CSV
    csv_path = base_dir / "amplitude_scaling_COM.csv"
    with csv_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["N", "A_com", "omega_c_com"])
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    print(f"\nCompleted cluster N-scan. Results in: {base_dir}")
    print(f"  - Spectra: cluster_N*.npz")
    print(f"  - Amplitudes: {csv_path}")