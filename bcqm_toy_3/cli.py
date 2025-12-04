# cli.py
from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from .config import load_config
from .simulate import run_ensemble_for_wcoh
from .spectra import fit_beta


def run_from_config(config_path: Path) -> None:
    cfg = load_config(config_path)

    out_dir = cfg.output.base_dir / cfg.scan.label
    out_dir.mkdir(parents=True, exist_ok=True)

    # RNG seeded once for reproducibility
    rng = np.random.default_rng(cfg.simulation.seed)

    w_values = []
    A_values = []
    omega_values = []

    for wcoh in cfg.scan.wcoh_values:
        print(f"Running ensemble for W_coh = {wcoh} ...")
        A, omega_c = run_ensemble_for_wcoh(
            wcoh=float(wcoh),
            sim_cfg=cfg.simulation,
            rng=rng,
            out_dir=out_dir,
        )
        print(f"  A = {A:.6g}, omega_c = {omega_c:.6g}")
        w_values.append(float(wcoh))
        A_values.append(float(A))
        omega_values.append(float(omega_c))

    # Write amplitude summary CSV
    csv_path = out_dir / "amplitude_scaling.csv"
    with csv_path.open("w", encoding="utf-8") as f:
        f.write("Wcoh,A,omega_c\n")
        for w, A, oc in zip(w_values, A_values, omega_values):
            f.write(f"{w},{A},{oc}\n")

    # Fit beta from the A(W_coh) values
    beta, beta_err = fit_beta(np.array(w_values), np.array(A_values))

    print()
    print("Wcoh:", w_values)
    print("A:", A_values)
    print(f"Fitted beta: {beta}")
    print(f"Estimated error: {beta_err}")


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="bcqm_toy_3",
        description="BCQM toy inertial-noise model (direction-free, with beta fit).",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser(
        "run", help="Run a W_coh scan from a YAML configuration file."
    )
    run_parser.add_argument("config", type=str, help="Path to YAML config file.")

    args = parser.parse_args()

    if args.command == "run":
        run_from_config(Path(args.config))


if __name__ == "__main__":
    main()
