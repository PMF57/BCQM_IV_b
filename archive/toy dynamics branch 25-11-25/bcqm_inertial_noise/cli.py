"""
Command-line interface for bcqm_inertial_noise.

Provides a simple "run" command that takes a configuration file, runs
simulations for each W_coh listed, estimates spectra, and writes results
to the configured output directory.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List

import numpy as np

from .config import load_config
from .simulate import simulate_ensemble
from .spectra import estimate_Sa, estimate_amplitude


def main(argv: List[str] | None = None) -> None:
    parser = argparse.ArgumentParser(prog="bcqm_inertial_noise")
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="run simulations from a configuration file")
    run_parser.add_argument("config", type=str, help="path to YAML configuration file")

    args = parser.parse_args(argv)

    if args.command == "run":
        run_from_config(Path(args.config))
    else:  # pragma: no cover - should not happen with required=True
        parser.error(f"Unknown command {args.command!r}")


def run_from_config(config_path: Path) -> None:
    cfg = load_config(config_path)
    rng = np.random.default_rng(cfg.random_seed)

    results = []

    for w_coh in cfg.simulation.w_coh_list:
        print(f"Running ensemble for W_coh = {w_coh} ...")
        positions, accs = simulate_ensemble(cfg, w_coh, rng)
        # Estimate spectrum using ensemble-averaged Welch method
        from .spectra import SpectrumResult, estimate_Sa  # local import to avoid cycles
        spec = estimate_Sa(
            accs,
            dt=cfg.simulation.dt,
            segment_length=cfg.analysis.segment_length,
            overlap=cfg.analysis.overlap,
            window=cfg.analysis.window,
        )
        A = estimate_amplitude(spec, cfg.analysis.reference_frequency_index)

        results.append((float(w_coh), spec, A))

        # Save per-W_coh results
        out_base = cfg.output_dir / f"Wcoh_{w_coh:g}"
        out_base.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            out_base.with_suffix(".npz"),
            omega=spec.omega,
            Sa=spec.Sa,
            Wcoh=float(w_coh),
            A=A,
        )
        print(f"Saved spectrum to {out_base.with_suffix('.npz')}")

    # Optionally, save a simple summary table
    if cfg.storage.save_summary_stats:
        summary_path = cfg.output_dir / "amplitude_scaling.csv"
        with summary_path.open("w", encoding="utf8") as f:
            # Header
            f.write("Wcoh,A\n")
            # One line per coherence horizon
            for w, spec, A in results:
                f.write(f"{w},{A}\n")
        print(f"Wrote amplitude summary to {summary_path}")

if __name__ == "__main__":  # pragma: no cover
    main()
