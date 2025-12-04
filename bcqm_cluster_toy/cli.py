# bcqm_cluster_toy/cli.py

from __future__ import annotations

import argparse
from pathlib import Path

from .cluster_simulate import run_cluster_scan


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="Cluster toy N-scan for BCQM IV_b."
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_p = subparsers.add_parser(
        "run",
        help="Run an N-scan over cluster sizes for the toy model.",
    )
    run_p.add_argument(
        "config",
        type=str,
        help="Path to YAML config file (e.g. bcqm_cluster_toy/configs/cluster_n_scan.yml)",
    )

    args = parser.parse_args(argv)

    if args.command == "run":
        config_path = Path(args.config)
        run_cluster_scan(config_path)
    else:
        parser.error(f"Unknown command: {args.command!r}")


if __name__ == "__main__":
    main()