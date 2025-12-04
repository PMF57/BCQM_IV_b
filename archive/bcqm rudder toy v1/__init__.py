"""
bcqm_inertial_noise
===================

Prototype code for the BCQM IV_b inertial-noise follow-up.

This package is deliberately small and modular. It provides:

- A minimal event-graph + probe/cluster model (model.py).
- Simulation routines for trajectories and ensembles (simulate.py).
- Spectral estimation and scaling tools (spectra.py).
- Configuration loading and validation (config.py).
- A simple command-line interface (cli.py).

The first working version focuses on a single probe degree of freedom,
implemented in a way that generalises naturally to small clusters.
"""
