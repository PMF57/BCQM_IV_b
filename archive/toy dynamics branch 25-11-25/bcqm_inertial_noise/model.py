"""
Model definitions for bcqm_inertial_noise.

This module defines the event-graph and probe/cluster state, together
with simple initialisation and update rules. The first version uses a
very simple 1D structure that is intended to be refined towards the
full BCQM event-graph picture.

The key design decision is that the probe state is always treated as a
vector in configuration space, even when there is only a single degree
of freedom. This makes it straightforward to extend to small clusters
later without changing the external interfaces.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np

from .config import Config


@dataclass
class GraphState:
    """
    Minimal graph state for a 1D chain of sites.

    Attributes
    ----------
    n_vertices:
        Number of sites in the chain.
    adjacency:
        Integer array of shape (n_vertices, 2) giving left/right
        neighbours (with simple reflecting boundary conditions).
    """

    n_vertices: int
    adjacency: np.ndarray  # shape (n_vertices, 2)


@dataclass
class ProbeState:
    """
    Probe or cluster configuration in a 1D effective coordinate.

    Attributes
    ----------
    positions:
        Array of shape (n_dof,) giving the current coordinate of each
        degree of freedom. For a single probe, n_dof = 1.
    vertex_indices:
        Integer array of shape (n_dof,) giving the current graph vertex
        index for each degree of freedom.
    last_step:
        The last step taken by the probe centre of mass (for a single
        probe this is just the previous Δx). Used to implement a simple
        persistence (memory) in the random walk.
    """

    positions: np.ndarray
    vertex_indices: np.ndarray
    last_step: float = 0.0
    
@dataclass
class Params:
    """
    Simulation parameters needed by the update rules.

    This mirrors the Config contents but keeps only what the low-level
    update functions actually require.
    """

    dt: float
    w_coh: float
    hop_radius: int
    drift_bias: float
    dx: float = 1.0  # spatial step; can be refined later


def initialise_graph(cfg: Config) -> GraphState:
    """
    Construct a simple 1D chain graph with reflecting boundaries.

    Parameters
    ----------
    cfg:
        Top-level configuration.

    Returns
    -------
    GraphState
        The initialised graph.
    """
    n = cfg.graph.n_vertices
    adjacency = np.zeros((n, 2), dtype=int)
    for i in range(n):
        left = i - 1 if i > 0 else 0
        right = i + 1 if i < n - 1 else n - 1
        adjacency[i, 0] = left
        adjacency[i, 1] = right
    return GraphState(n_vertices=n, adjacency=adjacency)


def initialise_probe(cfg: Config, graph: GraphState, rng: np.random.Generator) -> ProbeState:
    """
    Initialise a single-probe configuration near the centre of the graph.

    Parameters
    ----------
    cfg:
        Top-level configuration.
    graph:
        Graph state.
    rng:
        Random number generator.

    Returns
    -------
    ProbeState
        The initial probe state.
    """
    # For now, start at the centre with a small random offset.
    centre_idx = graph.n_vertices // 2
    x0 = float(centre_idx) * 1.0  # dx = 1 for now

    positions = np.array([x0], dtype=float)
    vertex_indices = np.array([centre_idx], dtype=int)

    return ProbeState(positions=positions, vertex_indices=vertex_indices, last_step=0.0)

def make_params(cfg: Config, w_coh: float) -> Params:
    """
    Extract the parameters needed by the low-level update rules.

    Parameters
    ----------
    cfg:
        Top-level configuration.
    w_coh:
        Coherence horizon for this run.

    Returns
    -------
    Params
    """
    return Params(
        dt=cfg.simulation.dt,
        w_coh=float(w_coh),
        hop_radius=cfg.graph.hop_radius,
        drift_bias=cfg.graph.drift_bias,
        dx=1.0,
    )


def current_position(probe: ProbeState) -> float:
    """
    Return the current centre-of-mass position of the probe/cluster.

    For a single probe this is just the single coordinate.
    """
    return float(np.mean(probe.positions))


def step_probe(
    probe: ProbeState,
    graph: GraphState,
    params: Params,
    rng: np.random.Generator,
) -> ProbeState:
    """
    Advance the probe/cluster state by one time step.

    Pseudocode/intended refinement:
    - determine a local hop kernel based on the current configuration,
      hop radius, and a phase-memory rule controlled by w_coh;
    - sample the next configuration from this kernel;
    - map configuration changes to effective position changes.

    The present minimal implementation uses a biased random-walk update
    for the centre-of-mass coordinate as a placeholder. This will be
    refined in later versions to reflect the full BCQM event-chain
    dynamics.

    Parameters
    ----------
    probe:
        Current probe/cluster state.
    graph:
        Graph state (currently used only for scale/reference).
    params:
        Simulation parameters.
    rng:
        Random number generator.

    Returns
    -------
    ProbeState
        Updated probe state.
    """
    x = current_position(probe)

    # Simple biased, persistent random step.
    bias = params.drift_bias
    p_right_base = 0.5 * (1.0 + bias)
    p_right_base = min(max(p_right_base, 0.0), 1.0)

    # Persistence parameter: rho ~ exp(-dt / W_coh), clipped to [0,1].
    # For large W_coh, rho ~ 1 -> strong persistence.
    # For small W_coh, rho ~ 0 -> nearly memoryless walk.
    rho = float(np.exp(-params.dt / params.w_coh))
    rho = min(max(rho, 0.0), 1.0)

    # Determine the probability of stepping right, taking into account
    # the previous step direction.
    if probe.last_step > 0:
        # Previous step was to the right: favour right again.
        p_right = rho + (1.0 - rho) * p_right_base
    elif probe.last_step < 0:
        # Previous step was to the left: disfavour right.
        p_right = (1.0 - rho) * p_right_base
    else:
        # No previous history: use the base bias.
        p_right = p_right_base

    p_right = min(max(p_right, 0.0), 1.0)

    # Decide step direction: +1 (right) or -1 (left) in vertex index.
    step_dir = 1 if rng.random() < p_right else -1

    # Current vertex index for the single probe
    current_idx = int(probe.vertex_indices[0])

    # Use graph adjacency: 0 = left neighbour, 1 = right neighbour
    if step_dir > 0:
        next_idx = int(graph.adjacency[current_idx, 1])
    else:
        next_idx = int(graph.adjacency[current_idx, 0])

    # Update vertex index and corresponding position
    new_vertex_indices = probe.vertex_indices.copy()
    new_vertex_indices[0] = next_idx

    new_positions = probe.positions.copy()
    new_positions[0] = float(next_idx) * params.dx

    # The effective last_step in coordinate space
    last_step = new_positions[0] - x

    return ProbeState(
        positions=new_positions,
        vertex_indices=new_vertex_indices,
        last_step=last_step,
    )
