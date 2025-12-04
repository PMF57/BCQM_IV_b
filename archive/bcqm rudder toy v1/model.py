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
    rudder_eps: float = 0.1


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

def compute_hop_probabilities_bcqm(
    probe: ProbeState,
    graph: GraphState,
    params: Params,
) -> tuple[float, float]:
    """
    BCQM IV_b rudder-only hop kernel on a 1D chain.

    State:
      - probe.last_step encodes the rudder as the last displacement
        along the chain; we use only its sign s_n ∈ {-1, 0, +1}.
    Parameters:
      - params.rudder_eps: small persistence parameter ε (0 < ε < 0.5)

    Rule (rudder-only part of the kernel):

        p_R^rud = 1/2 + ε * s_n
        p_L^rud = 1/2 - ε * s_n

    Interruptions with rate q = Δt / W_coh are handled in step_probe_bcqm;
    here we only return the rudder-biased probabilities.
    """
    # Rudder s_n ∈ {-1, 0, +1}. If no history, last_step should be 0.
    last_step = float(getattr(probe, "last_step", 0.0))
    if last_step > 0.0:
        s_n = 1.0
    elif last_step < 0.0:
        s_n = -1.0
    else:
        s_n = 0.0

    # Persistence parameter ε: how strongly the rudder tilts the bow.
    eps = float(getattr(params, "rudder_eps", 0.1))

    # Rudder-only probabilities
    p_R = 0.5 + eps * s_n
    p_L = 0.5 - eps * s_n

    # Clamp and normalise
    p_R = max(0.0, min(1.0, p_R))
    p_L = max(0.0, min(1.0, p_L))

    total = p_R + p_L
    if total <= 0.0:
        return 0.5, 0.5

    p_R /= total
    p_L /= total

    # p_left for left neighbour, p_right for right neighbour
    return p_L, p_R

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



def step_probe_bcqm(
    probe: ProbeState,
    graph: GraphState,
    params: Params,
    rng: np.random.Generator,
) -> ProbeState:
    """
    One BCQM IV_b step using the rudder + interruption hop kernel.

    We reuse the same movement / graph logic as step_probe, but the
    choice of left vs right is now driven by:

      - a rudder-only kernel (compute_hop_probabilities_bcqm), and
      - explicit interruptions at rate q = Δt / W_coh which override
        the rudder and produce an unbiased hop for that tick.
    """

    # Current position in coordinate space (for last_step update)
    x = current_position(probe)

    # Interruption probability: here we interpret W_coh in units of
    # hops, so q = 1 / W_coh is the probability that an interruption
    # occurs on a given step (clamped to [0, 1]).
    W_coh = float(params.w_coh)
    if W_coh <= 0.0:
        q = 1.0
    else:
        q = 1.0 / W_coh
        if q > 1.0:
            q = 1.0

    # Decide whether this tick is an interruption
    u_int = float(rng.random())
    if u_int < q:
        # Interruption: ignore rudder and use unbiased hop
        p_left, p_right = 0.5, 0.5
    else:
        # No interruption: use rudder-biased kernel
        p_left, p_right = compute_hop_probabilities_bcqm(probe, graph, params)

    # Choose a step direction based on these probabilities
    u = float(rng.random())
    if u < p_left:
        step_dir = -1
    else:
        step_dir = +1

    # Movement logic mirrors step_probe.

    # Current vertex index on the 1D chain
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
