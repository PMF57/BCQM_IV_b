## Code branches and toy models

We used two main code branches to validate the `bcqm_inertial_noise`
pipeline (trajectory → spectrum → amplitude → scaling):

### 1. `toy_dynamics_branch_21-11-25`

This branch contains two simple 1D random-walk models:

- **Explicitly W-coh–scaled toy**

  The step statistics are chosen so that the inertial-noise amplitude
  must scale as a power of the coherence horizon,
  e.g. via an explicit scaling of the step size with `W_coh`. Running
  the full pipeline gives a fitted exponent

  \[
    A(W_{\mathrm{coh}}) \propto W_{\mathrm{coh}}^{-\beta},
    \qquad \beta \approx 0.7 .
  \]

  This is a positive control: it shows the code can recover
  a genuine `W_coh`-dependence when it is built into the dynamics.

- **Persistent random walk**

  A nearest-neighbour random walk with a persistence parameter
  `rho = exp(-dt / W_coh)`, so `W_coh` only affects the short-range
  step correlations. In this case the measured amplitude
  `A(W_coh)` is essentially flat over the scan and the fitted exponent
  is

  \[
    \beta \approx 0 ,
  \]

  consistent with no `W_coh` scaling. This is a negative control:
  the pipeline does not artificially invent a power law when only weak
  local correlations are present.

### 2. `bcqm_rudder_toy_v1`

This branch implements a BCQM-motivated **rudder + interruption**
kernel on a 1D chain:

- The probe carries a `last_step` field; its sign defines a local
  rudder `s_n ∈ {−1,0,+1}` (“I was there, now I am here”).
- In the absence of interruptions, hop probabilities are tilted by the
  rudder:

  \[
    p_R = 1/2 + \varepsilon s_n, \quad
    p_L = 1/2 - \varepsilon s_n,
  \]

  with a small persistence parameter `epsilon`.

- On each tick, an **interruption** occurs with probability
  `q = 1 / W_coh`. On an interrupted tick the rudder is ignored and an
  unbiased hop is taken (`p_L = p_R = 0.5`), representing the
  coarse-grained influence of other threads/events.

For `W_coh` scanned over `[5, 160]` (in hops) the measured amplitude
`A(W_coh)` is approximately constant
(`~ 7×10⁻⁴` with small variations), and a log–log fit gives

\[
  \beta \approx -0.02 ,
\]

again consistent with no significant `W_coh` scaling.

Taken together with the toys in `toy_dynamics_branch_21-11-25`, this
confirms:

- the pipeline correctly recovers a non-trivial `W_coh` scaling
  when present, and
- a finite coherence horizon, by itself, does not guarantee a strong
  inertial-noise power law; the microdynamics of how `W_coh` enters
  the hop statistics is crucial.