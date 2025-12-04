# BCQM IV_b – Inertial Noise Toy Models: TESTING.md

## 1. Purpose

This file documents how we test the **inertial–noise pipeline** for BCQM IV_b:

- starting from a simple, thread-level toy kernel,
- generating acceleration time series,
- computing power spectra,
- extracting an effective **noise amplitude** \(A(W_{	ext{coh}})\),
- and fitting a scaling exponent \(eta\) in
  \[
  A(W_{	ext{coh}}) \sim C \; W_{	ext{coh}}^{-eta}.
  \]
- and, in a separate cluster branch, how the **centre-of-mass (COM) noise amplitude**
  \(A_{\mathrm{COM}}(N)\) scales with the number of probes \(N\).

At this stage **we are *not* trying to build the final physical kernel**. The current code provides:

- archived examples that deliberately *do* or *do not* scale with \(W_{	ext{coh}}\),  
- a **canonical control toy model** (`bcqm_toy_3`) that is **blind** to \(W_{	ext{coh}}\) and therefore should give **\(eta pprox 0\)**, and  
- an **independent-probe cluster toy** (`bcqm_cluster_toy`) whose COM amplitude should fall as \(A_{\mathrm{COM}}(N) \propto N^{-1/2}\) if the probes really are uncorrelated.

This validates the numerical pipeline itself:
simulation → spectrum → amplitude \& \(\omega_c\) → log–log fit for \(eta\) and, in the cluster case, log–log fits of \(A_{\mathrm{COM}}(N)\) vs.~\(N\).


---

## 2. Environment / prerequisites

You need a working Python 3 environment with (at minimum):

- `numpy`
- `pyyaml`
- `scipy` (if any later fits use it; current pipeline can run with `numpy` + `pyyaml` only)
- `matplotlib` (optional, for your own plotting)

On macOS, something like:

```bash
python3 -m venv venv
source venv/bin/activate
pip install numpy pyyaml matplotlib
```

is sufficient for the current toy.

---

## 3. Repository / folder layout

For BCQM IV_b tests we use a folder of the form:

```text
bcqm_toy_3/
    __init__.py
    cli.py
    config.py
    model.py
    spectra.py
    configs/
        wcoh_scan_phase1.yml
    outputs/
        (created automatically on first run)
```

You can keep this `bcqm_toy_3` folder on your Desktop (or anywhere convenient). The `outputs/` directory is created automatically when you run the scan.

> **Note**  
> Earlier experimental toys (Toy A and Toy B below) were kept in separate folders / branches (e.g. archived `toy_dynamics_21-11-25`, `bcqm_rudder_toy_v1_26-11-25`). They are not required to run the current canonical control model.

---

## 4. Toy models overview

### 4.1 Toy model A – Explicit \(W_{\text{coh}}\) scaling (archived)

**Status:** Archived; kept only as a pipeline sanity check.

- Implemented a deliberately **\(W_{\text{coh}}\)-dependent** kernel.
- Example behaviours:
  - A model with **explicit amplitude scaling** \(A \propto W_{\text{coh}}^{-1/2}\) that yielded a fitted \(\beta \approx 0.7\).
  - A persistent random walk with correlation parameter \(\rho = \exp(-\Delta t / W_{\text{coh}})\) that gave \(\beta \approx 0\).
- Purpose:
  - Prove that the pipeline can:
    - recover a non-zero \(\beta\) when it is put in by hand, and
    - return \(\beta \approx 0\) for a model with no clean amplitude scaling.

These models are useful as **historical tests**, but the notation and language pre-date the final IV_b primitives clean-up (events / edges / threads / kernels).

### 4.2 Toy model B – “rudder / bow” exploratory toy (archived)

**Status:** Archived; not used in the final narrative.

- Used more geometric language:
  - “rudder”, “bow”, “left/right”, etc., as if there were already a background space.
- Useful **exploration** at the time, but:
  - conceptually inconsistent with the **graph-first, pre-spacetime** minimalism we have now locked in for BCQM,
  - therefore excluded from the canonical presentation of IV_b.

You can keep this branch for your own records, but the recommended path for new readers and tests is to **start with `bcqm_toy_3`**.

### 4.3 Toy model C – `bcqm_toy_3` (canonical control, \(\beta \approx 0\))

**Status:** **Canonical toy** for BCQM IV_b testing.

- Lives in the `bcqm_toy_3/` folder described above.
- Implements a **direction-free**, **\(W_{\text{coh}}\)-blind** kernel for the acceleration on a thread.
- \(W_{\text{coh}}\) appears **only** as a label in the scan, not in the local dynamics itself.

This is your **control model**: if the pipeline is behaving properly, the fitted \(\beta\) should be statistically consistent with zero.

---

## 5. Canonical toy: `bcqm_toy_3` structure

### 5.1 Configuration (`config.py` + YAML)

- `config.py`:
  - Defines a lightweight dataclass that loads parameters from a YAML file.
  - Top-level sections:

    ```yaml
    simulation:
      dt: 1.0
      n_steps: 16384
      n_ensembles: 64
      seed: 12345
      sign_mode: 1    # +1 or -1

    scan:
      wcoh_values: [5.0, 10.0, 20.0, 40.0, 80.0, 160.0]
      label: wcoh_scan_phase1

    output:
      base_dir: outputs
    ```

- Key points:
  - `dt`, `n_steps` control the time resolution and length of each trajectory.
  - `n_ensembles` is the number of independent trajectories per \(W_{\text{coh}}\).
  - `seed` ensures reproducibility.
  - `sign_mode` can be `+1` or `-1`:
    - it flips the sign of the driving noise,
    - but **does not change** the power spectrum of the process (PSD uses \(|\mathrm{FFT}|^2\)), so it does **not** change \(A\) or \(\omega_c\).

### 5.2 Model (`model.py`)

Implements a simple OU-like update rule for acceleration:

\[
a_{n+1} = (1 - \gamma) a_n + \text{sign\_mode} \times \sigma \xi_n,
\]

with:

- \(\xi_n \sim \mathcal{N}(0, 1)\) i.i.d. Gaussian noise,
- \(\gamma\) and \(\sigma\) fixed (independent of \(W_{\text{coh}}\)).

Crucially:

- The **kernel is stationary and does *not* depend on \(W_{\text{coh}}\)**.
- \(W_{\text{coh}}\) only appears in the **scan loop** (see `simulate.py`) as a label so we can probe for spurious scaling.

### 5.3 Spectral analysis (`spectra.py`)

Provides three main functions:

1. `estimate_Sa(time_series, dt)`  
   - Compute a one-sided PSD for the acceleration series using an FFT-based estimator.

2. `estimate_amplitude_and_omega_c(freqs, Sa)`  
   - Compute:
     - the **total noise amplitude** \(A\) via
       \[
       A = \sqrt{\int S_a(\omega)\, d\omega},
       \]
     - a characteristic frequency \(\omega_c\) as the **spectral centroid** of \(S_a(\omega)\).

3. `fit_beta(Wcoh_values, amplitudes)`  
   - Perform a **log–log linear regression** to fit:
     \[
     A(W_{\text{coh}}) \sim C W_{\text{coh}}^{-\beta},
     \]
   - Return:
     - fitted \(\beta\),
     - uncertainty estimate (from the regression),
     - and the fitted intercept \(\log C\) if needed.

### 5.4 Simulation driver (`simulate.py`)

For each \(W_{\text{coh}}\) in the scan:

1. Run `n_ensembles` trajectories using the stationary kernel.
2. Compute and average the PSDs over the ensemble.
3. Extract `A` and `omega_c` using `spectra.py`.
4. Save per-\(W_{\text{coh}}\) results to `.npz` files:

   ```text
   outputs/wcoh_scan_phase1/Wcoh_5.0.npz
   outputs/wcoh_scan_phase1/Wcoh_10.0.npz
   ...
   ```

5. Collect all amplitudes and frequencies into a CSV file:

   ```text
   outputs/wcoh_scan_phase1/amplitude_scaling.csv
   ```

   with columns:

   ```text
   Wcoh, A, omega_c
   ```

### 5.5 Command-line interface (`cli.py`)

Top-level entry point:

```bash
python3 -m bcqm_toy_3.cli run bcqm_toy_3/configs/wcoh_scan_phase1.yml
```

This will:

- load the YAML config,
- run the scan over the specified `wcoh_values`,
- print out a line for each \(W_{\text{coh}}\):

  ```text
  Running ensemble for W_coh = 5.0 ...
    A = ...
    omega_c = ...
  ```

- fit \(\beta\) at the end and print something like:

  ```text
  Wcoh: [  5.  10.  20.  40.  80. 160.]
  A: [2.29 ...]
  Fitted beta: 3.8e-04
  Estimated error: 7.2e-04
  ```

---

## 6. Canonical test run and expected result

Using the example config:

```yaml
simulation:
  dt: 1.0
  n_steps: 16384
  n_ensembles: 64
  seed: 12345
  sign_mode: 1

scan:
  wcoh_values: [5.0, 10.0, 20.0, 40.0, 80.0, 160.0]
  label: wcoh_scan_phase1

output:
  base_dir: outputs
```

and running:

```bash
python3 -m bcqm_toy_3.cli run bcqm_toy_3/configs/wcoh_scan_phase1.yml
```

you observed:

- All amplitudes clustered tightly around:

  ```text
  A ≈ 2.29
  ```

- Fitted exponent:

  ```text
  beta ≈ 3.8e-4 ± 7.2e-4
  ```

which is statistically consistent with:

- **\(\beta = 0\)**
- i.e. **no detectable scaling with \(W_{\text{coh}}\)**.

Flipping the sign:

```yaml
sign_mode: -1
```

produced **identical** amplitudes and \(\omega_c\), as expected, because:

- \(a_n \to -a_n\) flips the sign of the acceleration series, but
- the PSD depends on \(|\mathrm{FFT}(a)|^2\) and is therefore unchanged.

This confirms that `bcqm_toy_3` behaves as a **W_coh-blind control model** and that the pipeline is not inventing a spurious \(\beta\).

---

## 7. How to add new kernels

To test new IV_b ideas:

1. **Clone `model.py` logic** into a new function (or a new module) that:
   - still exposes a clean update rule per time step,
   - but now makes the local dynamics depend on \(W_{\text{coh}}\) in a controlled way.
2. Expose a switch in the config (e.g. `kernel_type: "control"` vs `"wcoh_dependent"`).
3. Re-run the same `wcoh_scan_phase1.yml` (or variants) and observe:
   - how \(A(W_{\text{coh}})\) behaves,
   - what \(\beta\) you obtain,
   - whether the spectrum shows new features (e.g. changes in \(\omega_c\)).

The **key point** is that the **pipeline remains unchanged**:

- only the kernel changes,
- making it easy to attribute any change in \(\beta\) or spectral shape to the new physics encoded in the kernel.

---


## 7. Cluster toy: `bcqm_cluster_toy` COM scaling vs \(N\)

The `bcqm_cluster_toy` package implements a simple **independent-probe
cluster toy** whose only role is to test how the centre-of-mass (COM)
inertial-noise amplitude scales with the number of probes \(N\).

- Each probe follows the same \(W_{\text{coh}}\)-blind OU-like
  acceleration kernel as in `bcqm_toy_3`.
- A cluster of size \(N\) is modelled as \(N\) such probes evolved in
  parallel with independent Gaussian kicks.
- The COM acceleration is
  \[
  a_{\mathrm{COM},n} = \frac{1}{N} \sum_{i=1}^{N} a_{i,n},
  \]
  and we analyse \(a_{\mathrm{COM},n}\) with exactly the same
  spectrum–amplitude machinery as in the single-probe case.

### 7.1 Configuration and driver

The YAML configuration (e.g. `bcqm_cluster_toy/configs/cluster_n_scan.yml`)
has the structure

```yaml
simulation:
  dt: 1.0
  n_steps: 16384
  n_ensembles: 64
  seed: 12345

cluster:
  N_values: [2, 4, 8, 16, 32, 64, 128]
  label: cluster_n_scan_v1

output:
  base_dir: outputs_cluster
```

The driver

```bash
python3 -m bcqm_cluster_toy.cli run bcqm_cluster_toy/configs/cluster_n_scan.yml
```

will, for each \(N\) in `N_values`:

1. generate ensembles of COM trajectories,
2. compute and average the COM spectra,
3. extract \(A_{\mathrm{COM}}(N)\) and \(\omega_{c,\mathrm{COM}}\),
4. write `cluster_N{N}.npz` files under
   `outputs_cluster/cluster_n_scan_v1/`,
5. and collect the amplitudes into
   `outputs_cluster/cluster_n_scan_v1/amplitude_scaling_COM.csv`.

### 7.2 Expected result

For the configuration above we observe:

- a log--log fit of \(A_{\mathrm{COM}}(N)\) vs.~\(N\) with slope
  \(\alpha \approx -0.495\ and \(R^2 \\approx 0.99997\),
- the product \(A_{\\mathrm{COM}}(N)\\sqrt{N}\) is nearly constant,
  \(\langle A_{\\mathrm{COM}}(N)\\sqrt{N}\\rangle \\approx 17.6\) with
  a scatter below \(1\,\%\),
- the COM spectral centroid \(\\omega_{c,\\mathrm{COM}}\\) remains
  essentially independent of \(N\) at the few-per-mille level.

This is exactly what we expect for \(N\) independent probes: the COM
acceleration is an average over \(N\) identical and uncorrelated
contributions, so its variance (and therefore its noise amplitude)
should fall as \(1/\\sqrt{N}\). In the BCQM IV\_b context, this branch
therefore plays the same role for **cluster COM noise** that
`bcqm_toy_3` plays for **single-probe noise**: it provides a clean
independent-probe baseline, demonstrating that any departure from
\(N^{-1/2}\) suppression in future correlated or \(W_{\\text{coh}}\)-sensitive
cluster kernels can be attributed to the modified microdynamics rather
than to analysis artefacts.


## 8. Relation to the BCQM IV_b text

In the IV_b manuscript, `bcqm_toy_3` is described as:

- a **minimal, graph-level, W_coh-blind kernel** on a thread,
- a **control experiment** whose only job is to demonstrate that:
  - the numerical machinery is stable and reproducible,
  - a stationary kernel with no \(W_{\text{coh}}\) dependence yields \(\beta \approx 0\).

Later IV_b phases (and subsequent papers) are reserved for:

- genuinely \(W_{\text{coh}}\)-sensitive kernels,
- many-thread / entangled-cluster extensions,
- and eventually the full **BCQM inertial-noise prediction**.

For now, this TESTING.md records the status of the canonical `bcqm_toy_3` control toy, the independent-probe cluster baseline `bcqm_cluster_toy`, and the archived earlier toys that were used to validate the pipeline.