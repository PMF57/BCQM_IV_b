# BCQM_IV_b — Inertial Noise Toy Models for BCQM IV\_b

This repository contains the numerical toy models and testing pipeline
used in the **Boundary–Condition Quantum Mechanics IV\_b** work on
inertial noise and the coherence horizon \(W_{\mathrm{coh}}\).

The goal here is **not** to provide a final physical kernel, but to:

- validate the **simulation → spectrum → amplitude → scaling** pipeline,
- show that the code correctly detects a non-trivial
  \(W_{\mathrm{coh}}\)-dependence when it is present, and
- provide a clean **\(W_{\mathrm{coh}}\)-blind control model**
  (`bcqm_toy_3`) which yields an amplitude scaling exponent
  \(\beta \approx 0\).

For more detail, see `TESTING.md` and the “Numerical checks – code
validation” appendix of the BCQM IV\_b manuscript.

---

## Repository layout

At the time of writing the repository contains:

```text
BCQM_IV_b/
  bcqm_toy_3/              # canonical, W_coh-blind control toy
  bcqm_cluster_toy/        # independent-probe cluster COM toy
  archive/                 # archived toy models (positive/negative controls)
  outputs/
    wcoh_scan_phase1/      # example run outputs for bcqm_toy_3
  outputs_cluster/
    cluster_n_scan_v1/     # example cluster N-scan outputs
  TESTING.md               # detailed testing notes
  LICENSE                  # license information
  .gitattributes
  .DS_Store                # (macOS metadata; can be ignored)
```

### `bcqm_toy_3/` — canonical control toy

This is the **current** and **canonical** toy model used in BCQM IV\_b.

- Implements a simple, stationary acceleration kernel on a single
  thread:
  - an Ornstein–Uhlenbeck–type update
    \(a_{n+1} = (1-\gamma)a_n + \sigma s \xi_n\),
    with fixed \(\gamma\), \(\sigma\), and sign mode \(s\),
    and Gaussian kicks \(\xi_n\).
- Crucially, the kernel is **independent** of \(W_{\mathrm{coh}}\):
  - \(W_{\mathrm{coh}}\) appears only in the **scan**, as a label
    for different ensembles.
- The code:
  - generates ensembles of acceleration trajectories for each
    \(W_{\mathrm{coh}}\) in a scan list,
  - computes one-sided power spectral densities \(S_a(\omega)\),
  - extracts a total amplitude
    \(A(W_{\mathrm{coh}}) = (\int S_a(\omega)\,d\omega)^{1/2}\)
    and a characteristic frequency \(\omega_c\),
  - saves the results to `outputs/`,
  - and performs a log–log fit
    \(A(W_{\mathrm{coh}}) \propto W_{\mathrm{coh}}^{-\beta}\).

In a representative configuration (as described in `TESTING.md`), the
measured amplitudes are essentially constant and the fit returns
\(\beta\) consistent with zero, as expected for a
\(W_{\mathrm{coh}}\)-blind kernel.


### `bcqm_cluster_toy/` — cluster COM toy (independent baseline)

This sibling package implements a simple **independent-probe cluster toy**
used in BCQM IV\_b Section~6 and Appendix~A to test how the
centre-of-mass (COM) inertial-noise amplitude scales with the number of
probes.

- Each probe thread follows the **same** $W_{\mathrm{coh}}$-blind
  Ornstein–Uhlenbeck–type kernel as in `bcqm_toy_3`.
- A cluster of size $N$ is modelled as $N$ such probes evolved in
  parallel with independent Gaussian kicks.
- The only collective observable is the COM acceleration,
  $$
  a_{\mathrm{COM},n} = \frac{1}{N} \sum_{i=1}^{N} a_{i,n},
  $$
  so this branch provides a clean **independent-probe baseline**:
  any suppression of COM noise arises purely from averaging.

The `cluster_simulate.py` driver:

- takes a list of cluster sizes
  $N \in \{2,4,8,16,32,64,128\}$,
- generates COM trajectories using the same $(\Delta t, N_\text{steps},
  n_\text{ensembles})$ as the single-probe scan,
- computes ensemble-averaged COM spectra and extracts a COM amplitude
  $A_{\mathrm{COM}}(N)$ and centroid $\omega_{c,\mathrm{COM}}$,

and writes:

- `cluster_N{N}.npz` files with $(\omega, S_{a,\mathrm{COM}})$, and
- an `amplitude_scaling_COM.csv` summarising
  $(N, A_{\mathrm{COM}}, \omega_{c,\mathrm{COM}})$.

In a representative configuration (matching the IV\_b figures), a
log--log fit of $A_{\mathrm{COM}}(N)$ vs.~$N$ gives a slope
$\alpha \simeq -0.495$ with $R^2 \simeq 0.99997$, and the product
$A_{\mathrm{COM}}(N) \sqrt{N}$ is constant at the $\sim 1\,\%$
level across $N=2$--$128$, confirming the expected $N^{-1/2}$ COM
suppression for independent probes.

### `archive/` — historical toy models (optional controls)

The `archive/` directory is reserved for **historical toy models** that
were used to validate the pipeline, but are **not** part of the
canonical IV\_b narrative. Typical contents include:

- a *toy dynamics* branch in which \(W_{\mathrm{coh}}\) is built
  directly into the step statistics, producing a non-trivial power-law
  \(A(W_{\mathrm{coh}})\) scaling (positive control);
- a *rudder* toy that implements a minimal BCQM-style memory rule with
  interruptions, showing that not every appearance of \(W_{\mathrm{coh}}\)
  produces a clean amplitude power law (negative control).

These models are documented in `TESTING.md` and the IV\_b appendix, but
new users only need `bcqm_toy_3` to reproduce the main control result.

### `outputs/wcoh_scan_phase1/`

This directory contains example outputs from a standard scan run of the
canonical toy:

- individual `.npz` files for each \(W_{\mathrm{coh}}\) value
  (averaged spectra and summary statistics),
- an `amplitude_scaling.csv` file with columns

  ```text
  Wcoh, A, omega_c
  ```

You can delete or regenerate these files as needed.

---

## Installation

You need Python 3 with a minimal set of scientific packages. On macOS or
Linux, a typical setup might be:

```bash
python3 -m venv venv
source venv/bin/activate
pip install numpy pyyaml matplotlib
```

(Additional packages such as `scipy` can be installed if you wish to
extend the analysis, but they are not strictly required for the basic
pipeline.)

---

## Running the canonical control scan

From the repository root (`BCQM_IV_B/`), a typical run looks like:

```bash
python3 -m bcqm_toy_3.cli run bcqm_toy_3/configs/wcoh_scan_phase1.yml
```

This will:

1. Load the specified YAML configuration.
2. For each \(W_{\mathrm{coh}}\) in the `wcoh_values` list:
   - generate an ensemble of trajectories,
   - compute and average the acceleration power spectra,
   - extract \(A\) and \(\omega_c\),
   - save the results under `outputs/wcoh_scan_phase1/`.
3. Collect the amplitudes into `outputs/wcoh_scan_phase1/amplitude_scaling.csv`.
4. Perform a log–log fit
   \(A(W_{\mathrm{coh}}) \propto W_{\mathrm{coh}}^{-\beta}\)
   and print the fitted \(\beta\) and an uncertainty estimate.

Example output (numbers indicative only):

```text
Running ensemble for W_coh = 5.0 ...
  A = 2.30, omega_c = 0.26
...
Fitted beta: 3.8e-04
Estimated error: 7.2e-04
```

## Running the cluster N-scan

To reproduce the simple cluster results used in BCQM IV\_b Section~6,
run the cluster driver from the repository root:

```bash
python3 -m bcqm_cluster_toy.cli run bcqm_cluster_toy/configs/cluster_n_scan.yml
```

This will:

1. Load the specified YAML configuration for the cluster scan.
2. For each cluster size $N$ in the `N_values` list:
   - generate ensembles of COM acceleration trajectories,
   - compute and average the COM power spectra,
   - extract $A_{\mathrm{COM}}(N)$ and $\omega_{c,\mathrm{COM}}$,
   - save the spectra to `outputs_cluster/cluster_n_scan_v1/cluster_N{N}.npz`.
3. Collect the COM amplitudes into
   `outputs_cluster/cluster_n_scan_v1/amplitude_scaling_COM.csv`.

In the configuration used for the IV\_b figures, a log--log fit of
$A_{\mathrm{COM}}(N)$ against $N$ yields a slope
$\alpha \approx -0.495$ with $R^2 \approx 0.99997$, and the product
$A_{\mathrm{COM}}(N)\sqrt{N}$ is constant at the $\sim 1\,\%$
level across $N=2$--$128$, as expected for an independent-probe COM
baseline.


---

## Relation to BCQM IV\_b

This repository is intended as a **companion code base** for the
BCQM IV\_b manuscript on inertial noise and the coherence horizon:

- `bcqm_toy_3` corresponds to the **canonical control toy** described
  in the main text and in the appendix “Numerical checks – code
  validation”.
- The archived toys in `archive/` correspond to the **positive and
  negative controls** also described in that appendix.
- `TESTING.md` mirrors and expands on the appendix, providing a
  prose description of:
  - the purpose of each toy,
  - the pipeline structure,
  - and the expected qualitative behaviour of \(A(W_{\mathrm{coh}})\).

Once BCQM IV\_b is publicly available, you may wish to add its DOI and
citation here.

---

## License

See `LICENSE` for full details.

---

## Contact

For questions, comments, or bug reports related to this code or the
BCQM IV\_b work, please contact the author of the BCQM series.
