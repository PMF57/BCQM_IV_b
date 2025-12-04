# `bcqm_cluster_toy` — independent-probe cluster COM toy

This package implements a simple **independent-probe cluster toy** used in
BCQM IV_b (Section 6 and Appendix A) to test how the **centre-of-mass (COM)
inertial-noise amplitude** scales with the number of probes \(N\).

The key points:

- Each probe thread follows the **same** \(W_{\mathrm{coh}}\)-blind
  Ornstein–Uhlenbeck–type acceleration kernel as in the canonical single-probe
  toy `bcqm_toy_3`.
- A cluster of size \(N\) is modelled as \(N\) such probes evolved in
  parallel with independent Gaussian kicks.
- The only collective observable is the COM acceleration
  \[
  a_{\mathrm{COM},n} = \frac{1}{N} \sum_{i=1}^{N} a_{i,n},
  \]
  so any suppression of COM noise arises purely from averaging over independent
  contributions. This provides a **clean baseline**: if future, more physical
  cluster kernels show deviations from \(N^{-1/2}\) behaviour, the cause lies
  in the modified microdynamics rather than in the analysis pipeline.

---

## Layout

Inside the repository this package typically sits as:

```text
BCQM_IV_B/
  bcqm_toy_3/           # canonical W_coh-blind single-probe control toy
  bcqm_cluster_toy/     # this independent-probe cluster COM toy
  archive/              # earlier experimental toys (positive/negative controls)
  ...
```

The main files here are:

- `config.py` – YAML → `Config` loader (simulation, cluster scan, and output).
- `cluster_model.py` – N-body OU-like kernel for N probes.
- `cluster_spectra.py` – one-sided PSD and amplitude extraction for COM.
- `cluster_simulate.py` – top-level N-scan routine.
- `cli.py` – command-line entrypoint.
- `configs/cluster_n_scan.yml` – example configuration.

---

## Requirements

This package uses only standard scientific Python libraries:

- Python 3.9+
- `numpy`
- `pyyaml`

If you can run `bcqm_toy_3`, you should already have everything needed. From a
fresh environment:

```bash
pip install numpy pyyaml
```

---

## Running the cluster N-scan

From the **repository root** (the directory that contains `bcqm_cluster_toy/`),
run:

```bash
python3 -m bcqm_cluster_toy.cli run bcqm_cluster_toy/configs/cluster_n_scan.yml
```

The example configuration `configs/cluster_n_scan.yml` has the form:

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

This tells the driver to:

1. Use `dt = 1.0`, `n_steps = 16384`, `n_ensembles = 64`, `seed = 12345`.
2. Loop over cluster sizes \(N \in \{2,4,8,16,32,64,128\}\).
3. For each \(N\):
   - simulate ensembles of COM acceleration trajectories,
   - compute ensemble-averaged COM spectra,
   - extract the COM amplitude \(A_{\mathrm{COM}}(N)\) and COM centroid
     \(\omega_{c,\mathrm{COM}}\),
   - save the spectra as `cluster_N{N}.npz` under
     `outputs_cluster/cluster_n_scan_v1/`.

It also writes a summary CSV:

- `outputs_cluster/cluster_n_scan_v1/amplitude_scaling_COM.csv`

with columns:

- `N` – cluster size,
- `A_com` – COM amplitude \(A_{\mathrm{COM}}(N)\),
- `omega_c_com` – COM spectral centroid \(\omega_{c,\mathrm{COM}}\).

---

## Expected behaviour

For the configuration above (the one used to generate Fig. 4 in BCQM IV_b), a
log–log fit of \(A_{\mathrm{COM}}(N)\) vs. \(N\) yields:

- slope \(\alpha \approx -0.495\),
- coefficient of determination \(R^2 \approx 0.99997\),

and the product \(A_{\mathrm{COM}}(N)\sqrt{N}\) is nearly constant:

- \(\langle A_{\mathrm{COM}}(N)\sqrt{N} \rangle \approx 17.6\)
- with a scatter below \(1\,\%\) across \(N = 2\)–128.

This is exactly the expected **\(N^{-1/2}\) suppression** of COM noise for
independent probes. In the BCQM IV_b programme, `bcqm_cluster_toy` therefore
plays the same role for **cluster COM noise** that `bcqm_toy_3` plays for
**single-probe noise**:

- it is a **null / baseline model** that demonstrates the analysis pipeline
  recovers the correct scaling when it is present.

---

## Relation to the BCQM IV_b paper

The results from this package are used in:

- **Section 6** (“Results: simple cluster models and centre-of-mass noise”)  
  — Fig. 4 (COM amplitude vs. N on log–log axes),  
- **Appendix A** (“Numerical checks – code validation”)  
  — subsection describing the `bcqm_cluster_toy` branch and its outputs.

If you reproduce the N-scan with the example configuration, you should obtain
numbers compatible with those quoted there (within expected stochastic
variation).

---

## Licence

This package inherits the licence of the parent repository (see `LICENSE` at
the repository root).
