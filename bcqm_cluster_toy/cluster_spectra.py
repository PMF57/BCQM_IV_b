# bcqm_cluster_toy/cluster_spectra.py

from __future__ import annotations

from typing import Tuple
import numpy as np


def estimate_psd_one_sided(
    timeseries: np.ndarray, dt: float
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Estimate one-sided PSD S_a(omega) from a real-valued time series.

    Parameters
    ----------
    timeseries : np.ndarray
        1D array of length n_steps.
    dt : float
        Time step (dimensionless or physical).

    Returns
    -------
    omega : np.ndarray
        Angular frequencies (same units as 1/dt).
    Sa : np.ndarray
        One-sided power spectral density.
    """
    x = np.asarray(timeseries, dtype=float)
    x = x - x.mean()
    n = x.size

    # Real FFT
    fft_vals = np.fft.rfft(x)
    freqs = np.fft.rfftfreq(n, d=dt)  # in cycles per unit time
    omega = 2.0 * np.pi * freqs      # angular frequency

    # Simple one-sided PSD normalisation:
    # S(omega_k) ≈ (2 dt^2 / T) |X_k|^2 with T = n * dt
    T = n * dt
    Sa = (2.0 * dt**2 / T) * np.abs(fft_vals) ** 2

    return omega, Sa


def estimate_amplitude_and_omega_c(
    omega: np.ndarray, Sa: np.ndarray
) -> Tuple[float, float]:
    """
    Estimate overall amplitude A and spectral centroid omega_c
    from a one-sided PSD S_a(omega).
    """
    # Avoid omega=0 in centroid to prevent division issues
    mask = omega > 0
    omega_nz = omega[mask]
    Sa_nz = Sa[mask]

    # Amplitude: sqrt of integrated PSD
    A = float(np.sqrt(np.trapz(Sa_nz, omega_nz)))

    # Spectral centroid
    num = np.trapz(omega_nz * Sa_nz, omega_nz)
    den = np.trapz(Sa_nz, omega_nz)
    omega_c = float(num / den) if den > 0 else 0.0

    return A, omega_c


def ensemble_com_spectrum(
    a_com: np.ndarray, dt: float
) -> Tuple[np.ndarray, np.ndarray, float, float]:
    """
    Compute average one-sided PSD for COM acceleration over ensembles.

    Parameters
    ----------
    a_com : np.ndarray
        Shape (n_ensembles, n_steps).
    dt : float

    Returns
    -------
    omega : np.ndarray
    Sa_avg : np.ndarray
    A_com : float
    omega_c_com : float
    """
    n_ensembles, _ = a_com.shape
    psd_list = []
    omega = None

    for e in range(n_ensembles):
        omega_e, Sa_e = estimate_psd_one_sided(a_com[e, :], dt)
        if omega is None:
            omega = omega_e
        psd_list.append(Sa_e)

    Sa_stack = np.vstack(psd_list)
    Sa_avg = Sa_stack.mean(axis=0)

    A_com, omega_c_com = estimate_amplitude_and_omega_c(omega, Sa_avg)

    return omega, Sa_avg, A_com, omega_c_com