# spectra.py
from __future__ import annotations

import math
from typing import Tuple

import numpy as np


def estimate_Sa(acceleration: np.ndarray, dt: float) -> Tuple[np.ndarray, np.ndarray]:
    """
    Estimate a one-sided power spectral density for the acceleration series.

    Returns (frequencies, Sa).
    """
    n = len(acceleration)
    # Remove DC component
    a = acceleration - np.mean(acceleration)

    freqs = np.fft.rfftfreq(n, d=dt)
    fft_vals = np.fft.rfft(a)

    # Simple PSD estimate; normalisation chosen to be consistent across runs
    Sa = (2.0 * dt**2 / (n * dt)) * (np.abs(fft_vals) ** 2)

    return freqs, Sa


def estimate_amplitude_and_omega_c(
    freqs: np.ndarray,
    Sa: np.ndarray,
) -> Tuple[float, float]:
    """
    Estimate a single amplitude parameter A and a characteristic frequency omega_c
    from the PSD.

    Here we use:
      - A = sqrt( integral Sa df )  (overall scale)
      - omega_c = spectral centroid of Sa, converted to angular frequency
    """
    if np.all(Sa <= 0):
        return 0.0, 0.0

    df = freqs[1] - freqs[0] if len(freqs) > 1 else 1.0

    # Overall scale
    integral = float(np.sum(Sa) * df)
    A = math.sqrt(max(integral, 0.0))

    # Spectral centroid (avoid zero frequency in the weight)
    omega = 2.0 * math.pi * freqs
    weight = Sa.copy()
    weight[0] = 0.0
    denom = float(np.sum(weight) * df)
    if denom <= 0.0:
        omega_c = 0.0
    else:
        omega_c = float(np.sum(weight * omega) * df / denom)

    return A, omega_c


def fit_beta(Wcoh: np.ndarray, A: np.ndarray) -> Tuple[float, float]:
    """
    Fit A(W_coh) ~ C * W_coh^{-beta} on a log-log scale.

    Returns (beta, beta_error).
    """
    Wcoh = np.asarray(Wcoh, dtype=float)
    A = np.asarray(A, dtype=float)

    mask = (Wcoh > 0.0) & (A > 0.0)
    if not np.any(mask):
        return float("nan"), float("nan")

    W = Wcoh[mask]
    Am = A[mask]

    logW = np.log(W)
    logA = np.log(Am)

    # Linear fit: logA = intercept + slope * logW
    # => A ~ exp(intercept) * W^{slope}, so beta = -slope
    p, cov = np.polyfit(logW, logA, 1, cov=True)
    slope = p[0]
    beta = -float(slope)

    # Uncertainty on slope -> uncertainty on beta
    if cov is None:
        beta_err = float("nan")
    else:
        beta_err = float(math.sqrt(cov[0, 0]))

    return beta, beta_err
