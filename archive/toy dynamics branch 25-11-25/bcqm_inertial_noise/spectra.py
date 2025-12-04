"""
Spectral estimation and scaling tools for bcqm_inertial_noise.

The goal is to estimate one-sided acceleration-noise spectra
S_a(omega; W_coh) from discrete acceleration time series, and to provide
helpers for amplitude scaling and data collapse.

This module intentionally uses only NumPy and basic windowing to keep
dependencies light. The numerical normalisation can be refined later;
for IV_b the emphasis is on consistent shapes and relative scaling
rather than absolute calibration.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple, Sequence

import numpy as np


@dataclass
class SpectrumResult:
    omega: np.ndarray  # angular frequencies
    Sa: np.ndarray  # one-sided acceleration noise spectrum
    window: str
    segment_length: int
    overlap: float


def _make_window(kind: str, n: int) -> np.ndarray:
    kind = kind.lower()
    if kind == "hann":
        return np.hanning(n)
    elif kind in {"rect", "rectangular", "boxcar"}:
        return np.ones(n)
    else:
        raise ValueError(f"Unknown window type: {kind!r}")


def estimate_Sa(
    acc: np.ndarray,
    dt: float,
    segment_length: int,
    overlap: float = 0.5,
    window: str = "hann",
) -> SpectrumResult:
    """
    Estimate a one-sided acceleration-noise spectrum from time series
    data using a simple Welch-style averaging procedure.

    Parameters
    ----------
    acc:
        Acceleration time series. Can be 1D (single realisation) or 2D
        with shape (n_real, n_samples).
    dt:
        Time step between samples.
    segment_length:
        Length of segments used in the Welch estimate.
    overlap:
        Fraction of overlap between successive segments (0 < overlap < 1).
    window:
        Name of the window function to use (currently 'hann' or
        'rect').

    Returns
    -------
    SpectrumResult
        Frequencies (omega) and estimated one-sided spectrum Sa(omega).
    """
    acc = np.asarray(acc, dtype=float)
    if acc.ndim == 1:
        acc = acc[None, :]  # add realisation axis

    n_real, n_samples = acc.shape
    if segment_length > n_samples:
        segment_length = n_samples

    step = int(segment_length * (1.0 - overlap))
    if step <= 0:
        raise ValueError("segment_length * (1-overlap) must be >= 1")

    win = _make_window(window, segment_length)
    win_power = np.sum(win**2)

    # Accumulate power spectra
    psd_accum = None
    n_segments_total = 0

    for r in range(n_real):
        x = acc[r]
        start = 0
        while start + segment_length <= n_samples:
            seg = x[start:start+segment_length]
            seg = seg - np.mean(seg)
            seg_win = seg * win

            # FFT and power
            fft_vals = np.fft.rfft(seg_win)
            psd = (np.abs(fft_vals)**2) / (win_power)
            if psd_accum is None:
                psd_accum = psd
            else:
                psd_accum += psd
            n_segments_total += 1
            start += step

    if psd_accum is None or n_segments_total == 0:
        raise ValueError("Not enough data to form at least one segment")

    psd_mean = psd_accum / n_segments_total

    # Frequency axis: f_k = k / (dt * segment_length), omega = 2*pi*f
    freqs = np.fft.rfftfreq(segment_length, d=dt)
    omega = 2.0 * np.pi * freqs

    # Convert to one-sided PSD for real signals: factor 2 on non-DC, non-Nyquist
    Sa = psd_mean.copy()
    if Sa.size > 2:
        Sa[1:-1] *= 2.0

    Sa = Sa * dt  # approximate scaling to S(omega) with units consistent with dt

    return SpectrumResult(omega=omega, Sa=Sa, window=window, segment_length=segment_length, overlap=overlap)


def estimate_amplitude(
    spectrum: SpectrumResult,
    ref_index: int,
) -> float:
    """
    Extract a simple amplitude A_a from a spectrum by sampling Sa at a
    reference index.

    This is a placeholder for more sophisticated definitions (e.g.
    band-integrated power), but is sufficient to define scaling tests
    of A_a(W_coh) as a function of the coherence horizon.
    """
    if not (0 <= ref_index < spectrum.Sa.size):
        raise IndexError(f"ref_index {ref_index} out of range for spectrum of length {spectrum.Sa.size}")
    return float(spectrum.Sa[ref_index])


def fit_amplitude_scaling(
    Wcoh_list: Sequence[float],
    A_list: Sequence[float],
) -> Tuple[float, float]:
    """
    Fit A_a(W_coh) ~ W_coh^{-beta} in log-log space.

    Parameters
    ----------
    Wcoh_list:
        Sequence of coherence horizons.
    A_list:
        Corresponding amplitudes extracted from spectra.

    Returns
    -------
    beta, beta_err
        Fitted exponent and a rough uncertainty estimate, obtained from
        the covariance matrix of a linear least-squares fit in log space.
    """
    W = np.asarray(Wcoh_list, dtype=float)
    A = np.asarray(A_list, dtype=float)
    if W.size != A.size:
        raise ValueError("Wcoh_list and A_list must have the same length")

    logW = np.log(W)
    logA = np.log(A)

    # Linear fit: logA = c - beta * logW
    A_mat = np.vstack([np.ones_like(logW), -logW]).T
    coeffs, residuals, rank, s = np.linalg.lstsq(A_mat, logA, rcond=None)
    c, beta = coeffs

    # Rough uncertainty: from covariance matrix
    if W.size > 2 and residuals.size == 1:
        sigma2 = residuals[0] / (W.size - 2)
        cov = sigma2 * np.linalg.inv(A_mat.T @ A_mat)
        beta_err = float(np.sqrt(cov[1, 1]))
    else:
        beta_err = np.nan

    return float(beta), beta_err
