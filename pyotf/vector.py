#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Vectorial Richards-Wolf / Debye reference solvers and helpers."""

import numpy as np

from .zernike import noll2degrees, zernike


def zernike_phase_map(rho, phi, pcoefs=None):
    """Return phase aberration map from Noll-ordered Zernike coefficients."""
    if pcoefs is None:
        return np.zeros_like(rho)
    if np.isscalar(pcoefs):
        raise TypeError("`pcoefs` must be a 1D sequence, not a scalar")
    pcoefs = np.atleast_1d(np.asarray(pcoefs))
    if pcoefs.ndim != 1:
        raise TypeError("`pcoefs` must be a 1D sequence")
    if pcoefs.size == 0:
        return np.zeros_like(rho)
    zerns = zernike(rho, phi, *noll2degrees(np.arange(pcoefs.size) + 1))
    return (zerns * pcoefs[:, None, None]).sum(0)


def richards_wolf_linearly_polarized_psf(
    *,
    wl,
    na,
    ni,
    x,
    y,
    z,
    pcoefs=None,
    n_theta=80,
    n_phi=120,
):
    """Compute a vectorial PSF by direct Debye integration (reference implementation).

    Notes
    -----
    This implementation is intentionally explicit and slower than FFT-based pupil methods,
    and is intended for validation / benchmark use.
    """
    if wl <= 0:
        raise ValueError("`wl` must be positive")
    if ni <= 0:
        raise ValueError("`ni` must be positive")
    if na <= 0 or na > ni:
        raise ValueError("`na` must satisfy 0 < na <= ni")
    if n_theta < 2:
        raise ValueError("`n_theta` must be >= 2")
    if n_phi < 1:
        raise ValueError("`n_phi` must be >= 1")

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    xx, yy = np.meshgrid(x, y, indexing="xy")

    k = 2 * np.pi * ni / wl
    alpha = np.arcsin(na / ni)

    theta = np.linspace(0, alpha, n_theta)
    phi = np.linspace(0, 2 * np.pi, n_phi, endpoint=False)
    dtheta = theta[1] - theta[0]
    dphi = phi[1] - phi[0] if n_phi > 1 else 2 * np.pi

    th, ph = np.meshgrid(theta, phi, indexing="ij")
    st, ct = np.sin(th), np.cos(th)
    cp, sp = np.cos(ph), np.sin(ph)

    # x-linearly polarized input transformed by Richards-Wolf matrix (aplanatic objective)
    ax = np.sqrt(ct) * (ct * cp**2 + sp**2)
    ay = np.sqrt(ct) * (ct - 1) * cp * sp
    az = -np.sqrt(ct) * st * cp

    # normalized pupil radius and aberration phase
    rho = st / np.sin(alpha)
    phase_ab = zernike_phase_map(rho, ph, pcoefs=pcoefs)

    jac = st * dtheta * dphi
    common = jac * np.exp(1j * phase_ab)

    out = np.zeros((len(z), len(y), len(x)), dtype=float)
    for iz, zv in enumerate(z):
        kz = np.exp(1j * k * zv * ct)
        phase_xy = np.exp(
            1j
            * k
            * st[..., None, None]
            * (
                xx[None, None] * cp[..., None, None]
                + yy[None, None] * sp[..., None, None]
            )
        )
        wx = (common * ax * kz)[..., None, None] * phase_xy
        wy = (common * ay * kz)[..., None, None] * phase_xy
        wz = (common * az * kz)[..., None, None] * phase_xy
        ex = wx.sum(axis=(0, 1))
        ey = wy.sum(axis=(0, 1))
        ez = wz.sum(axis=(0, 1))
        out[iz] = (np.abs(ex) ** 2 + np.abs(ey) ** 2 + np.abs(ez) ** 2).real
    total = out.sum()
    if total <= 0 or not np.isfinite(total):
        raise ValueError(
            "Computed PSF normalization is invalid; check input parameters"
        )
    return out / total


def axial_profile(psf):
    """Return axial intensity profile by integrating each z plane."""
    return np.asarray(psf).sum(axis=(1, 2))


def optical_sectioning_ratio(psf, dz=2):
    """Compute off-focus to in-focus ratio using integrated z planes."""
    profile = axial_profile(psf)
    z0 = len(profile) // 2
    zidx = z0 + dz
    if not (0 <= zidx < len(profile)):
        lo = -z0
        hi = len(profile) - 1 - z0
        raise ValueError(f"`dz` is out of range. valid range is [{lo}, {hi}]")
    return profile[zidx] / profile[z0]
