#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from pyotf.vector import optical_sectioning_ratio, richards_wolf_linearly_polarized_psf


def test_richards_wolf_reference_psf_is_normalized():
    x = np.linspace(-200, 200, 17)
    y = np.linspace(-200, 200, 17)
    z = np.linspace(-300, 300, 9)
    psf = richards_wolf_linearly_polarized_psf(
        wl=520,
        na=1.1,
        ni=1.33,
        x=x,
        y=y,
        z=z,
        n_theta=36,
        n_phi=48,
    )
    assert psf.shape == (len(z), len(y), len(x))
    np.testing.assert_allclose(psf.sum(), 1.0, rtol=1e-5, atol=1e-8)


def test_richards_wolf_aberration_modifies_psf():
    x = np.linspace(-200, 200, 15)
    y = np.linspace(-200, 200, 15)
    z = np.linspace(-250, 250, 7)
    psf0 = richards_wolf_linearly_polarized_psf(
        wl=520, na=1.0, ni=1.33, x=x, y=y, z=z, n_theta=30, n_phi=40
    )
    pcoefs = np.zeros(12)
    pcoefs[3] = 0.25
    psf1 = richards_wolf_linearly_polarized_psf(
        wl=520, na=1.0, ni=1.33, x=x, y=y, z=z, pcoefs=pcoefs, n_theta=30, n_phi=40
    )
    assert not np.allclose(psf0, psf1)


def test_optical_sectioning_ratio_positive():
    x = np.linspace(-150, 150, 11)
    y = np.linspace(-150, 150, 11)
    z = np.linspace(-200, 200, 9)
    psf = richards_wolf_linearly_polarized_psf(
        wl=520, na=1.0, ni=1.33, x=x, y=y, z=z, n_theta=24, n_phi=32
    )
    ratio = optical_sectioning_ratio(psf, dz=2)
    assert ratio > 0
