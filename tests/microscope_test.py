#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np

from pyotf.microscope import ConfocalMicroscope, FourPiConfocalMicroscope


BASE_KWARGS = dict(
    model="hanser",
    na=1.1,
    ni=1.33,
    wl=520,
    size=32,
    pixel_size=80,
    oversample_factor=1,
)


def test_confocal_small_pinhole_reduces_out_of_focus_signal():
    confocal_open = ConfocalMicroscope(wl_exc=488, pinhole_size=2.0, **BASE_KWARGS)
    confocal_small = ConfocalMicroscope(wl_exc=488, pinhole_size=0.4, **BASE_KWARGS)

    psf_open = confocal_open.PSF
    psf_small = confocal_small.PSF

    z_center = psf_open.shape[0] // 2
    center_open = psf_open[z_center].sum()
    center_small = psf_small[z_center].sum()
    off_axis_open = psf_open[z_center + 2].sum()
    off_axis_small = psf_small[z_center + 2].sum()

    # Smaller pinhole should improve optical sectioning (higher on/off contrast).
    assert off_axis_small / center_small < off_axis_open / center_open


def test_confocal_excitation_psf_matches_grid_shape():
    confocal = ConfocalMicroscope(wl_exc=488, pinhole_size=1.0, **BASE_KWARGS)
    assert confocal.excitation_psf.shape == confocal.model.PSFi.shape


def test_fourpi_confocal_produces_normalized_psf():
    model = FourPiConfocalMicroscope(wl_exc=488, pinhole_size=0.8, **BASE_KWARGS)
    assert model.PSF.shape == model.model.PSFi.shape
    np.testing.assert_allclose(model.PSF.sum(), 1.0)
