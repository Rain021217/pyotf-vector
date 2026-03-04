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


def test_confocal_supports_named_aberration_api():
    confocal = ConfocalMicroscope(
        wl_exc=488,
        pinhole_size=1.0,
        aberrations_em={"defocus": 0.25},
        aberrations_exc={"vertical astigmatism": 0.12},
        **BASE_KWARGS,
    )
    assert confocal.PSF.shape == confocal.model.PSFi.shape


def test_confocal_detector_and_object_modes_match_normalization():
    confocal_obj = ConfocalMicroscope(wl_exc=488, pinhole_size=1.2, pinhole_mode="object", **BASE_KWARGS)
    confocal_det = ConfocalMicroscope(wl_exc=488, pinhole_size=1.2, pinhole_mode="detector", **BASE_KWARGS)
    np.testing.assert_allclose(confocal_obj.PSF.sum(), 1.0)
    np.testing.assert_allclose(confocal_det.PSF.sum(), 1.0)


def test_fourpi_confocal_produces_normalized_psf():
    model = FourPiConfocalMicroscope(wl_exc=488, pinhole_size=0.8, **BASE_KWARGS)
    assert model.PSF.shape == model.model.PSFi.shape
    np.testing.assert_allclose(model.PSF.sum(), 1.0)


def test_fourpi_phase_scan_returns_stack():
    model = FourPiConfocalMicroscope(wl_exc=488, pinhole_size=0.8, **BASE_KWARGS)
    phases = np.linspace(0, 2 * np.pi, 5)
    stack = model.phase_scan(phases, channel="both")
    assert stack.shape[0] == len(phases)
    assert stack.shape[1:] == model.model.PSFi.shape


def test_fourpi_phase_changes_axial_profile():
    model = FourPiConfocalMicroscope(wl_exc=488, pinhole_size=0.8, **BASE_KWARGS)
    model.phase_exc = 0.0
    psf0 = model.model_psf
    model.phase_exc = np.pi
    psf1 = model.model_psf
    assert not np.allclose(psf0, psf1)


def test_fourpi_arm_specific_aberrations_and_amplitude_ratio_work():
    model = FourPiConfocalMicroscope(
        wl_exc=488,
        pinhole_size=0.8,
        amp_ratio_exc=0.7,
        amp_ratio_det=1.3,
        aberrations_exc_arm2={"vertical coma": 0.18},
        aberrations_det_arm1={"oblique astigmatism": 0.1},
        **BASE_KWARGS,
    )
    np.testing.assert_allclose(model.PSF.sum(), 1.0)
