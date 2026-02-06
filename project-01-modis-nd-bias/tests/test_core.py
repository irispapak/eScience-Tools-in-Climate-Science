import numpy as np
import xarray as xr

from modis_nd_bias.core import (
    calculate_bias_nd_and_homogeneity,
    calculate_droplet_concentration,
)


def test_calculate_bias_nd_and_homogeneity():
    nd = np.array([[10.0, 12.0], [20.0, 18.0]])
    lwp_mean = np.array([[100.0, 100.0], [200.0, 200.0]])
    lwp_std = np.array([[10.0, 10.0], [40.0, 40.0]])

    bias, homog = calculate_bias_nd_and_homogeneity(
        nresol=2,
        nreg=2,
        nd_regmean=nd,
        lwp_regmean=lwp_mean,
        lwp_regstd=lwp_std,
    )

    assert bias.shape == (2, 2)
    assert homog.shape == (2,)
    # baseline (resolution index 0) bias should be 0
    assert np.allclose(bias[:, 0], 0.0)
    # homogeneity = std/mean at 1km
    assert np.allclose(homog, np.array([0.1, 0.2]))


def test_calculate_droplet_concentration_shapes():
    re = xr.DataArray(np.full((2, 2), 10.0))
    tau = xr.DataArray(np.full((2, 2), 8.0))

    nd = calculate_droplet_concentration(re, tau)
    assert nd.shape == re.shape
    assert np.all(np.isfinite(nd.values))
    assert np.all(nd.values > 0)
