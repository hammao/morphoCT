import numpy as np
import pytest

from morphoct import network, properties


def _slit(nz=8, H=12, h=10, nx=4):
    b = np.zeros((nz, H, nx), np.uint8)
    b[:, 1:1 + h, :] = 1
    return b, H, h


def test_lbm_slit_poiseuille():
    b, H, h = _slit()
    r = properties.permeability(b, 1.0, axes=(0,))
    exact = h**2 / 12 * h / H  # um^2 with 1 um voxels
    assert abs(r["axes"]["z"]["k_um2"] / exact - 1) < 0.03


def test_formation_factor_slit():
    b, H, h = _slit()
    r = properties.formation_factor(b, 1.0, axes=(0,))
    assert abs(r["axes"]["z"]["F"] - H / h) < 0.02 * H / h


def test_kozeny_carman_positive():
    assert properties.kozeny_carman(0.2, 0.05) > 0


@pytest.mark.slow
def test_network_blobs(blobs):
    r = network.permeability(blobs, 5.0, axes=(0,))
    assert r["n_pores"] > 10 and r["axes"]["z"]["k_mD"] > 0
