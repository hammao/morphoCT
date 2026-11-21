import numpy as np

from morphoct import morphometrics, rve


def test_morphometrics_blobs(blobs):
    m = morphometrics.compute(blobs, 2.0)
    assert abs(m["porosity"] - blobs.mean()) < 1e-12
    assert m["minkowski"]["S_V_per_um"] > 0
    assert m["connectivity"]["connected_porosity"] <= m["porosity"] + 1e-12
    assert m["two_point_correlation"]["S2"][0] == m["porosity"] or abs(m["two_point_correlation"]["S2"][0] - m["porosity"]) < 1e-9
    # stereology: mean chord 4 phi / S_V within 25 % of the measured mean chord
    assert abs(m["chord_stereology_um"] / m["mean_chord_length_pooled_um"] - 1) < 0.25


def test_specific_surface_of_slab():
    b = np.zeros((40, 40, 40), np.uint8)
    b[:20] = 1  # one flat interface of area 40x40 in a 40^3 box
    m = morphometrics.compute(b, 1.0, pore_sizes=False, chords=False)
    assert abs(m["minkowski"]["S_V_per_um"] - 1 / 40) / (1 / 40) < 0.1


def test_rev_iid_noise_gamma_one():
    rng = np.random.default_rng(0)
    b = (rng.random((64, 64, 64)) < 0.2).astype(np.uint8)
    r = rve.rev_porosity(b, 1.0, windows=[2, 3, 4, 6, 8, 12, 16], fit_from_edge_voxels=2)
    assert abs(r["fit"]["gamma"] - 1) < 0.1
    # integral range of iid voxels is one voxel: K = phi (1 - phi)
    assert abs(r["fit"]["K"] / (0.2 * 0.8) - 1) < 0.2
    assert r["sample_is_rev"]
