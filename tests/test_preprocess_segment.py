import numpy as np

from morphoct import preprocess, segment
from morphoct.io import Volume


def test_mirror_and_crop():
    a = np.arange(8).reshape(2, 2, 2)
    m = preprocess.mirror(a, (0,))
    assert m.shape == (4, 2, 2) and np.array_equal(m[3], a[0])
    v = preprocess.crop(Volume(np.zeros((10, 10, 10))), size=(4, 4, 4), center=True)
    assert v.shape == (4, 4, 4) and v.meta["crop_start"] == [3, 3, 3]


def test_otsu_recovers_two_phase(blobs):
    rng = np.random.default_rng(0)
    gray = np.where(blobs == 1, 60.0, 180.0) + rng.normal(0, 25, blobs.shape)
    vol = preprocess.denoise(Volume(gray.astype(np.float32), 1.0), "median", size=3)
    mask, t = segment.otsu(vol)
    assert 90 < t < 150
    assert segment.dice(mask, blobs) > 0.95
    sens = segment.threshold_sensitivity(vol, t)
    assert sens[0]["porosity"] < sens[-1]["porosity"]


def test_nlm_runs():
    rng = np.random.default_rng(0)
    g = rng.normal(100, 10, (12, 16, 16)).astype(np.float32)
    out = preprocess.denoise(Volume(g), "nlm")
    assert out.data.std() < g.std()


def test_multi_otsu_and_watershed(blobs):
    gray = np.where(blobs == 1, 20.0, 200.0)
    gray[:, :, :8] = 110.0
    labels, ts = segment.multi_otsu(Volume(gray), 3)
    assert len(ts) == 2 and set(np.unique(labels)) <= {0, 1, 2}
    grains = segment.watershed_grains(blobs == 0)
    assert grains.max() > 1


def test_qc_figure(tmp_path, blobs):
    p = segment.qc_figure(blobs.astype(float), blobs, tmp_path / "qc.png", thresholds=[0.5])
    assert p.exists()
