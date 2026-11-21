import json

import numpy as np
import pytest
import yaml
from scipy import ndimage as ndi

from morphoct import foram
from morphoct.cli import main
from morphoct.io import Volume, save_volume


@pytest.mark.slow
def test_cli_end_to_end_grey(tmp_path, blobs):
    rng = np.random.default_rng(0)
    gray = (np.where(blobs == 1, 60, 180) + rng.normal(0, 20, blobs.shape)).clip(0, 255).astype(np.uint8)
    save_volume(tmp_path / "g.tif", Volume(gray, 5.0))
    cfg = {"name": "t", "output": str(tmp_path / "out"),
           "input": {"path": str(tmp_path / "g.tif"), "voxel_size": 5.0},
           "preprocess": {"denoise": "median", "size": 3}, "segment": {"method": "otsu"},
           "rve": {"windows": [4, 6, 8, 12, 16, 24]},
           "properties": {"subvolume": 24, "axes": [0], "max_steps": 3000, "tol": 1e-3},
           "network": {"axes": [0]}}
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(cfg))
    assert main(["run", str(tmp_path / "c.yaml")]) == 0
    res = json.loads((tmp_path / "out" / "results.json").read_text())
    assert abs(res["morphometrics"]["porosity"] - blobs.mean()) < 0.03
    for f in ("report.html", "report.pdf", "qc_segmentation.png", "rev_porosity.png", "binary.npz"):
        assert (tmp_path / "out" / f).exists()


def _hollow_cube(n=20, outer=8, inner=4, origin=4, voxel=2.0):
    """Closed cubic shell of edge ``outer`` with a cubic cavity of edge ``inner``."""
    m = np.zeros((n, n, n), dtype=np.uint8)
    o0, o1 = origin, origin + outer
    i0, i1 = origin + (outer - inner) // 2, origin + (outer - inner) // 2 + inner
    m[o0:o1, o0:o1, o0:o1] = 1
    m[i0:i1, i0:i1, i0:i1] = 0
    n_shell = outer**3 - inner**3
    return m, n_shell, voxel


def test_foram_helpers_and_stubs():
    cal = foram.linear_density_calibration([10, 200], [0.0, 2.71])
    assert abs(cal(200) - 2.71) < 1e-9
    m = np.zeros((10, 10, 10), bool)
    m[2:8, 2:8, 2:8] = True
    assert foram.shell_volume(m, 2.0) == 216 * 8
    g = np.where(m, 150.0, 0.0)
    assert foram.density_index(g, m) == 150.0
    with pytest.raises(NotImplementedError):
        foram.segment_chambers(m)
    with pytest.raises(NotImplementedError):
        foram.separate_interior(m)
    with pytest.raises(NotImplementedError):
        foram.wall_thickness(m, 1.0)


def test_hollow_cube_volumes_and_surface():
    m, n_shell, vs = _hollow_cube()
    assert foram.shell_volume(m, vs) == n_shell * vs**3
    filled = n_shell + 4**3
    assert foram.filled_test_volume(m, vs, close_radius=0) == filled * vs**3
    assert foram.convex_hull_volume(m, vs) == 8**3 * vs**3
    # exposed faces: outer 6*8^2 plus inner 6*4^2
    assert foram.shell_surface_area(m, vs, method="voxel_faces") == (6 * 8**2 + 6 * 4**2) * vs**2
    sa = foram.shell_surface_area(m, vs, method="marching_cubes")
    # marching cubes is smoother than the voxel-face count but of the same order
    vf = (6 * 8**2 + 6 * 4**2) * vs**2
    assert 0.5 * vf < sa < vf


def test_segment_test_and_analyse_on_synthetic_shell():
    m, n_shell, vs = _hollow_cube()
    rng = np.random.default_rng(0)
    grey = np.where(m, 180.0, 20.0) + rng.normal(0, 5, m.shape)
    shell, thr, info = foram.segment_test(grey, min_size=8, largest_only=True)
    assert 20 < thr < 180
    # allow a few partial-volume / noise voxels
    assert abs(int(shell.sum()) - n_shell) <= 10
    res = foram.analyse(grey, vs, min_size=8, close_radius=0)
    assert res.threshold == thr
    assert res.filled_volume_um3 > res.shell_volume_um3 > 0
    assert 0.05 < res.test_porosity < 0.25
    assert 0.75 < res.calcification_fraction < 0.95
    assert res.surface_area_um2 > 0
    assert res.density_index is not None
    d = res.as_dict()
    assert "shell_volume_um3" in d
    # fill-holes on the clean mask recovers the outer cube
    assert ndi.binary_fill_holes(m.astype(bool)).sum() == 8**3


def test_aperture_needs_closing_and_hull():
    """A one-face opening stops hole filling; closing seals it. Convex hull is the outer cube."""
    m, n_shell, vs = _hollow_cube()
    # open a 2x2 window on the +x face of the outer cube (origin 4, edge 8)
    m[6:8, 6:8, 10:12] = 0  # punch through the two-voxel +x wall
    assert foram.filled_test_volume(m, vs, close_radius=0) == foram.shell_volume(m, vs)
    sealed = foram.filled_test_volume(m, vs, close_radius=2)
    assert sealed >= 8**3 * vs**3 * 0.95
    assert sealed > foram.shell_volume(m, vs)
    assert foram.convex_hull_volume(m, vs) == 8**3 * vs**3
    # grey analyse with default closing should report a non-trivial interior
    grey = np.where(m, 180.0, 20.0)
    res = foram.analyse(grey, vs, min_size=8, close_radius=2)
    assert res.test_porosity > 0.05
    assert res.convex_hull_volume_um3 == 8**3 * vs**3
    assert 0.7 < res.shell_to_hull_fraction < 0.95
