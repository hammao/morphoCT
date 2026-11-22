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
    # ambient occlusion remains Phase 3; closing is implemented
    with pytest.raises(NotImplementedError):
        foram.separate_interior(m, method="ambient_occlusion")
    interior, info = foram.separate_interior(m, method="closing", close_radius=0)
    assert info["method"] == "closing"
    assert interior.dtype == np.uint8


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


def _ball(shape, centre, radius):
    zz, yy, xx = np.ogrid[:shape[0], :shape[1], :shape[2]]
    z, y, x = centre
    return (zz - z) ** 2 + (yy - y) ** 2 + (xx - x) ** 2 <= radius ** 2


def _cylinder_x(shape, y0, z0, x0, x1, radius):
    zz, yy, xx = np.ogrid[:shape[0], :shape[1], :shape[2]]
    along = (xx >= min(x0, x1)) & (xx <= max(x0, x1))
    return along & ((yy - y0) ** 2 + (zz - z0) ** 2 <= radius ** 2)


def multi_chamber_phantom(vs=1.0, wall=4, air=10.0, calcite=200.0):
    """Three spherical lumens (radii 7, 6, 5) joined by thin necks, coated by a
    uniform wall of ``wall`` voxels. Greys are constant in air and in the shell."""
    from skimage.morphology import ball

    shape = (40, 40, 52)
    # centres along x, with a 1-voxel gap plus a neck of radius 2
    c1, r1 = (20, 20, 12), 7
    c2, r2 = (20, 20, 26), 6
    c3, r3 = (20, 20, 38), 5
    lumen = _ball(shape, c1, r1) | _ball(shape, c2, r2) | _ball(shape, c3, r3)
    lumen |= _cylinder_x(shape, 20, 20, 12, 26, 2)
    lumen |= _cylinder_x(shape, 20, 20, 26, 38, 2)
    shell = ndi.binary_dilation(lumen, ball(wall)) & ~lumen
    grey = np.full(shape, air, dtype=np.float32)
    grey[shell] = calcite
    truth = {
        "n_chambers": 3,
        "radii": (r1, r2, r3),
        "wall_diameter_um": wall * vs,
        "air": air,
        "calcite": calcite,
        "calcite_density": 2.71,
        "voxel_size": vs,
        "shell_voxels": int(shell.sum()),
        "shell_volume_um3": int(shell.sum()) * vs**3,
        "centres": (c1, c2, c3),
    }
    return grey, shell.astype(np.uint8), lumen.astype(np.uint8), truth


def test_segment_chambers_on_three_chamber_phantom():
    grey, shell, lumen, truth = multi_chamber_phantom()
    vs = truth["voxel_size"]
    labels, info = foram.segment_chambers(lumen, vs, smooth_sigma=0.0, h_ratio=0.25,
                                          min_chamber_voxels=30, grey=grey)
    assert info["n_chambers"] == 3, info
    # growth-path order starts at the smallest chamber (radius 5)
    vols = info["volumes_um3"]
    assert vols[0] < vols[1] < vols[2] or vols[0] == min(vols)
    assert abs(vols[0] - (4 / 3 * np.pi * 5**3)) / vols[0] < 0.25
    # labels 1..3 cover the lumen
    assert int(labels.max()) == 3
    assert set(np.unique(labels[lumen.astype(bool)])) <= {1, 2, 3}


def test_wall_thickness_known_diameter():
    grey, shell, lumen, truth = multi_chamber_phantom(wall=4, vs=1.0)
    thick, stats = foram.wall_thickness(shell, 1.0)
    # Hildebrand diameter of a 4-voxel coating is 4 um; allow a voxel of discrete error
    assert abs(stats["median"] - 4.0) <= 1.0, stats
    assert abs(stats["mean"] - 4.0) <= 1.2, stats
    assert 3.0 <= stats["max"] <= 5.0, stats
    assert thick[shell.astype(bool)].max() > 0
    # a spherical shell of radial thickness 4 (independent construction)
    z, y, x = np.ogrid[:32, :32, :32]
    r = np.sqrt((z - 16) ** 2 + (y - 16) ** 2 + (x - 16) ** 2)
    sph = (r <= 12) & (r >= 8)
    t2, s2 = foram.wall_thickness(sph, 1.0)
    assert abs(s2["median"] - 4.0) <= 1.0, s2
    assert abs(s2["max"] - 4.0) <= 1.0, s2


def test_calibrated_and_relative_density():
    grey, shell, lumen, truth = multi_chamber_phantom(air=10.0, calcite=200.0)
    vs = 1.0
    # uncalibrated: no calcite standard
    u = foram.calibrate_shell_density(grey, shell, {"estimate_air": True}, erode=0)
    assert u["density_calibrated"] is False
    assert u["density_method"] == "uncalibrated_mean_grey"
    assert abs(u["density_index"] - 200.0) < 1e-6
    assert any("UNCALIBRATED" in n for n in u["notes"])
    # Iwasaki two-point
    cal = {"method": "iwasaki_ct", "air_grey": 10.0, "calcite_grey": 200.0,
           "calcite_density_g_cm3": 2.71}
    d = foram.calibrate_shell_density(grey, shell, cal, erode=0)
    assert d["density_calibrated"] is True
    assert abs(d["relative_ct_mean"] - 1000.0) < 1e-6
    assert abs(d["shell_density_g_cm3"] - 2.71) < 1e-6
    # linear standards, same two points
    lin = {"method": "linear", "grey_standards": [10.0, 200.0], "density_standards": [0.0, 2.71]}
    d2 = foram.calibrate_shell_density(grey, shell, lin, erode=0)
    assert abs(d2["shell_density_g_cm3"] - 2.71) < 1e-6
    # end-to-end analyse: mass and flag
    res = foram.analyse(grey, vs, min_size=8, close_radius=0, calibration=cal,
                        min_chamber_voxels=30, smooth_sigma=0.0, h_ratio=0.25)
    assert res.n_chambers == 3
    assert res.density_calibrated is True
    assert abs(res.shell_density_g_cm3 - 2.71) < 1e-3
    expected_mass = 2.71 * res.shell_volume_um3 * 1e-6
    assert abs(res.shell_mass_ug - expected_mass) < 1e-9
    assert res.test_bulk_density_g_cm3 is not None
    un = foram.analyse(grey, vs, min_size=8, close_radius=0, calibration=None,
                       chambers=False, wall_thickness_map=False)
    assert un.density_calibrated is False
    assert un.shell_mass_ug is None
    assert any("UNCALIBRATED" in n for n in un.notes)


def test_cli_foram_on_phantom(tmp_path):
    grey, shell, lumen, truth = multi_chamber_phantom()
    path = tmp_path / "ph.tif"
    save_volume(path, Volume(grey.astype(np.uint8), 1.0))
    out = tmp_path / "out"
    cal = {"air_grey": 10.0, "calcite_grey": 200.0, "calcite_density_g_cm3": 2.71}
    (tmp_path / "cal.yaml").write_text(yaml.safe_dump(cal))
    rc = main([
        "foram", "--input", str(path), "--voxel-size", "1.0", "--dtype", "uint8",
        "--output", str(out), "--name", "phantom",
        "--close-radius", "0", "--min-chamber-voxels", "30",
        "--h-ratio", "0.25", "--smooth-sigma", "0",
        "--calibration", str(tmp_path / "cal.yaml"),
    ])
    assert rc == 0
    res = json.loads((out / "results.json").read_text())["foram"]
    assert res["n_chambers"] == 3
    assert res["density_calibrated"] is True
    assert (out / "qc_chambers.png").exists()
    assert (out / "qc_thickness.png").exists()
    assert (out / "chamber_labels.npy").exists()
    assert (out / "wall_thickness_um.npy").exists()
