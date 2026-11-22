import numpy as np
import pytest

from morphoct import fibsem
from morphoct.io import Volume, save_volume


def _fei_xml(path, dx_nm=5.0, dy_nm=5.0, dz_nm=10.0):
    path.write_text(
        """<?xml version="1.0"?>
<Metadata>
  <PixelSize unit="nm">5.0</PixelSize>
  <PixelSizeX unit="nm">{dx}</PixelSizeX>
  <PixelSizeY unit="nm">{dy}</PixelSizeY>
  <SliceThickness unit="nm">{dz}</SliceThickness>
  <HighVoltage unit="V">5000</HighVoltage>
  <WorkingDistance unit="m">0.004</WorkingDistance>
</Metadata>
""".format(dx=dx_nm, dy=dy_nm, dz=dz_nm)
    )
    return path


def test_parse_fei_xml_voxel_size(tmp_path):
    p = _fei_xml(tmp_path / "meta.xml", dx_nm=5, dy_nm=5, dz_nm=10)
    meta = fibsem.parse_metadata_xml(p)
    vs = meta["voxel_size_um"]
    assert vs == pytest.approx((0.01, 0.005, 0.005))
    assert meta["high_voltage_v"] == 5000


def test_parse_text_sidecar(tmp_path):
    p = tmp_path / "info.txt"
    p.write_text("pixel size: 8 nm\nslice thickness: 16 nm\n")
    meta = fibsem.parse_metadata_text(p)
    assert meta["voxel_size_um"] == pytest.approx((0.016, 0.008, 0.008))


def test_load_stack_from_frames_and_xml(tmp_path):
    import tifffile

    _fei_xml(tmp_path / "Atlas.xml", 4, 4, 8)
    for i in range(5):
        tifffile.imwrite(tmp_path / f"slice_{i:03d}.tif", (np.full((16, 20), i + 10, np.uint8)))
    vol = fibsem.load_stack(tmp_path)
    assert vol.shape == (5, 16, 20)
    assert vol.voxel_size == pytest.approx((0.008, 0.004, 0.004))
    assert vol.data[3].mean() == 13


def test_register_recovers_known_shifts():
    rng = np.random.default_rng(0)
    base = ndi_blob((32, 32), rng)
    stack = np.zeros((6, 32, 32), np.float32)
    true = [(0, 0)]
    stack[0] = base
    accy = accx = 0
    for z in range(1, 6):
        accy += 1
        accx += 1 if z % 2 else 0
        stack[z] = np.roll(np.roll(base, accy, axis=0), accx, axis=1)
        true.append((accy, accx))
    vol = Volume(stack, (0.01, 0.005, 0.005))
    aligned, shifts = fibsem.register_slices(vol, upsample=1)
    # recovered shifts should be close to the known integer rolls
    for z, (ty, tx) in enumerate(true):
        assert abs(shifts[z, 0] - ty) <= 1.1
        assert abs(shifts[z, 1] - tx) <= 1.1
    # after alignment, later slices should match the first on the overlap
    a, b = aligned.data[0], aligned.data[-1]
    # correlation on the overlap should be high
    a = a - a.mean()
    b = b - b.mean()
    corr = float((a * b).mean() / (a.std() * b.std() + 1e-12))
    assert corr > 0.8


def ndi_blob(shape, rng):
    y, x = np.ogrid[:shape[0], :shape[1]]
    img = np.zeros(shape, np.float32)
    for _ in range(8):
        cy, cx = rng.uniform(6, shape[0] - 6), rng.uniform(6, shape[1] - 6)
        r = rng.uniform(2.5, 4.5)
        img += np.exp(-(((y - cy) ** 2 + (x - cx) ** 2) / (2 * r * r)))
    return img * 100 + 20


def test_curtaining_reduces_vertical_stripes():
    rng = np.random.default_rng(1)
    clean = np.stack([ndi_blob((40, 48), rng) for _ in range(4)])
    stripes = 25 * np.sin(np.arange(48) * 0.8)[None, None, :]
    dirty = clean + stripes
    vol = Volume(dirty, (0.01, 0.01, 0.01))
    out = fibsem.correct_curtaining(vol, width=2)
    # column-mean variation should drop
    dirty_col = dirty.mean(axis=(0, 1)).std()
    clean_col = out.data.mean(axis=(0, 1)).std()
    assert clean_col < 0.5 * dirty_col


def test_shading_flattens_ramp():
    base = np.ones((4, 32, 32), np.float32) * 50
    ramp = np.linspace(0, 40, 32)[None, None, :]
    vol = Volume(base + ramp, 1.0)
    out = fibsem.correct_shading(vol, sigma=6)
    assert out.data[0].std() < 0.4 * vol.data[0].std()


def test_crop_milled_region_finds_window():
    d = np.full((8, 40, 50), 10.0, np.float32)
    rng = np.random.default_rng(2)
    d[:, 8:30, 10:40] = 40 + rng.normal(0, 8, (8, 22, 30))
    # add slice-to-slice change inside the window only
    d[:, 8:30, 10:40] += np.arange(8)[:, None, None] * 3
    vol = Volume(d, 1.0)
    crop = fibsem.crop_milled_region(vol, pad=0)
    # the crop should be close to the true window
    assert crop.shape[1] <= 26
    assert crop.shape[2] <= 34
    assert crop.shape[1] >= 18
    assert crop.shape[2] >= 26


def test_resample_isotropic():
    d = np.zeros((4, 8, 8), np.float32)
    vol = Volume(d, (0.02, 0.01, 0.01))
    iso = fibsem.resample_isotropic(vol)
    assert iso.voxel_size == pytest.approx((0.01, 0.01, 0.01))
    assert iso.shape[0] == 8  # z doubled
    assert iso.shape[1:] == (8, 8)


def test_prepare_and_cli(tmp_path):
    import tifffile
    import yaml
    from morphoct.cli import main

    _fei_xml(tmp_path / "meta.xml", 10, 10, 10)
    rng = np.random.default_rng(3)
    # a small pore-like dark blob in a brighter matrix, 16^3, isotropic
    z, y, x = np.ogrid[:16, :16, :16]
    dark = ((z - 8) ** 2 + (y - 8) ** 2 + (x - 8) ** 2) < 16
    stack = np.where(dark, 40, 180).astype(np.uint8)
    stack = stack + rng.integers(0, 8, stack.shape, dtype=np.uint16)
    stack = stack.clip(0, 255).astype(np.uint8)
    for i in range(16):
        tifffile.imwrite(tmp_path / f"f_{i:03d}.tif", stack[i])
    out = tmp_path / "out"
    cfg = {
        "name": "fib_phantom",
        "output": str(out),
        "input": {"path": str(tmp_path)},
        "fibsem": {"register": False, "curtaining": False, "shading": False,
                   "crop_milled": False, "isotropic": False},
        "preprocess": {"denoise": "none"},
        "segment": {"method": "otsu", "pore": "dark"},
        "rve": False,
        "properties": {},
    }
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(cfg))
    assert main(["fibsem", str(tmp_path / "c.yaml")]) == 0
    import json
    res = json.loads((out / "results.json").read_text())
    assert 0.05 < res["morphometrics"]["porosity"] < 0.5

