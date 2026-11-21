import numpy as np

from morphoct.io import Volume, load_volume, save_volume, to_voxel_image


def test_raw_roundtrip_orders(tmp_path):
    a = np.arange(2 * 3 * 4, dtype=np.uint16).reshape(2, 3, 4)
    a.astype("<u2").tofile(tmp_path / "v.raw")
    v = load_volume(tmp_path / "v.raw", shape=(2, 3, 4), dtype="uint16", voxel_size=2.0)
    assert np.array_equal(v.data, a) and v.voxel_size == (2.0, 2.0, 2.0)
    # the same bytes declared as x-fastest-last (xyz) come back transposed
    w = load_volume(tmp_path / "v.raw", shape=(2, 3, 4), dtype="uint16", order="xyz")
    assert w.shape == (4, 3, 2)


def test_raw_big_endian_and_gz(tmp_path):
    import gzip

    a = np.arange(24, dtype=">i2").reshape(2, 3, 4)
    with gzip.open(tmp_path / "v.raw.gz", "wb") as f:
        f.write(a.tobytes())
    v = load_volume(tmp_path / "v.raw.gz", shape=(2, 3, 4), dtype="int16", byteorder="big")
    assert np.array_equal(v.data, a.astype(np.int16))


def test_tiff_npz_npy(tmp_path):
    a = (np.random.default_rng(0).random((5, 6, 7)) * 255).astype(np.uint8)
    vol = Volume(a, 3.5)
    save_volume(tmp_path / "a.tif", vol)
    t = load_volume(tmp_path / "a.tif", voxel_size=3.5)
    assert np.array_equal(t.data, a)
    save_volume(tmp_path / "a.npz", vol)
    assert load_volume(tmp_path / "a.npz").voxel_size == (3.5, 3.5, 3.5)
    save_volume(tmp_path / "a.npy", vol)
    assert np.array_equal(load_volume(tmp_path / "a.npy").data, a)


def test_tiff_directory(tmp_path):
    import tifffile

    a = (np.arange(4 * 5 * 6) % 251).astype(np.uint8).reshape(4, 5, 6)
    d = tmp_path / "stack"
    d.mkdir()
    for i in range(4):
        tifffile.imwrite(d / f"s_{i:04d}.tif", a[i])
    assert np.array_equal(load_volume(d).data, a)


def test_itk_nrrd(tmp_path):
    import pytest

    sitk = pytest.importorskip("SimpleITK")
    a = np.random.default_rng(1).integers(0, 255, (4, 5, 6)).astype(np.uint8)
    img = sitk.GetImageFromArray(a)
    img.SetSpacing((1.0, 2.0, 3.0))
    sitk.WriteImage(img, str(tmp_path / "a.nrrd"))
    v = load_volume(tmp_path / "a.nrrd")
    assert np.array_equal(v.data, a) and v.voxel_size == (3.0, 2.0, 1.0)


def test_to_voxel_image_axes():
    a = np.zeros((2, 3, 4), np.uint8)
    img = to_voxel_image(a, (1.0, 2.0, 3.0))
    assert img.data.shape == (4, 3, 2) and tuple(img.spacing) == (3.0, 2.0, 1.0)
