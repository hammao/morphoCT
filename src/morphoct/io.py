"""Reading and writing microCT volumes.

Every reader returns a :class:`Volume`: a 3D numpy array indexed ``[z, y, x]`` (slice
order, as stacks are stored) with an isotropic or anisotropic voxel size in micrometres.
:func:`to_voxel_image` converts it to a morphoGen ``VoxelImage`` (indexed ``[x, y, z]``).
"""

from __future__ import annotations

import glob
import json
import os
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


@dataclass
class Volume:
    data: np.ndarray
    voxel_size: tuple[float, float, float] = (1.0, 1.0, 1.0)  # (dz, dy, dx) in micrometres
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.data = np.asarray(self.data)
        if self.data.ndim != 3:
            raise ValueError(f"Volume: expected a 3D array, got shape {self.data.shape}")
        vs = self.voxel_size
        if np.isscalar(vs):
            vs = (float(vs),) * 3
        self.voxel_size = tuple(float(v) for v in vs)
        if any(v <= 0 for v in self.voxel_size):
            raise ValueError("voxel_size must be positive")

    @property
    def shape(self) -> tuple[int, int, int]:
        return tuple(self.data.shape)

    @property
    def isotropic_voxel(self) -> float:
        if not np.allclose(self.voxel_size, self.voxel_size[0], rtol=1e-6):
            raise ValueError(f"anisotropic voxels {self.voxel_size}: resample first")
        return self.voxel_size[0]

    def with_data(self, data: np.ndarray, **meta) -> "Volume":
        return Volume(data, self.voxel_size, {**self.meta, **meta})


_ITK_EXT = (".nii", ".nii.gz", ".nrrd", ".nhdr", ".mha", ".mhd")


def _detect(path: str) -> str:
    p = path.lower()
    if os.path.isdir(path) or any(c in path for c in "*?["):
        return "tiff_dir"
    if p.endswith((".tif", ".tiff")):
        return "tiff"
    if p.endswith(".npy"):
        return "npy"
    if p.endswith(".npz"):
        return "npz"
    if p.endswith((".raw", ".raw.gz", ".bin", ".vol", ".dat")):
        return "raw"
    if p.endswith(_ITK_EXT):
        return "itk"
    raise ValueError(f"cannot guess the format of {path!r}; pass format=")


def read_raw(path, shape, dtype="uint8", byteorder="little", offset=0, order="zyx") -> np.ndarray:
    """A headerless binary volume. ``shape`` is given in the order of ``order``: ``zyx``
    means x varies fastest on disk (ImageJ, Avizo, most scanners)."""
    dt = np.dtype(dtype).newbyteorder("<" if byteorder == "little" else ">")
    shape = tuple(int(s) for s in shape)
    if str(path).endswith(".gz"):
        import gzip

        with gzip.open(path, "rb") as f:
            f.read(offset)
            arr = np.frombuffer(f.read(), dtype=dt)
    else:
        arr = np.fromfile(path, dtype=dt, offset=offset)
    n = int(np.prod(shape))
    if arr.size != n:
        raise ValueError(f"{path}: {arr.size} values on disk, shape {shape} needs {n}")
    arr = arr.reshape(shape)
    if order == "xyz":
        arr = arr.transpose(2, 1, 0)
    elif order != "zyx":
        raise ValueError("order: 'zyx' or 'xyz'")
    return np.ascontiguousarray(arr.astype(dt.newbyteorder("=")))


def read_tiff_stack(path) -> tuple[np.ndarray, dict]:
    """A multi-page TIFF, or a directory/glob of single-slice TIFFs sorted by name."""
    import tifffile

    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, "*.tif")) + glob.glob(os.path.join(path, "*.tiff")))
    elif any(c in str(path) for c in "*?["):
        files = sorted(glob.glob(str(path)))
    else:
        files = [str(path)]
    if not files:
        raise FileNotFoundError(path)
    if len(files) == 1:
        with tifffile.TiffFile(files[0]) as tif:
            arr = tif.asarray()
            info = {"imagej": tif.imagej_metadata or {}}
        if arr.ndim == 2:
            arr = arr[None]
    else:
        arr = np.stack([tifffile.imread(f) for f in files])
        info = {"n_files": len(files)}
    return arr, info


def load_volume(path, *, format: str | None = None, voxel_size=None, shape=None,
                dtype="uint8", byteorder="little", offset=0, order="zyx", crop=None) -> Volume:
    """Load a microCT volume.

    ``format``: ``tiff``, ``tiff_dir``, ``raw``, ``npy``, ``npz`` or ``itk`` (NIfTI, NRRD,
    MetaImage via SimpleITK); guessed from the extension when omitted. RAW needs ``shape``
    and ``dtype``. ``voxel_size`` (micrometres) overrides what the file says. ``crop`` is
    ``[[z0, z1], [y0, y1], [x0, x1]]`` applied right after reading.
    """
    path = str(path)
    fmt = format or _detect(path)
    meta = {"source": path, "format": fmt}
    vs = None
    if fmt in ("tiff", "tiff_dir"):
        data, info = read_tiff_stack(path)
        sp = (info.get("imagej") or {}).get("spacing")
        if sp:
            meta["tiff_spacing"] = sp
    elif fmt == "raw":
        if shape is None:
            raise ValueError("RAW volumes need shape=(nz, ny, nx) and dtype=")
        data = read_raw(path, shape, dtype, byteorder, offset, order)
    elif fmt == "npy":
        data = np.load(path)
    elif fmt == "npz":
        z = np.load(path)
        data = z[z.files[0]]
        if "voxel_size" in z.files:
            vs = tuple(z["voxel_size"])
    elif fmt == "itk":
        try:
            import SimpleITK as sitk
        except ImportError as e:  # pragma: no cover
            raise ImportError("NIfTI/NRRD/MHA need SimpleITK: pip install 'morphoct[itk]'") from e
        img = sitk.ReadImage(path)
        data = sitk.GetArrayFromImage(img)  # z, y, x
        sx, sy, sz = img.GetSpacing()
        vs = (sz, sy, sx)
        meta["itk_spacing_units"] = "as stored (usually mm for NIfTI): check"
    else:
        raise ValueError(f"unknown format {fmt!r}")
    if crop is not None:
        (z0, z1), (y0, y1), (x0, x1) = crop
        data = data[z0:z1, y0:y1, x0:x1]
    if voxel_size is not None:
        vs = voxel_size
    if vs is None:
        vs = (1.0, 1.0, 1.0)
        meta["voxel_size_warning"] = "voxel size unknown, set voxel_size (micrometres)"
    return Volume(np.ascontiguousarray(data), vs, meta)


def save_volume(path, vol: Volume) -> Path:
    """Write ``.npz`` (with voxel size), ``.npy``, ``.raw`` (+ JSON sidecar) or ``.tif``."""
    path = Path(path)
    s = str(path).lower()
    if s.endswith(".npz"):
        np.savez_compressed(path, data=vol.data, voxel_size=np.array(vol.voxel_size))
    elif s.endswith(".npy"):
        np.save(path, vol.data)
    elif s.endswith(".raw"):
        vol.data.tofile(path)
        side = {"shape_zyx": list(vol.shape), "dtype": str(vol.data.dtype),
                "byteorder": "little", "voxel_size_um": list(vol.voxel_size)}
        path.with_suffix(".json").write_text(json.dumps(side, indent=2))
    elif s.endswith((".tif", ".tiff")):
        import tifffile

        d = vol.data.astype(np.uint8) if vol.data.dtype == bool else vol.data
        dz, dy, dx = vol.voxel_size
        tifffile.imwrite(path, d, imagej=True, resolution=(1 / dx, 1 / dy),
                         metadata={"spacing": dz, "unit": "um", "axes": "ZYX"})
    else:
        raise ValueError(f"unsupported output {path}")
    return path


def to_voxel_image(binary: np.ndarray, voxel_size, *, periodic=False, name="phase"):
    """A morphoGen VoxelImage from a ``[z, y, x]`` label array (phase 1 = pore)."""
    from morphogen.image import VoxelImage

    vs = voxel_size if not np.isscalar(voxel_size) else (voxel_size,) * 3
    dz, dy, dx = vs
    return VoxelImage(np.ascontiguousarray(np.asarray(binary).transpose(2, 1, 0)),
                      spacing=(dx, dy, dz), periodic=periodic, name=name)
