"""Cropping, denoising, and mirroring for periodic solvers."""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from morphoct.io import Volume


def crop(vol: Volume, start=(0, 0, 0), size=None, *, center=False) -> Volume:
    """Sub-volume of ``size`` voxels (z, y, x) from ``start``, or centred when ``center``."""
    shape = np.array(vol.shape)
    size = shape if size is None else np.minimum(np.array(size, dtype=int), shape)
    start = (shape - size) // 2 if center else np.array(start, dtype=int)
    if np.any(start < 0) or np.any(start + size > shape):
        raise ValueError(f"crop {start}+{size} outside volume {tuple(shape)}")
    sl = tuple(slice(int(a), int(a + b)) for a, b in zip(start, size))
    return vol.with_data(np.ascontiguousarray(vol.data[sl]), crop_start=start.tolist(), crop_size=size.tolist())


def denoise(vol: Volume, method: str = "median", **kw) -> Volume:
    """``median`` (size, default 3), ``gaussian`` (sigma, default 1), or ``nlm`` (non-local
    means of scikit-image; ``h`` in units of the estimated noise sigma, default 0.8;
    ``patch_size`` 3, ``patch_distance`` 5). NLM on large volumes is slow: crop first."""
    d = vol.data
    if method in (None, "none"):
        return vol
    if method == "median":
        out = ndi.median_filter(d, size=int(kw.get("size", 3)))
    elif method == "gaussian":
        out = ndi.gaussian_filter(d.astype(np.float32), float(kw.get("sigma", 1.0)))
    elif method == "nlm":
        from skimage.restoration import denoise_nl_means, estimate_sigma

        f = d.astype(np.float32)
        lo, hi = np.percentile(f, [0.1, 99.9])
        f = np.clip((f - lo) / max(hi - lo, 1e-12), 0, 1)
        sigma = float(np.mean(estimate_sigma(f[f.shape[0] // 2], channel_axis=None)))
        out = denoise_nl_means(f, h=float(kw.get("h", 0.8)) * sigma, sigma=sigma, fast_mode=True,
                               patch_size=int(kw.get("patch_size", 3)),
                               patch_distance=int(kw.get("patch_distance", 5)), channel_axis=None)
        out = (out * (hi - lo) + lo).astype(np.float32)
    else:
        raise ValueError(f"unknown denoise method {method!r}")
    return vol.with_data(out, denoise=method)


def mirror(binary: np.ndarray, axes=(0,)) -> np.ndarray:
    """Concatenate the volume with its reflection along each axis in ``axes``.

    A non-periodic sample becomes periodic along those axes without the artificial
    junction a plain wrap would create, so periodic FFT and lattice Boltzmann solvers can
    be used. The size doubles along each mirrored axis; along a mirrored flow axis the
    apparent permeability and conductivity of a symmetric doubled sample equal those of
    the series combination of the original and its reflection.
    """
    out = np.asarray(binary)
    for a in axes:
        out = np.concatenate([out, np.flip(out, axis=a)], axis=a)
    return np.ascontiguousarray(out)


def downsample(vol: Volume, factor: int = 2) -> Volume:
    """Block-average by an integer factor (grey) or majority vote (binary)."""
    f = int(factor)
    d = vol.data
    nz, ny, nx = (s // f * f for s in d.shape)
    d = d[:nz, :ny, :nx].reshape(nz // f, f, ny // f, f, nx // f, f).mean(axis=(1, 3, 5))
    if vol.data.dtype == bool or set(np.unique(vol.data[:2]).tolist()) <= {0, 1}:
        d = (d >= 0.5).astype(vol.data.dtype)
    return Volume(d, tuple(v * f for v in vol.voxel_size), {**vol.meta, "downsample": f})
