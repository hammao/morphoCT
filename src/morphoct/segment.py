"""Baseline segmentation: global thresholds and marker-based watershed, with QC plots.

These are deliberately simple, reproducible baselines (Otsu 1979; multi-Otsu; Andrä et
al. 2013 Part I show segmentation is the largest single source of spread in digital rock
results). A threshold-sensitivity sweep is reported so that every downstream number can
be quoted with its segmentation uncertainty.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage as ndi

from morphoct.io import Volume


def _sample(d: np.ndarray, n: int = 2_000_000) -> np.ndarray:
    flat = d.ravel()
    if flat.size <= n:
        return flat
    rng = np.random.default_rng(0)
    return flat[rng.integers(0, flat.size, n)]


def otsu(vol: Volume, *, pore="dark") -> tuple[np.ndarray, float]:
    """Binary pore mask (uint8, 1 = pore) and the threshold. ``pore``: ``dark`` (X-ray CT
    of rock: pores attenuate least) or ``bright``."""
    from skimage.filters import threshold_otsu

    t = float(threshold_otsu(_sample(vol.data)))
    mask = vol.data <= t if pore == "dark" else vol.data > t
    return mask.astype(np.uint8), t


def multi_otsu(vol: Volume, classes: int = 3) -> tuple[np.ndarray, list[float]]:
    """Labels 0..classes-1 in order of grey value (e.g. pore, matrix/microporosity, grain)."""
    from skimage.filters import threshold_multiotsu

    ts = threshold_multiotsu(_sample(vol.data), classes=classes)
    return np.digitize(vol.data, ts).astype(np.uint8), [float(t) for t in ts]


def clean(mask: np.ndarray, min_size: int = 0) -> np.ndarray:
    """Remove pore clusters and solid specks smaller than ``min_size`` voxels."""
    if min_size <= 0:
        return mask
    from skimage.morphology import remove_small_objects

    m = remove_small_objects(mask.astype(bool), min_size=min_size)
    s = remove_small_objects(~m, min_size=min_size)
    return (~s).astype(np.uint8)


def watershed_grains(mask_solid: np.ndarray, *, min_distance: int = 3, footprint: int = 7) -> np.ndarray:
    """Marker-based watershed of the solid phase into grains (distance-transform peaks as
    markers). Returns int32 grain labels (0 = pore)."""
    from skimage.segmentation import watershed

    solid = mask_solid.astype(bool)
    dt = ndi.distance_transform_edt(solid)
    dt_s = ndi.gaussian_filter(dt, 0.5)
    peaks = (dt_s == ndi.maximum_filter(dt_s, size=footprint)) & (dt_s >= min_distance)
    markers, _ = ndi.label(peaks)
    return watershed(-dt_s, markers, mask=solid).astype(np.int32)


def threshold_sensitivity(vol: Volume, threshold: float, rel=(-0.05, -0.025, 0.0, 0.025, 0.05), *, pore="dark"):
    """Porosity at thresholds shifted by fractions of the grey range: the segmentation
    uncertainty band."""
    s = _sample(vol.data)
    lo, hi = np.percentile(s, [0.5, 99.5])
    rows = []
    for r in rel:
        t = threshold + r * (hi - lo)
        phi = float(np.mean(s <= t) if pore == "dark" else np.mean(s > t))
        rows.append({"shift": r, "threshold": float(t), "porosity": phi})
    return rows


def dice(a: np.ndarray, b: np.ndarray) -> float:
    a, b = a.astype(bool), b.astype(bool)
    return float(2 * np.logical_and(a, b).sum() / max(a.sum() + b.sum(), 1))


def qc_figure(gray: np.ndarray | None, mask: np.ndarray, path, *, thresholds=(), title="Segmentation QC"):
    """Three orthogonal mid-slices with the pore mask overlaid, plus the grey histogram
    with the thresholds."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = 4 if gray is not None else 3
    fig, axs = plt.subplots(1, n, figsize=(4 * n, 4))
    c = [s // 2 for s in mask.shape]
    sl = [(np.s_[c[0], :, :], "z mid"), (np.s_[:, c[1], :], "y mid"), (np.s_[:, :, c[2]], "x mid")]
    for ax, (s, lab) in zip(axs, sl):
        base = gray[s] if gray is not None else mask[s]
        ax.imshow(base, cmap="gray", interpolation="nearest")
        ax.imshow(np.ma.masked_where(mask[s] == 0, mask[s]), cmap="autumn", alpha=0.35, interpolation="nearest")
        ax.set_title(lab)
        ax.axis("off")
    if gray is not None:
        ax = axs[-1]
        ax.hist(_sample(gray, 500_000), bins=256, color="0.3")
        for t in thresholds:
            ax.axvline(t, color="r")
        ax.set_title("grey histogram")
        ax.set_yscale("log")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path
