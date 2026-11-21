"""Foraminifera module (Phase 3 starting point).

Implemented here: grey-scale shell segmentation, calcite volume, enclosed test
volume, a surface-area estimate, and simple density proxies. Chamber labelling,
ambient-occlusion interior separation and wall-thickness maps remain stubs.

Only generic, self-written helpers are included. No code is copied from ForamJ
(no licence file in its repository) or MFSE (GPL-3.0, MATLAB; compatible, but we
port the method, not the code).

Method sources (DOIs verified):

* Trend et al. 2026, ForamJ, Scientific Reports, doi:10.1038/s41598-026-43276-3:
  threshold the test, remove homogeneous infill, fill chambers, Euclidean distance
  transform (EDT) to separate septa from the outer wall, seeded 3D watershed of the
  inner chamber space, per-chamber volume, centroid and mean intensity.
* Siccha et al. 2026, MFSE, J. Open Res. Softw., doi:10.5334/jors.662, and Siccha et al.
  2026, Mar. Micropaleontol., doi:10.1016/j.marmicro.2026.102624: chamber lumen vs
  exterior by ambient occlusion (Titschack et al. 2018, doi:10.1016/j.cageo.2018.03.007),
  threshold chosen to maximise the volume-to-surface ratio; growth-stage ("in vivo")
  delimitation.
* Johnstone et al. 2010, Mar. Micropaleontol., doi:10.1016/j.marmicro.2010.07.004:
  CT-number based shell density and dissolution assessment, calibrated with standards.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field

import numpy as np
from scipy import ndimage as ndi


@dataclass
class ForamResult:
    shell_volume_um3: float
    filled_volume_um3: float | None = None
    convex_hull_volume_um3: float | None = None
    surface_area_um2: float | None = None
    test_porosity: float | None = None
    calcification_fraction: float | None = None
    shell_to_hull_fraction: float | None = None
    close_radius_voxels: int | None = None
    density_index: float | None = None
    threshold: float | None = None
    n_shell_voxels: int | None = None
    chamber_volumes_um3: list[float] = field(default_factory=list)
    chamber_centroids_um: list[tuple[float, float, float]] = field(default_factory=list)
    wall_thickness_um: dict = field(default_factory=dict)
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def shell_volume(shell_mask: np.ndarray, voxel_size: float) -> float:
    """Calcite volume (um^3) of a binary shell mask."""
    return float(np.count_nonzero(shell_mask)) * voxel_size**3


def filled_test_volume(shell_mask: np.ndarray, voxel_size: float, *, close_radius: int = 0) -> float:
    """Volume (um^3) of the test after filling enclosed chambers (binary fill holes).

    If ``close_radius`` is positive, a spherical morphological closing of that
    radius (voxels) is applied first so that foram apertures smaller than the
    closing diameter can be sealed. Without closing, chambers that open to the
    exterior through the aperture are not filled.
    """
    m = np.asarray(shell_mask).astype(bool)
    if close_radius:
        from skimage.morphology import ball

        try:
            from skimage.morphology import closing as _close
        except ImportError:
            from skimage.morphology import binary_closing as _close
        m = _close(m, ball(close_radius))
    filled = ndi.binary_fill_holes(m)
    return float(np.count_nonzero(filled)) * voxel_size**3


def convex_hull_volume(shell_mask: np.ndarray, voxel_size: float) -> float:
    """Volume (um^3) of the 3D convex hull of the shell.

    Computed with Qhull on the corners of surface voxels, so a solid block of
    n x n x n voxels has hull volume n^3 (not the smaller volume of voxel centres).
    """
    from itertools import product

    from scipy.spatial import ConvexHull, QhullError

    m = np.asarray(shell_mask).astype(bool)
    if not np.any(m):
        return 0.0
    eroded = ndi.binary_erosion(m)
    surf = m & ~eroded
    idx = np.argwhere(surf if np.any(surf) else m)
    if idx.shape[0] < 2:
        return float(np.count_nonzero(m)) * voxel_size**3
    offs = np.array(list(product((-0.5, 0.5), repeat=3)))
    pts = (idx[:, None, :] + offs[None, :, :]).reshape(-1, 3)
    try:
        hull = ConvexHull(pts)
    except QhullError:
        return float("nan")
    return float(hull.volume) * voxel_size**3


def shell_surface_area(shell_mask: np.ndarray, voxel_size: float, *, method: str = "marching_cubes") -> float:
    """Shell surface area (um^2), including the inner chamber walls of a closed test.

    ``marching_cubes`` uses scikit-image (Lewiner et al. 2003). ``voxel_faces`` counts
    exposed voxel faces (an upper bound used as a sanity check).
    """
    m = np.asarray(shell_mask).astype(bool)
    if method == "voxel_faces":
        area = 0
        for ax in range(3):
            diff = np.diff(m.astype(np.int8), axis=ax)
            area += int(np.count_nonzero(diff))
        return float(area) * voxel_size**2
    if method != "marching_cubes":
        raise ValueError("method: 'marching_cubes' or 'voxel_faces'")
    from skimage import measure

    if m.ndim != 3 or not np.any(m) or np.all(m):
        return 0.0
    try:
        verts, faces, *_ = measure.marching_cubes(m.astype(np.uint8), level=0.5, spacing=(voxel_size,) * 3)
    except (ValueError, RuntimeError):
        return 0.0
    return float(measure.mesh_surface_area(verts, faces))


def linear_density_calibration(grey_standards, density_standards):
    """Fit grey = a * density + b on scanned standards (e.g. air, water, calcite crystal)
    and return a function grey -> density (g/cm^3). Needs at least two standards scanned
    with the same settings as the specimens."""
    g = np.asarray(grey_standards, float)
    d = np.asarray(density_standards, float)
    if g.size < 2:
        raise ValueError("at least two standards")
    a, b = np.polyfit(d, g, 1)
    return lambda grey: (np.asarray(grey, float) - b) / a


def density_index(grey: np.ndarray, shell_mask: np.ndarray, calibration=None, *, erode: int = 1) -> float:
    """Mean grey value (or calibrated density) of the shell interior, eroded by ``erode``
    voxels to limit partial-volume bias at the surfaces."""
    m = shell_mask.astype(bool)
    if erode:
        m = ndi.binary_erosion(m, iterations=erode)
    vals = grey[m]
    if vals.size == 0:
        raise ValueError("shell too thin for the erosion")
    v = float(vals.mean())
    return float(calibration(v)) if calibration else v


def segment_test(grey: np.ndarray, *, method: str = "otsu", min_size: int = 64,
                 largest_only: bool = True) -> tuple[np.ndarray, float, dict]:
    """Segment the calcite shell from a grey microCT volume.

    Baseline (not the full ForamJ or MFSE pipeline): Otsu threshold with the shell
    treated as the bright phase, optional removal of small objects, and optionally
    keeping only the largest connected component (the specimen). Returns the uint8
    shell mask (1 = calcite), the threshold, and a small info dict.

    Infill removal, aperture-aware interior separation and chamber watershed are
    not done here.
    """
    import inspect

    from skimage.filters import threshold_otsu
    from skimage.morphology import remove_small_objects

    g = np.asarray(grey)
    if g.ndim != 3:
        raise ValueError(f"expected a 3D volume, got shape {g.shape}")
    if method != "otsu":
        raise ValueError("only method='otsu' is implemented")
    flat = g.ravel()
    if flat.size > 2_000_000:
        rng = np.random.default_rng(0)
        sample = flat[rng.integers(0, flat.size, 2_000_000)]
    else:
        sample = flat
    thr = float(threshold_otsu(sample))
    shell = g > thr
    if min_size:
        params = inspect.signature(remove_small_objects).parameters
        if "max_size" in params:
            # skimage 0.26+: max_size removes objects smaller than or equal to the value
            shell = remove_small_objects(shell, max_size=max(min_size - 1, 0))
        else:
            shell = remove_small_objects(shell, min_size=min_size)
    info = {"method": method, "threshold": thr, "min_size": min_size, "largest_only": largest_only}
    if largest_only and np.any(shell):
        labels, n = ndi.label(shell)
        if n:
            counts = np.bincount(labels.ravel())
            counts[0] = 0
            keep = int(np.argmax(counts))
            shell = labels == keep
            info["n_components_before"] = int(n)
    mask = shell.astype(np.uint8)
    info["n_shell_voxels"] = int(mask.sum())
    return mask, thr, info


def separate_interior(shell_mask: np.ndarray, *, method: str = "ambient_occlusion", n_rays: int = 64):
    """Planned: chamber lumen vs exterior by ambient occlusion with the threshold that
    maximises volume/surface (MFSE method, Titschack et al. 2018), reimplemented."""
    raise NotImplementedError("Phase 3: ambient-occlusion interior separation")


def segment_chambers(interior_mask: np.ndarray, seeds=None, *, smooth_sigma: float = 1.0):
    """Planned: EDT + seeded 3D watershed of the chamber space into chambers (ForamJ
    method from the paper), returning int labels ordered by growth sequence."""
    raise NotImplementedError("Phase 3: chamber watershed")


def wall_thickness(shell_mask: np.ndarray, voxel_size: float, *, split_septa: bool = True):
    """Planned: local thickness of outer wall and septa (PoreSpy local thickness, checked
    against BoneJ2)."""
    raise NotImplementedError("Phase 3: wall and septa thickness")


def analyse(grey: np.ndarray, voxel_size: float, **kw) -> ForamResult:
    """End-to-end baseline: segment the shell, then volume, enclosed test volume,
    surface area and density proxies. Chamber-level metrics remain future work.

    Keyword arguments are passed to :func:`segment_test`. ``density_erode`` (default 1)
    is used by :func:`density_index`.
    """
    erode = kw.pop("density_erode", 1)
    close_radius = int(kw.pop("close_radius", 8))
    shell, thr, info = segment_test(grey, **kw)
    vs = float(voxel_size)
    sv = shell_volume(shell, vs)
    fv = filled_test_volume(shell, vs, close_radius=close_radius)
    hv = convex_hull_volume(shell, vs)
    sa = shell_surface_area(shell, vs, method="marching_cubes")
    test_phi = float((fv - sv) / fv) if fv > 0 else float("nan")
    calc = float(sv / fv) if fv > 0 else float("nan")
    hull_frac = float(sv / hv) if hv > 0 else float("nan")
    notes = [
        "Baseline Otsu shell segmentation; chambers are not labelled.",
        f"Filled test volume applies a morphological closing of radius {close_radius} voxels and then fills holes, to seal apertures smaller than the closing diameter.",
        "test_porosity is enclosed void volume over filled test volume; calcification_fraction is shell over filled test volume; shell_to_hull_fraction is shell over the convex hull (a parameter-light packing proxy).",
    ]
    dens = None
    try:
        dens = density_index(np.asarray(grey), shell, erode=erode)
    except ValueError as e:
        notes.append(f"density_index skipped: {e}")
    return ForamResult(
        shell_volume_um3=sv,
        filled_volume_um3=fv,
        surface_area_um2=sa,
        test_porosity=test_phi,
        calcification_fraction=calc,
        shell_to_hull_fraction=hull_frac,
        density_index=dens,
        threshold=thr,
        n_shell_voxels=info.get("n_shell_voxels"),
        convex_hull_volume_um3=hv,
        close_radius_voxels=close_radius,
        notes=notes,
    )
