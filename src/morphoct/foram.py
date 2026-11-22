"""Foraminifera module: shell segmentation, chambers, wall thickness and density.

Implemented from published method descriptions, with self-written code. No code is
copied from ForamJ (no licence file in its repository) or MFSE (GPL-3.0, MATLAB;
compatible, but we port the method, not the code).

Method sources (DOIs verified):

* Speijer, Van Loo, Masschaele, Vlassenbroeck, Cnudde and Jacobs 2008, Geosphere,
  doi:10.1130/GES00176.1: HRXCT chamber segmentation and volumes; chambers are
  closed at the foramen; ontogenetic volume increase after the proloculus.
* Briguglio, Hohenegger and Less 2013, J. Foraminifer. Res., doi:10.2113/gsjfr.43.1.72
  (see also Briguglio, Metscher and Hohenegger 2011, doi:10.3906/yer-0910-44):
  three-dimensional chamber-lumen volumes from microCT.
* Trend et al. 2026, ForamJ, Scientific Reports, doi:10.1038/s41598-026-43276-3:
  threshold the test, fill chambers, smooth the inner volume, Euclidean distance
  transform (EDT), seed the chamber centres, 3D watershed of the lumen, then
  per-chamber volume, centroid and mean intensity. BoneJ-style local thickness
  of the shell.
* Siccha et al. 2026, MFSE, J. Open Res. Softw., doi:10.5334/jors.662, and
  Siccha et al. 2026, Mar. Micropaleontol., doi:10.1016/j.marmicro.2026.102624:
  chamber lumen vs exterior (ambient occlusion, Titschack et al. 2018,
  doi:10.1016/j.cageo.2018.03.007, is the MFSE default; morphoCT uses
  morphological closing plus hole filling for the same split unless AO is requested).
* Hildebrand and Ruegsegger 1997, J. Microsc., doi:10.1046/j.1365-2818.1997.1340694.x:
  local thickness as the diameter of the largest inscribed sphere.
* Iwasaki et al. 2015, Paleoceanography, doi:10.1002/2014PA002639: relative CT
  number with air = 0 and a calcite standard = 1000:
  CT = 1000 * (mu - mu_air) / (mu_calcite - mu_air).
* Johnstone, Schulz, Barker and Elderfield 2010, Mar. Micropaleontol.,
  doi:10.1016/j.marmicro.2010.07.004: CT imaging of test calcite and dissolution;
  greyscale can mark preservation but is easily biased by infill, which is why they
  also defined the empirical XDX index. morphoCT uses the quantitative greyscale
  path, with an explicit uncalibrated flag when no standard is given.
* Todd, Schmidt, Robinson and De Schepper 2019, Paleoceanography and Paleoclimatology,
  doi:10.1029/2019PA003738: microCT calcite volume and test thickness in planktic
  foraminifera.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

CALCITE_DENSITY_G_CM3 = 2.71  # end-member calcite, used when a standard is scaled to 1000


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
    density_calibrated: bool = False
    density_method: str | None = None
    shell_density_g_cm3: float | None = None
    shell_mass_ug: float | None = None
    test_bulk_density_g_cm3: float | None = None
    relative_ct_mean: float | None = None
    air_grey: float | None = None
    calcite_grey: float | None = None
    threshold: float | None = None
    n_shell_voxels: int | None = None
    n_chambers: int | None = None
    chamber_volumes_um3: list[float] = field(default_factory=list)
    chamber_equivalent_radius_um: list[float] = field(default_factory=list)
    chamber_centroids_um: list[tuple[float, float, float]] = field(default_factory=list)
    chamber_mean_grey: list[float] = field(default_factory=list)
    chamber_order: str | None = None
    wall_thickness_um: dict = field(default_factory=dict)
    chamber_wall_thickness_um: list[dict] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return asdict(self)


def shell_volume(shell_mask: np.ndarray, voxel_size: float) -> float:
    """Calcite volume (um^3) of a binary shell mask."""
    return float(np.count_nonzero(shell_mask)) * voxel_size**3


def _close_mask(mask: np.ndarray, radius: int) -> np.ndarray:
    """Morphological closing by a Euclidean ball of ``radius`` voxels.

    Implemented with two Euclidean distance transforms so large radii stay
    cheap on ~200^3 volumes. The inequalities (dilate: dist <= R, erode:
    dist > R) match scikit-image ``closing`` with ``ball(R)``.
    """
    m = np.asarray(mask).astype(bool)
    if radius <= 0:
        return m
    dil = ndi.distance_transform_edt(~m) <= radius
    return ndi.distance_transform_edt(dil) > radius


def _border_touching(mask: np.ndarray) -> np.ndarray:
    """Connected components of ``mask`` that touch any face of the volume."""
    m = np.asarray(mask).astype(bool)
    labels, n = ndi.label(m)
    if n == 0:
        return np.zeros(m.shape, bool)
    border = np.zeros(m.shape, bool)
    border[0] = border[-1] = True
    border[:, 0] = border[:, -1] = True
    border[:, :, 0] = border[:, :, -1] = True
    keep = np.unique(labels[border])
    keep = keep[keep != 0]
    return np.isin(labels, keep) if keep.size else np.zeros(m.shape, bool)


def exterior_erosion_interior(shell_mask: np.ndarray, radius: int) -> dict:
    """Enclosed lumen by cutting narrow apertures from the outside.

    The background (everything that is not shell) is eroded by ``radius``
    voxels. Components that still touch the volume border are the true
    exterior; they are dilated back by the same radius, stopping at the shell.
    Narrow tunnels (apertures and foramina whose half-width is <= radius) are
    cut and do not refill into the chambers, while wide bays such as a deep
    umbilicus that remain connected to the border stay exterior. This avoids
    the wrapping artefact of a large morphological closing, which can swallow
    exterior indentations into the test envelope.

    Small chambers whose radius is less than ``radius`` are still kept: they
    vanish during the erosion but are never marked as exterior, so they return
    in ``void & ~exterior``.
    """
    shell = np.asarray(shell_mask).astype(bool)
    void = ~shell
    if radius <= 0:
        filled = ndi.binary_fill_holes(shell)
        interior = filled & void
        return {"shell": shell, "closed": shell, "filled": filled, "interior": interior,
                "method": "fill_holes"}
    # Tightly cropped volumes may have less than ``radius`` voxels of
    # background at the faces, which would pinch the exterior away from the
    # border. Pad with background so the true exterior stays connected.
    pad = int(radius) + 2
    shell_p = np.pad(shell, pad, constant_values=False)
    void_p = ~shell_p
    core = ndi.distance_transform_edt(void_p) > radius
    exterior_core = _border_touching(core)
    exterior_p = ndi.distance_transform_edt(~exterior_core) <= radius
    exterior_p &= void_p
    sl = (slice(pad, -pad),) * shell.ndim
    exterior = exterior_p[sl]
    interior = void & ~exterior
    filled = interior | shell
    return {"shell": shell, "closed": shell, "filled": filled, "interior": interior,
            "method": "exterior_erosion", "exterior": exterior}



def _convex_hull_mask(shell: np.ndarray, *, step: int = 2) -> np.ndarray:
    """Rasterise the convex hull of the shell on a downsampled grid and upsample."""
    from scipy.spatial import ConvexHull, QhullError

    m = np.asarray(shell).astype(bool)
    eroded = ndi.binary_erosion(m)
    surf = m & ~eroded
    pts = np.argwhere(surf if np.any(surf) else m)
    if pts.shape[0] < 4:
        return m.copy()
    try:
        hull = ConvexHull(pts[:: max(int(step), 1)].astype(float))
    except QhullError:
        return m.copy()
    eq = hull.equations.astype(np.float32)
    grid = np.mgrid[0:m.shape[0]:step, 0:m.shape[1]:step, 0:m.shape[2]:step]
    coords = np.column_stack([g.ravel() for g in grid]).astype(np.float32)
    inside = np.ones(len(coords), bool)
    for start in range(0, len(eq), 64):
        E = eq[start:start + 64]
        inside &= (coords @ E[:, :3].T + E[:, 3] <= (step + 0.5)).all(axis=1)
        if not inside.any():
            break
    small = inside.reshape(grid[0].shape)
    out = np.repeat(np.repeat(np.repeat(small, step, 0), step, 1), step, 2)
    return out[: m.shape[0], : m.shape[1], : m.shape[2]]


def _enclosure_ratio(component: np.ndarray, shell: np.ndarray) -> float:
    """Shell-contact over void-contact of a dilated component (higher is more enclosed)."""
    dil = ndi.binary_dilation(component)
    shell_c = int((dil & shell).sum())
    void_c = int((dil & ~shell & ~component).sum())
    return shell_c / max(void_c, 1)


def hull_enclosure_interior(shell_mask: np.ndarray, radius: int = 8, *,
                            h_ratio: float = 0.2, min_chamber_voxels: int = 80,
                            min_ratio: float = 2.0) -> dict:
    """Enclosed lumen: tight exterior-erosion plus last chambers from the convex hull.

    Early chambers whose foramina are narrower than 2*radius are recovered by
    :func:`exterior_erosion_interior`. Chambers that stay open to the exterior
    (typically the last chamber, with a wide primary aperture) are recovered by
    watershed of the void inside the convex hull, keeping only basins whose
    enclosure ratio (shell contact / void contact) is at least ``min_ratio``.
    The umbilicus and other open bays have a low ratio and are discarded.
    """
    shell = np.asarray(shell_mask).astype(bool)
    hull = _convex_hull_mask(shell)
    cand = hull & ~shell
    kept = np.zeros(shell.shape, bool)
    dropped = np.zeros(shell.shape, bool)
    if np.any(cand):
        lab, _ = segment_chambers(cand, 1.0, smooth_sigma=1.0, h_ratio=h_ratio,
                                  min_chamber_voxels=min_chamber_voxels)
        for i in range(1, int(lab.max()) + 1):
            C = lab == i
            if int(C.sum()) < min_chamber_voxels:
                dropped |= C
                continue
            if _enclosure_ratio(C, shell) >= min_ratio:
                kept |= C
            else:
                dropped |= C
    # Early, fully enclosed chambers (narrow foramina) that the hull watershed
    # may have merged or missed, excluding the discarded umbilicus / open bays.
    early = exterior_erosion_interior(shell, radius)["interior"]
    interior = kept | (early & ~dropped)
    filled = interior | shell
    exterior = ~shell & ~interior
    return {"shell": shell, "closed": shell, "filled": filled, "interior": interior,
            "method": "hull_enclosure", "exterior": exterior,
            "early_interior": early}


def sealed_test(shell_mask: np.ndarray, close_radius: int = 8, *, method: str = "hull_enclosure", **kw) -> dict:
    """Seal apertures and return closed, filled and interior masks.

    ``method='hull_enclosure'`` (default): exterior erosion plus last chambers
    recovered inside the convex hull (see :func:`hull_enclosure_interior`).
    ``method='exterior_erosion'``: cut narrow apertures from the outside only.
    ``method='closing'``: morphological closing of the shell plus hole filling
    (Speijer et al. 2008 / Briguglio-style sealing). A large closing can wrap
    exterior bays into the envelope; prefer exterior_erosion for chambers.
    """
    shell = np.asarray(shell_mask).astype(bool)
    if method == "hull_enclosure":
        return hull_enclosure_interior(shell, close_radius, **{k: kw[k] for k in
            ("h_ratio", "min_chamber_voxels", "min_ratio") if k in kw})
    if method == "exterior_erosion":
        return exterior_erosion_interior(shell, close_radius)
    if method != "closing":
        raise ValueError("method: 'hull_enclosure', 'exterior_erosion' or 'closing'")
    closed = _close_mask(shell, close_radius)
    filled = ndi.binary_fill_holes(closed)
    # Subtract the original shell, not the closed mask: a large closing can
    # temporarily fill small lumens, and those voxels must stay in the interior.
    interior = filled & ~shell
    return {"shell": shell, "closed": closed, "filled": filled, "interior": interior,
            "method": "closing"}


def filled_test_volume(shell_mask: np.ndarray, voxel_size: float, *, close_radius: int = 0) -> float:
    """Volume (um^3) of the test after filling enclosed chambers (binary fill holes).

    If ``close_radius`` is positive, a spherical morphological closing of that
    radius (voxels) is applied first so that foram apertures smaller than the
    closing diameter can be sealed. Without closing, chambers that open to the
    exterior through the aperture are not filled.
    """
    filled = sealed_test(shell_mask, close_radius, method="closing")["filled"]
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
    if abs(a) < 1e-18:
        raise ValueError("grey values of the standards are degenerate")
    return lambda grey: (np.asarray(grey, float) - b) / a


def estimate_air_grey(grey: np.ndarray, shell_mask: np.ndarray, *, dilate: int = 3) -> float:
    """Median grey of background voxels (dilated away from the shell), a proxy for air."""
    m = np.asarray(shell_mask).astype(bool)
    if dilate:
        m = ndi.binary_dilation(m, iterations=dilate)
    bg = np.asarray(grey)[~m]
    if bg.size == 0:
        raise ValueError("no background voxels to estimate air")
    return float(np.median(bg))


def relative_ct(grey, air_grey: float, calcite_grey: float) -> np.ndarray:
    """Iwasaki et al. 2015 relative CT number: air = 0, calcite standard = 1000."""
    denom = float(calcite_grey) - float(air_grey)
    if abs(denom) < 1e-12:
        raise ValueError("air and calcite grey values must differ")
    return 1000.0 * (np.asarray(grey, float) - float(air_grey)) / denom


def density_from_relative_ct(ct, *, calcite_density: float = CALCITE_DENSITY_G_CM3) -> np.ndarray:
    """Linear map of Iwasaki CT to density, with the calcite standard at ``calcite_density``."""
    return np.asarray(ct, float) / 1000.0 * float(calcite_density)


def density_index(grey: np.ndarray, shell_mask: np.ndarray, calibration=None, *, erode: int = 1) -> float:
    """Mean grey value (or calibrated density) of the shell interior, eroded by ``erode``
    voxels to limit partial-volume bias at the surfaces."""
    m = shell_mask.astype(bool)
    if erode:
        m = ndi.binary_erosion(m, iterations=erode)
    vals = np.asarray(grey)[m]
    if vals.size == 0:
        raise ValueError("shell too thin for the erosion")
    v = float(vals.mean())
    return float(calibration(v)) if calibration else v


def segment_test(grey: np.ndarray, *, method: str = "otsu", min_size: int = 64,
                 largest_only: bool = True) -> tuple[np.ndarray, float, dict]:
    """Segment the calcite shell from a grey microCT volume.

    Baseline: Otsu threshold with the shell treated as the bright phase, optional
    removal of small objects, and optionally keeping only the largest connected
    component (the specimen). Returns the uint8 shell mask (1 = calcite), the
    threshold, and a small info dict.
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


def separate_interior(shell_mask: np.ndarray, *, method: str = "hull_enclosure", close_radius: int = 8,
                      n_rays: int = 64):
    """Split chamber lumen from the exterior.

    ``method='hull_enclosure'`` (default): exterior erosion plus hull-recovered last chambers.
    ``method='closing'``: morphological closing plus hole filling.
    ``method='ambient_occlusion'`` is reserved for the MFSE/Titschack reimplementation.
    Returns (interior_mask uint8, info).
    """
    if method == "ambient_occlusion":
        raise NotImplementedError("Phase 3: ambient-occlusion interior separation")
    if method not in ("closing", "exterior_erosion", "hull_enclosure"):
        raise ValueError("method: 'closing', 'exterior_erosion', 'hull_enclosure' or 'ambient_occlusion'")
    parts = sealed_test(shell_mask, close_radius, method=method)
    info = {"method": method, "close_radius": close_radius,
            "n_interior_voxels": int(parts["interior"].sum()),
            "n_filled_voxels": int(parts["filled"].sum())}
    return parts["interior"].astype(np.uint8), info


def _growth_path(centroids: np.ndarray, volumes: np.ndarray) -> list[int]:
    """Order chambers by walking from the smallest (proloculus candidate) to its
    nearest unused neighbour. A first approximation of growth sequence for a
    trochospiral test (Speijer et al. 2008); a kummerform last chamber is then
    last spatially rather than by volume."""
    n = len(volumes)
    if n == 0:
        return []
    start = int(np.argmin(volumes))
    order = [start]
    left = set(range(n)) - {start}
    while left:
        last = order[-1]
        nxt = min(left, key=lambda i: float(np.linalg.norm(centroids[i] - centroids[last])))
        order.append(nxt)
        left.remove(nxt)
    return order


def segment_chambers(interior_mask: np.ndarray, voxel_size: float = 1.0, seeds=None,
                     *, smooth_sigma: float = 1.0, h_ratio: float = 0.3,
                     min_chamber_voxels: int = 50, grey: np.ndarray | None = None):
    """Label individual chambers in an enclosed lumen.

    Reimplementation of the ForamJ description (Trend et al. 2026): optional
    Gaussian smoothing of the interior, EDT, seeds at extended EDT maxima,
    3D watershed. Labels are renumbered 1..n along a nearest-neighbour path
    that starts at the smallest chamber (a proloculus candidate). Returns
    (labels int32, info dict).
    """
    from skimage.morphology import h_maxima
    from skimage.segmentation import watershed

    m = np.asarray(interior_mask).astype(bool)
    if not np.any(m):
        return np.zeros(np.asarray(interior_mask).shape, np.int32), {
            "n_chambers": 0, "volumes_um3": [], "centroids_um": [],
            "equivalent_radius_um": [], "mean_grey": [],
            "order": "growth_path_from_smallest", "notes": ["empty interior"],
        }
    work = m
    if smooth_sigma and smooth_sigma > 0:
        sm = ndi.gaussian_filter(m.astype(np.float32), smooth_sigma)
        work = sm > 0.5
        work &= m  # seeds stay inside the lumen; watershed uses the original mask
    edt = ndi.distance_transform_edt(work)
    if seeds is None:
        h = max(float(h_ratio) * float(edt.max()), 1.0)
        peaks = h_maxima(edt, h=h)
        peaks &= work
        seed_labels, n_seed = ndi.label(peaks)
        if n_seed == 0:
            # fallback: single seed at the EDT peak
            seed_labels = np.zeros(edt.shape, np.int32)
            peak = tuple(int(i) for i in np.unravel_index(int(np.argmax(edt)), edt.shape))
            seed_labels[peak] = 1
            n_seed = 1
        # put each seed at the EDT maximum of its plateau
        seeds = np.zeros(edt.shape, np.int32)
        for i in range(1, n_seed + 1):
            pts = seed_labels == i
            if not np.any(pts):
                continue
            loc = np.argwhere(pts)
            vals = edt[pts]
            p = tuple(int(x) for x in loc[int(np.argmax(vals))])
            seeds[p] = i
    else:
        seeds = np.asarray(seeds)
        n_seed = int(seeds.max()) if seeds.size else 0
    # Watershed the original lumen (not the smoothed mask) so chamber volumes
    # are not shrunk by the Gaussian. Seeds come from the smoothed EDT.
    edt_lumen = ndi.distance_transform_edt(m)
    labels = watershed(-edt_lumen, seeds, mask=m)
    # drop tiny fragments
    lab_ids = [i for i in range(1, int(labels.max()) + 1) if int((labels == i).sum()) >= min_chamber_voxels]
    cleaned = np.zeros(labels.shape, np.int32)
    for new_i, old in enumerate(lab_ids, 1):
        cleaned[labels == old] = new_i
    labels = cleaned
    leftover = m & (labels == 0)
    if np.any(leftover) and int(labels.max()) > 0:
        _, inds = ndi.distance_transform_edt(labels == 0, return_indices=True)
        labels[leftover] = labels[tuple(inds)][leftover]
    n = int(labels.max())
    vs = float(voxel_size)
    volumes, centroids, radii, mean_grey = [], [], [], []
    for i in range(1, n + 1):
        pts = labels == i
        nv = int(pts.sum())
        vol = nv * vs**3
        volumes.append(vol)
        com = ndi.center_of_mass(pts)
        centroids.append(tuple(float(c) * vs for c in com))
        radii.append(float((3.0 * vol / (4.0 * np.pi)) ** (1.0 / 3.0)))
        if grey is not None:
            mean_grey.append(float(np.asarray(grey)[pts].mean()) if nv else float("nan"))
    if n:
        order = _growth_path(np.asarray(centroids), np.asarray(volumes))
        remap = {old + 1: new for new, old in enumerate(order, 1)}
        relabeled = np.zeros_like(labels)
        for old, new in remap.items():
            relabeled[labels == old] = new
        labels = relabeled
        volumes = [volumes[i] for i in order]
        centroids = [centroids[i] for i in order]
        radii = [radii[i] for i in order]
        if mean_grey:
            mean_grey = [mean_grey[i] for i in order]
    info = {
        "n_chambers": n,
        "volumes_um3": volumes,
        "centroids_um": centroids,
        "equivalent_radius_um": radii,
        "mean_grey": mean_grey,
        "order": "growth_path_from_smallest",
        "h_ratio": h_ratio,
        "smooth_sigma": smooth_sigma,
        "min_chamber_voxels": min_chamber_voxels,
        "notes": [
            "Watershed of the enclosed lumen (ForamJ method, reimplemented). "
            "Labels 1..n follow a nearest-neighbour path starting at the smallest chamber.",
        ],
    }
    return labels.astype(np.int32), info


def wall_thickness(shell_mask: np.ndarray, voxel_size: float, *, split_septa: bool = True,
                   method: str = "dt"):
    """Local thickness of the shell as a per-voxel diameter map (um).

    Hildebrand and Ruegsegger (1997) define local thickness as the diameter of the
    largest inscribed sphere. PoreSpy's ``local_thickness`` returns a radius in
    voxels; this function multiplies by 2 and by the voxel size. ``split_septa``
    is reserved for a later outer-wall vs septa split and is currently ignored.
    Returns (thickness_um ndarray, stats dict).
    """
    import porespy as ps

    m = np.asarray(shell_mask).astype(bool)
    out = np.zeros(m.shape, np.float32)
    stats = {"method": "porespy.local_thickness (Hildebrand diameter)", "n": 0,
             "mean": float("nan"), "median": float("nan"), "std": float("nan"),
             "min": float("nan"), "max": float("nan"), "unit": "um, diameter",
             "split_septa": bool(split_septa)}
    if not np.any(m):
        return out, stats
    # radius in voxels -> diameter in um
    radius_vox = ps.filters.local_thickness(m, method=method)
    thick = (2.0 * np.asarray(radius_vox, float) * float(voxel_size)).astype(np.float32)
    out[m] = thick[m]
    vals = out[m]
    vals = vals[vals > 0]
    if vals.size:
        stats.update({
            "n": int(vals.size),
            "mean": float(vals.mean()),
            "median": float(np.median(vals)),
            "std": float(vals.std()),
            "min": float(vals.min()),
            "max": float(vals.max()),
            "p10": float(np.percentile(vals, 10)),
            "p90": float(np.percentile(vals, 90)),
        })
    return out, stats


def chamber_wall_thickness(thickness_um: np.ndarray, shell_mask: np.ndarray,
                           chamber_labels: np.ndarray, *, dilate: int = 2) -> list[dict]:
    """Thickness statistics on shell voxels adjacent to each chamber label."""
    shell = np.asarray(shell_mask).astype(bool)
    labels = np.asarray(chamber_labels)
    rows = []
    for i in range(1, int(labels.max()) + 1 if labels.size else 1):
        adj = ndi.binary_dilation(labels == i, iterations=dilate) & shell
        vals = thickness_um[adj]
        vals = vals[vals > 0]
        rows.append({
            "chamber": i,
            "n": int(vals.size),
            "mean": float(vals.mean()) if vals.size else float("nan"),
            "median": float(np.median(vals)) if vals.size else float("nan"),
            "std": float(vals.std()) if vals.size else float("nan"),
        })
    return rows


def _shell_sample(grey, shell_mask, erode: int):
    m = np.asarray(shell_mask).astype(bool)
    if erode:
        eroded = ndi.binary_erosion(m, iterations=erode)
        if np.any(eroded):
            m = eroded
    return m, np.asarray(grey)[m]


def calibrate_shell_density(grey, shell_mask, calibration: dict | None, *, erode: int = 1) -> dict:
    """Calibrate shell greys to density, or flag the result as uncalibrated.

    ``calibration`` keys (all optional except that a calcite reference is required
    for a calibrated result):

    * ``method``: ``iwasaki_ct`` (default when air and calcite greys are given) or
      ``linear`` (two or more standards).
    * ``air_grey``, ``calcite_grey``, ``calcite_density_g_cm3`` (default 2.71).
    * ``grey_standards``, ``density_standards`` for a general linear fit.
    * ``estimate_air`` (default True): median background grey if ``air_grey`` is omitted.

    Without a calcite (or second) standard the mean shell grey is returned as
    ``density_index`` and ``density_calibrated`` is False.
    """
    notes = []
    cal = dict(calibration or {})
    mask, vals = _shell_sample(grey, shell_mask, erode)
    if vals.size == 0:
        raise ValueError("shell too thin for the erosion")
    mean_grey = float(vals.mean())
    out = {
        "density_index": mean_grey,
        "density_calibrated": False,
        "density_method": "uncalibrated_mean_grey",
        "shell_density_g_cm3": None,
        "relative_ct_mean": None,
        "air_grey": None,
        "calcite_grey": None,
        "notes": notes,
        "sample_mask": mask,
    }
    air = cal.get("air_grey")
    if air is None and cal.get("estimate_air", True):
        try:
            air = estimate_air_grey(grey, shell_mask)
            notes.append(f"air_grey estimated from background median ({air:.4g}).")
        except ValueError as e:
            notes.append(f"air_grey not estimated: {e}")
    calcite = cal.get("calcite_grey")
    rho_std = float(cal.get("calcite_density_g_cm3", CALCITE_DENSITY_G_CM3))
    method = cal.get("method")
    std_g = cal.get("grey_standards")
    std_d = cal.get("density_standards")
    if std_g is not None and std_d is not None:
        method = method or "linear"
        fn = linear_density_calibration(std_g, std_d)
        dens = fn(vals)
        out.update({
            "density_calibrated": True,
            "density_method": "linear_standards",
            "shell_density_g_cm3": float(np.mean(dens)),
            "density_index": float(np.mean(dens)),
        })
        notes.append("Density from a linear fit on supplied grey and density standards.")
    elif air is not None and calcite is not None:
        method = method or "iwasaki_ct"
        ct = relative_ct(vals, air, calcite)
        dens = density_from_relative_ct(ct, calcite_density=rho_std)
        out.update({
            "density_calibrated": True,
            "density_method": "iwasaki_ct",
            "relative_ct_mean": float(ct.mean()),
            "shell_density_g_cm3": float(np.mean(dens)),
            "density_index": float(np.mean(dens)),
            "air_grey": float(air),
            "calcite_grey": float(calcite),
        })
        notes.append(
            f"Iwasaki relative CT (air={float(air):.4g}, calcite standard={float(calcite):.4g} "
            f"at {rho_std} g/cm3). CT_mean={float(ct.mean()):.4g}."
        )
    else:
        out["air_grey"] = float(air) if air is not None else None
        notes.append(
            "UNCALIBRATED: no calcite (or second) standard was given. density_index is the "
            "mean shell grey value, not g/cm3. Shell mass and test bulk density are omitted."
        )
    out["notes"] = notes
    return out


def _mass_ug(density_g_cm3: float, volume_um3: float) -> float:
    # 1 um^3 = 1e-12 cm^3; 1 g = 1e6 ug
    return float(density_g_cm3) * float(volume_um3) * 1e-6


def qc_chambers(grey, labels, path, *, title="Chambers"):
    """Mid-slice (the slice that cuts the most labels) of grey and coloured chambers."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lab = np.asarray(labels)
    g = np.asarray(grey)
    # slice that contains the most labelled voxels (more representative than unique-count)
    counts = [int((lab[z] > 0).sum()) for z in range(lab.shape[0])]
    z = int(np.argmax(counts)) if counts else lab.shape[0] // 2
    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].imshow(g[z], cmap="gray")
    ax[0].set_title("grey")
    cmap = plt.get_cmap("nipy_spectral").copy().with_extremes(under="black")
    vmax = max(int(lab.max()), 1) + 0.5
    ax[1].imshow(g[z], cmap="gray")
    im = ax[1].imshow(np.ma.masked_where(lab[z] == 0, lab[z]), cmap=cmap, vmin=0.5,
                      vmax=vmax, alpha=0.7)
    ax[1].set_title(f"chambers (n={int(lab.max())}, z={z})")
    ticks = list(range(1, int(lab.max()) + 1)) if int(lab.max()) <= 20 else None
    fig.colorbar(im, ax=ax[1], shrink=0.8, ticks=ticks)
    for a in ax:
        a.axis("off")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def qc_thickness(thickness_um, shell_mask, path, *, title="Wall thickness"):
    """Mid-slice thickness map and a histogram of shell diameters (um)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    t = np.asarray(thickness_um)
    m = np.asarray(shell_mask).astype(bool)
    counts = [int((m[z]).sum()) for z in range(m.shape[0])]
    z = int(np.argmax(counts)) if counts else t.shape[0] // 2
    vals = t[m]
    vals = vals[vals > 0]
    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    sl = np.ma.masked_where(~m[z], t[z])
    im = ax[0].imshow(sl, cmap="viridis")
    ax[0].set_title(f"diameter (um), z={z}")
    ax[0].axis("off")
    fig.colorbar(im, ax=ax[0], shrink=0.8)
    if vals.size:
        ax[1].hist(vals, bins=min(40, max(8, int(np.sqrt(vals.size)))), color="0.3")
        ax[1].axvline(float(np.median(vals)), color="C1", label=f"median {np.median(vals):.2f}")
        ax[1].axvline(float(vals.mean()), color="C0", linestyle="--", label=f"mean {vals.mean():.2f}")
        ax[1].legend(fontsize=8)
    ax[1].set_xlabel("local thickness diameter (um)")
    ax[1].set_ylabel("voxels")
    fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def analyse(grey: np.ndarray, voxel_size: float, **kw) -> ForamResult:
    """End-to-end analysis: segment the shell, seal apertures, measure volume and
    surface, label chambers, map wall thickness, and apply density calibration
    when standards are given.

    Keyword arguments include those of :func:`segment_test`, plus ``close_radius``
    (default 8), ``density_erode`` (default 1), ``calibration`` (dict or None),
    ``chambers`` (default True), ``h_ratio``, ``smooth_sigma``,
    ``min_chamber_voxels``, ``wall_thickness_map`` (default True), and ``output``
    (optional directory for QC figures and maps).
    """
    erode = kw.pop("density_erode", 1)
    close_radius = int(kw.pop("close_radius", 8))
    seal_method = kw.pop("seal_method", "hull_enclosure")
    calibration = kw.pop("calibration", None)
    do_chambers = bool(kw.pop("chambers", True))
    do_thick = bool(kw.pop("wall_thickness_map", True))
    h_ratio = float(kw.pop("h_ratio", 0.2))
    smooth_sigma = float(kw.pop("smooth_sigma", 1.0))
    min_ch = int(kw.pop("min_chamber_voxels", 50))
    output = kw.pop("output", None)
    name = kw.pop("name", "foram")
    shell, thr, info = segment_test(grey, **kw)
    vs = float(voxel_size)
    sv = shell_volume(shell, vs)
    parts = sealed_test(shell, close_radius, method=seal_method,
                        h_ratio=h_ratio, min_chamber_voxels=min_ch)
    fv = float(np.count_nonzero(parts["filled"])) * vs**3
    hv = convex_hull_volume(shell, vs)
    sa = shell_surface_area(shell, vs, method="marching_cubes")
    test_phi = float((fv - sv) / fv) if fv > 0 else float("nan")
    calc = float(sv / fv) if fv > 0 else float("nan")
    hull_frac = float(sv / hv) if hv > 0 else float("nan")
    notes = [
        f"Test envelope and chamber lumen use seal_method={seal_method} with radius {close_radius} voxels.",
        "test_porosity is enclosed void volume over filled test volume; calcification_fraction is shell over filled test volume; shell_to_hull_fraction is shell over the convex hull.",
    ]
    dens = calibrate_shell_density(np.asarray(grey), shell, calibration, erode=erode)
    notes.extend(dens["notes"])
    mass = None
    bulk = None
    if dens["density_calibrated"] and dens["shell_density_g_cm3"] is not None:
        mass = _mass_ug(dens["shell_density_g_cm3"], sv)
        bulk = float(dens["shell_density_g_cm3"] * calc) if calc == calc else None  # NaN-safe
        notes.append(
            f"Shell mass {mass:.4g} ug from calibrated mean density "
            f"{dens['shell_density_g_cm3']:.4g} g/cm3 times shell volume."
        )

    labels = np.zeros(shell.shape, np.int32)
    ch_info = {"n_chambers": 0, "volumes_um3": [], "centroids_um": [],
               "equivalent_radius_um": [], "mean_grey": [], "order": None}
    if do_chambers:
        labels, ch_info = segment_chambers(
            parts["interior"], vs, smooth_sigma=smooth_sigma, h_ratio=h_ratio,
            min_chamber_voxels=min_ch, grey=np.asarray(grey),
        )
        notes.append(
            f"Chamber count {ch_info['n_chambers']}; labels follow a growth-path order "
            "starting at the smallest chamber (proloculus candidate)."
        )

    thick_map = np.zeros(shell.shape, np.float32)
    thick_stats: dict = {}
    ch_thick: list[dict] = []
    if do_thick:
        thick_map, thick_stats = wall_thickness(shell, vs)
        if ch_info["n_chambers"]:
            ch_thick = chamber_wall_thickness(thick_map, shell, labels)
        notes.append(
            "Wall thickness is the Hildebrand diameter (um) from PoreSpy local thickness."
        )

    if output:
        out = Path(output)
        out.mkdir(parents=True, exist_ok=True)
        if do_chambers:
            qc_chambers(grey, labels, out / "qc_chambers.png", title=f"{name}: chambers")
            np.save(out / "chamber_labels.npy", labels)
        if do_thick:
            qc_thickness(thick_map, shell, out / "qc_thickness.png", title=f"{name}: wall thickness")
            np.save(out / "wall_thickness_um.npy", thick_map)

    return ForamResult(
        shell_volume_um3=sv,
        filled_volume_um3=fv,
        surface_area_um2=sa,
        test_porosity=test_phi,
        calcification_fraction=calc,
        shell_to_hull_fraction=hull_frac,
        density_index=dens["density_index"],
        density_calibrated=bool(dens["density_calibrated"]),
        density_method=dens["density_method"],
        shell_density_g_cm3=dens["shell_density_g_cm3"],
        shell_mass_ug=mass,
        test_bulk_density_g_cm3=bulk,
        relative_ct_mean=dens["relative_ct_mean"],
        air_grey=dens["air_grey"],
        calcite_grey=dens["calcite_grey"],
        threshold=thr,
        n_shell_voxels=info.get("n_shell_voxels"),
        n_chambers=ch_info.get("n_chambers"),
        chamber_volumes_um3=list(ch_info.get("volumes_um3") or []),
        chamber_equivalent_radius_um=list(ch_info.get("equivalent_radius_um") or []),
        chamber_centroids_um=list(ch_info.get("centroids_um") or []),
        chamber_mean_grey=list(ch_info.get("mean_grey") or []),
        chamber_order=ch_info.get("order"),
        wall_thickness_um=thick_stats,
        chamber_wall_thickness_um=ch_thick,
        convex_hull_volume_um3=hv,
        close_radius_voxels=close_radius,
        notes=notes,
    )


def run_config(config: dict | str | Path) -> dict:
    """Run :func:`analyse` from a YAML config or dict (used by ``morphoct foram``)."""
    import json

    import yaml

    from morphoct.io import load_volume
    from morphoct.report import write_json

    if not isinstance(config, dict):
        config = yaml.safe_load(Path(config).read_text())
    cfg = dict(config)
    inp = dict(cfg["input"])
    path = inp.pop("path")
    vol = load_volume(path, **inp)
    vs = vol.isotropic_voxel
    out = Path(cfg.get("output", f"morphoct_out/{cfg.get('name', 'foram')}"))
    out.mkdir(parents=True, exist_ok=True)
    fcfg = dict(cfg.get("foram") or {})
    res = analyse(vol.data, vs, output=out, name=cfg.get("name", "foram"), **fcfg)
    payload = {"name": cfg.get("name", "foram"), "config": cfg,
               "voxel_size_um": vs, "shape_zyx": list(vol.shape),
               "foram": res.as_dict()}
    write_json(payload, out / "results.json")
    (out / "notes.txt").write_text("\n".join(res.notes) + "\n")
    return payload
