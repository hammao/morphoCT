"""Foraminifera module (planned API; Phase 3).

Status: API stub. Only generic, self-written helpers are implemented
(:func:`shell_volume`, :func:`linear_density_calibration`, :func:`density_index`).
The chamber and interior methods are specified from the published papers and will be
written from scratch; no code is copied from ForamJ (no licence file in its repository)
or MFSE (GPL-3.0, MATLAB; compatible, but we port the method, not the code).

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

from dataclasses import dataclass, field

import numpy as np


@dataclass
class ForamResult:
    shell_volume_um3: float
    filled_volume_um3: float | None = None
    chamber_volumes_um3: list[float] = field(default_factory=list)
    chamber_centroids_um: list[tuple[float, float, float]] = field(default_factory=list)
    wall_thickness_um: dict = field(default_factory=dict)
    density_index: float | None = None
    notes: list[str] = field(default_factory=list)


def shell_volume(shell_mask: np.ndarray, voxel_size: float) -> float:
    """Calcite volume (um^3) of a binary shell mask."""
    return float(np.count_nonzero(shell_mask)) * voxel_size**3


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
    from scipy import ndimage as ndi

    m = shell_mask.astype(bool)
    if erode:
        m = ndi.binary_erosion(m, iterations=erode)
    vals = grey[m]
    if vals.size == 0:
        raise ValueError("shell too thin for the erosion")
    v = float(vals.mean())
    return float(calibration(v)) if calibration else v


def segment_test(grey: np.ndarray, *, method: str = "otsu"):
    """Planned: shell vs void vs infill segmentation (ForamJ steps 1 to 6)."""
    raise NotImplementedError("Phase 3: shell, infill and background segmentation")


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
    """Planned end-to-end foram analysis: segment_test -> separate_interior ->
    segment_chambers -> wall_thickness -> density_index."""
    raise NotImplementedError("Phase 3")
