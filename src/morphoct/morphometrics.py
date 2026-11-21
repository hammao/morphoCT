"""Morphometrics of a segmented pore space (1 = pore, 0 = solid), in physical units.

Delegates to morphoGen (Minkowski densities, covariance, chords, connectivity) and
PoreSpy (radial two-point correlation, pore size distribution). Conventions follow
Torquato (2002) for correlation functions and Ohser and Muecklich (2000) / Vogel et al.
(2010) for Minkowski functionals.
"""

from __future__ import annotations

import numpy as np

from morphoct.io import to_voxel_image


def porosity(binary: np.ndarray) -> float:
    return float(np.mean(np.asarray(binary) == 1))


def compute(binary: np.ndarray, voxel_size: float, *, max_lag: int | None = None,
            pore_sizes: bool = True, chords: bool = True) -> dict:
    """All MVP morphometrics. Lengths in micrometres, densities per micrometre (S_V),
    per micrometre squared (M_V) and per cubic micrometre (chi_V)."""
    import morphogen as mg

    b = (np.asarray(binary) == 1).astype(np.uint8)
    img = to_voxel_image(b, voxel_size, periodic=False)
    out: dict = {"voxel_size_um": float(voxel_size), "shape_zyx": list(b.shape)}
    phi = porosity(b)
    out["porosity"] = phi

    mk = mg.analysis.minkowski(img, phase=1, periodic=False)
    out["minkowski"] = {
        "V_V": float(mk.volume_fraction),
        "S_V_per_um": float(mk.specific_surface),
        "M_V_per_um2": float(mk.mean_curvature),
        "chi_V_per_um3": float(mk.euler_density),
    }
    out["specific_surface_per_um"] = float(mk.specific_surface)
    # per unit solid volume, the petrophysical convention (S/V_solid)
    out["specific_surface_solid_per_um"] = float(mk.specific_surface / max(1 - phi, 1e-12))

    comp = mg.analysis.connected_components(img, phase=1, connectivity=6, periodic=False)
    out["connectivity"] = {
        "n_pore_clusters": int(comp.n_components),
        "spanning_fraction": float(comp.spanning_fraction),
        "connected_porosity": float(phi * comp.spanning_fraction),
        # morphoGen axes are x, y, z; report in z, y, x order of the volume
        "percolating_zyx": [bool(comp.percolating.get(k, False)) for k in ("z", "y", "x")],
    }

    lag = max_lag or max(4, min(b.shape) // 4)
    cov = mg.analysis.covariance(img, phase=1, max_lag=lag, periodic=False)
    lags, s2 = cov.mean()
    out["two_point_correlation"] = {"r_um": [float(x) for x in lags], "S2": [float(x) for x in s2]}
    # correlation length: first lag where (S2 - phi^2) falls below 1/e of (phi - phi^2)
    norm = (np.asarray(s2) - phi**2) / max(phi - phi**2, 1e-12)
    below = np.nonzero(norm < np.exp(-1))[0]
    out["correlation_length_um"] = float(lags[below[0]]) if below.size else float("nan")

    if chords:
        cl = {}
        for ax_mg, name in ((2, "z"), (1, "y"), (0, "x")):
            c = mg.analysis.chord_lengths(img, phase=1, axis=ax_mg, periodic=False)
            cl[name] = float(np.mean(c)) if c.size else float("nan")
        out["mean_chord_length_um"] = cl
        out["mean_chord_length_pooled_um"] = float(np.nanmean(list(cl.values())))
        # stereological check: mean pore chord = 4 phi / S_V
        out["chord_stereology_um"] = float(4 * phi / max(mk.specific_surface, 1e-12))

    if pore_sizes:
        import porespy as ps

        lt = ps.filters.local_thickness(b.astype(bool), sizes=25)
        vals = lt[b == 1] * voxel_size
        if vals.size:
            out["pore_size_um"] = {
                "d10": float(np.percentile(vals, 10)),
                "d50": float(np.percentile(vals, 50)),
                "d90": float(np.percentile(vals, 90)),
                "note": "local thickness (maximal inscribed sphere diameter), volume weighted",
            }
    return out
