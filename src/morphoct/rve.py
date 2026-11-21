"""Representative volume element (REV/RVE) of a single scan by the window method.

Kanit et al. (2003): the variance of a property over sub-volumes of volume V scales as
D^2(V) = K V^(-gamma) (gamma = 1 with a finite integral range); the volume for which n
realisations give the mean within a relative error eps (95 %) is
V_REV = (4 K / (eps^2 Z^2 n))^(1/gamma). The window statistics and the fit are
morphoGen's (``morphogen.rve``); morphoCT adds the confidence interval of V_REV by Monte
Carlo over the fitted (ln K, gamma) and its expression as a cube edge.
"""

from __future__ import annotations

import numpy as np

from morphoct.io import to_voxel_image


def default_windows(shape, n: int = 7, smallest: int = 8) -> list[int]:
    largest = min(shape) // 2
    w = np.unique(np.geomspace(smallest, max(largest, smallest + 1), n).astype(int))
    return [int(x) for x in w if x >= smallest]


def rev_porosity(binary: np.ndarray, voxel_size: float, *, windows=None, eps_rel=0.05,
                 n_real: int = 1, fit_from_edge_voxels: int | None = None, seed: int = 0) -> dict:
    from morphogen import rve as mrve

    b = (np.asarray(binary) == 1).astype(np.uint8)
    img = to_voxel_image(b, voxel_size, periodic=False)
    windows = windows or default_windows(b.shape)
    stats = mrve.window_statistics(img, windows, phase=1)
    vmin = None
    if fit_from_edge_voxels:
        vmin = (fit_from_edge_voxels * voxel_size) ** 3
    else:
        # leave out windows smaller than ~3 correlation lengths: assume the smallest third
        vmin = float(stats.volume[min(2, len(stats.volume) - 3)]) if len(stats.volume) >= 5 else None
    vr = (vmin, None) if vmin else None
    law = mrve.fit_statistics(stats, volume_range=vr)
    mean = float(np.mean(b))
    v_rev = float(law.rve_volume(eps_rel, mean, n_real))
    rng = np.random.default_rng(seed)
    if np.isfinite(law.gamma_se) and np.isfinite(law.log_k_se):
        g = rng.normal(law.gamma, law.gamma_se, 4000)
        lk = rng.normal(np.log(law.K), law.log_k_se, 4000)
        g = g[g > 0.2]
        lk = lk[: g.size]
        vs = (4.0 * np.exp(lk) / (eps_rel**2 * mean**2 * n_real)) ** (1.0 / g)
        lo, hi = np.percentile(vs, [2.5, 97.5])
    else:
        lo = hi = float("nan")
    edge = lambda v: float(np.cbrt(v))  # noqa: E731
    return {
        "property": "porosity",
        "eps_rel": eps_rel,
        "n_realisations": n_real,
        "windows_voxels": windows,
        "table": stats.rows(),
        "fit": {"K": float(law.K), "gamma": float(law.gamma), "gamma_ci95": list(map(float, law.gamma_interval())),
                "r2": float(law.r2), "n_points": int(law.n_points)},
        "rev_volume_um3": v_rev,
        "rev_edge_um": edge(v_rev),
        "rev_edge_um_ci95": [edge(lo), edge(hi)] if np.isfinite(lo) else [float("nan")] * 2,
        "rev_edge_voxels": edge(v_rev) / voxel_size,
        "sample_edge_um": float(np.cbrt(np.prod(b.shape)) * voxel_size),
        "sample_is_rev": bool(np.prod(b.shape) * voxel_size**3 >= v_rev),
    }


def rev_figure(rev: dict, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = rev["table"]
    v = np.array([r["volume"] for r in t])
    var = np.array([r["variance"] for r in t])
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    ax.loglog(v, var, "o", label="windows")
    vv = np.geomspace(v.min(), max(v.max(), rev["rev_volume_um3"]) * 1.2, 50)
    ax.loglog(vv, rev["fit"]["K"] * vv ** (-rev["fit"]["gamma"]), "-", label=f"fit, gamma={rev['fit']['gamma']:.2f}")
    ax.set_xlabel("window volume (um^3)")
    ax.set_ylabel("variance of porosity")
    ax.legend()
    ax.set_title(f"REV edge {rev['rev_edge_um']:.0f} um (eps {rev['eps_rel']:.0%})")
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path
