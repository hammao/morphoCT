"""Transport properties on voxels with morphoGen's solvers.

* Absolute permeability: morphoGen's two-relaxation-time lattice Boltzmann Stokes solver
  (D3Q19, bounce-back walls; Ginzburg 2005), ``morphogen.homogenization.permeability_tensor``.
* Formation factor: FFT solution of the electrical conduction problem with a conducting
  brine (sigma = 1) and an insulating solid (sigma = 0): F = sigma_w / sigma_eff
  (Archie 1942); tortuosity factor tau = F * phi.

Both solvers are periodic. A cropped rock sample is not, so by default the sample is
mirrored along the flow axis (:func:`morphoct.preprocess.mirror`) before solving. The
transverse faces remain wrapped; this is the usual digital-rock compromise and is
reported in the results.
"""

from __future__ import annotations

import time

import numpy as np

from morphoct.io import to_voxel_image
from morphoct.preprocess import mirror

UM2_PER_MD = 9.869233e-4  # 1 mD = 9.869233e-16 m^2 = 9.869233e-4 um^2
_MG_AXIS = {0: 2, 1: 1, 2: 0}  # volume axis (z, y, x) -> morphoGen axis (x, y, z)
_NAMES = {0: "z", 1: "y", 2: "x"}


def permeability(binary: np.ndarray, voxel_size: float, *, axes=(0,), mirror_flow_axis=True,
                 tol: float = 1e-4, max_steps: int | None = None, relaxation: float = 4.0) -> dict:
    """Permeability (mD) along each axis in ``axes`` (volume order z=0, y=1, x=2)."""
    from morphogen.homogenization import permeability_tensor

    b = (np.asarray(binary) == 1).astype(np.uint8)
    res = {"method": "morphoGen TRT lattice Boltzmann D3Q19", "mirrored": bool(mirror_flow_axis), "axes": {}}
    for a in axes:
        t0 = time.time()
        arr = mirror(b, (a,)) if mirror_flow_axis else b
        img = to_voxel_image(arr, voxel_size, periodic=True)
        m = _MG_AXIS[a]
        out = permeability_tensor(img, pores=1, axes=(m,), tol=tol, max_steps=max_steps, relaxation=relaxation)
        k_um2 = float(out["tensor"][m, m])
        res["axes"][_NAMES[a]] = {
            "k_um2": k_um2,
            "k_mD": k_um2 / UM2_PER_MD,
            "steps": int(out["steps"]),
            "error_estimate": float(out["error"]),
            "seconds": round(time.time() - t0, 1),
        }
    vals = [v["k_mD"] for v in res["axes"].values()]
    res["mean_mD"] = float(np.mean(vals))
    return res


def formation_factor(binary: np.ndarray, voxel_size: float, *, axes=(0,), mirror_flow_axis=True,
                     tol: float = 1e-6) -> dict:
    """Formation factor along each axis by the FFT conduction solver (pore sigma = 1,
    solid sigma = 0)."""
    from morphogen.homogenization import Material, homogenize

    b = (np.asarray(binary) == 1).astype(np.uint8)
    phi = float(b.mean())
    mats = {0: Material(electrical_conductivity=0.0), 1: Material(electrical_conductivity=1.0)}
    res = {"method": "morphoGen FFT, electrical conduction, insulating solid", "mirrored": bool(mirror_flow_axis), "axes": {}}
    for a in axes:
        t0 = time.time()
        arr = mirror(b, (a,)) if mirror_flow_axis else b
        img = to_voxel_image(arr, voxel_size, periodic=True)
        out = homogenize(img, mats, properties=["electrical_conductivity_tensor"], tol=tol)
        tens = np.asarray(out["electrical_conductivity_tensor.periodic"]["tensor"])
        m = _MG_AXIS[a]
        sig = float(tens[m, m])
        F = 1.0 / sig if sig > 0 else float("inf")
        axis = {"sigma_eff": sig, "F": F, "tortuosity_factor": F * phi,
                "seconds": round(time.time() - t0, 1)}
        if sig <= 0:
            axis["note"] = "non-positive effective conductivity; the subvolume may not percolate along this axis"
        res["axes"][_NAMES[a]] = axis
    res["mean_F"] = float(np.mean([v["F"] for v in res["axes"].values()]))
    # Archie cementation exponent from F = phi^-m
    res["archie_m"] = float(np.log(res["mean_F"]) / -np.log(phi)) if 0 < phi < 1 else float("nan")
    return res


def kozeny_carman(phi: float, specific_surface_per_um: float, c: float = 5.0) -> float:
    """Kozeny-Carman estimate k = phi^3 / (c S_V^2) in mD (S_V per unit bulk volume),
    a sanity check on the solvers."""
    return phi**3 / (c * specific_surface_per_um**2) / UM2_PER_MD
