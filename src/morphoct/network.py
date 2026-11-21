"""Pore network extraction (PoreSpy SNOW2; Gostick 2017) and single-phase flow in
OpenPNM (Gostick et al. 2016): network permeability and statistics, a fast cross-check of
the voxel lattice Boltzmann result.
"""

from __future__ import annotations

import time

import numpy as np

UM2_PER_MD = 9.869233e-4
_FACES = {0: ("xmin", "xmax"), 1: ("ymin", "ymax"), 2: ("zmin", "zmax")}  # PoreSpy axis i = array axis i
_NAMES = {0: "z", 1: "y", 2: "x"}


def extract_dict(binary: np.ndarray, voxel_size: float) -> dict:
    """SNOW2 network of the pore phase as PoreSpy's dict (SI units, metres)."""
    import porespy as ps

    return ps.networks.snow2((np.asarray(binary) == 1), voxel_size=voxel_size * 1e-6).network


def extract(binary: np.ndarray, voxel_size: float, snow_network: dict | None = None):
    """SNOW2 network of the pore phase as an OpenPNM Network (SI units, metres)."""
    import openpnm as op

    return op.io.network_from_porespy(snow_network if snow_network is not None else extract_dict(binary, voxel_size))


GEOMETRIES = ("spheres_and_cylinders", "cones_and_cylinders", "pyramids_and_cuboids", "cubes_and_cuboids")


def _prepare(pn, geometry="spheres_and_cylinders", pore_diameter="inscribed", throat_diameter="inscribed"):
    import openpnm as op

    pn["pore.diameter"] = pn[f"pore.{pore_diameter}_diameter"]
    # conduit models need throats narrower than both adjacent pores
    dmin = pn["pore.diameter"][pn["throat.conns"]].min(axis=1)
    pn["throat.diameter"] = np.minimum(pn[f"throat.{throat_diameter}_diameter"], 0.99 * dmin)
    mods = getattr(op.models.collections.geometry, geometry)
    for name in ("pore.diameter", "throat.diameter", "pore.seed", "throat.seed", "pore.max_size", "throat.max_size"):
        mods.pop(name, None)
    pn.add_model_collection(mods)
    pn.regenerate_models()
    return pn


def permeability(binary: np.ndarray, voxel_size: float, *, axes=(0, 1, 2), snow_network: dict | None = None,
                 geometry: str = "spheres_and_cylinders", pore_diameter: str = "inscribed",
                 throat_diameter: str = "inscribed") -> dict:
    """Network permeability (mD) along each axis with water (Hagen-Poiseuille conduits).

    The conduit geometry and the diameters taken from SNOW2 (``inscribed``,
    ``equivalent``, ``extended``) change the result several-fold (on the Fontainebleau
    benchmark from 340 to 4300 mD): treat the network value as an uncalibrated estimate
    and use the lattice Boltzmann value as the reference."""
    import openpnm as op

    t0 = time.time()
    pn = extract(binary, voxel_size, snow_network)
    t_ext = time.time() - t0
    n_raw = int(pn.Np)
    h = op.utils.check_network_health(pn)
    trim = list(h["disconnected_pores"]) + list(h["isolated_pores"])
    if trim:
        op.topotools.trim(network=pn, pores=np.unique(trim))
    _prepare(pn, geometry, pore_diameter, throat_diameter)
    water = op.phase.Water(network=pn)
    water.add_model_collection(op.models.collections.physics.basic)
    water.regenerate_models()
    mu = float(np.mean(water["pore.viscosity"]))
    shape = np.asarray(binary).shape
    L_side = np.array(shape) * voxel_size * 1e-6
    res = {"method": f"PoreSpy snow2 + OpenPNM StokesFlow ({geometry}, {pore_diameter} pores, {throat_diameter} throats)", "n_pores_raw": n_raw,
           "n_pores": int(pn.Np), "n_throats": int(pn.Nt),
           "mean_coordination": float(2 * pn.Nt / max(pn.Np, 1)), "extract_seconds": round(t_ext, 1), "axes": {}}
    for a in axes:
        lo, hi = _FACES[a]
        inlet, outlet = pn.pores(lo), pn.pores(hi)
        if inlet.size == 0 or outlet.size == 0:
            res["axes"][_NAMES[a]] = {"k_mD": 0.0, "note": "no spanning pores"}
            continue
        sf = op.algorithms.StokesFlow(network=pn, phase=water)
        sf.set_value_BC(pores=inlet, values=1.0)
        sf.set_value_BC(pores=outlet, values=0.0)
        sf.run()
        Q = float(np.abs(sf.rate(pores=inlet))[0])
        L = L_side[a]
        A = np.prod(np.delete(L_side, a))
        k = Q * mu * L / (A * 1.0)
        res["axes"][_NAMES[a]] = {"k_m2": k, "k_mD": k * 1e12 / UM2_PER_MD}
    res["mean_mD"] = float(np.mean([v["k_mD"] for v in res["axes"].values()]))
    res["seconds"] = round(time.time() - t0, 1)
    return res
