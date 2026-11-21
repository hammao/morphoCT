"""Phase 1 benchmark against Andrä et al. (2013) digital rock physics benchmarks.

Data: https://github.com/fkrzikalla/drp-benchmarks (supplement of Andrä et al. 2013,
Part I doi:10.1016/j.cageo.2012.09.005 and Part II doi:10.1016/j.cageo.2012.09.008).
Published values: results/simulations.ods of that repository (sheet "2.2 Overview" and
"Details").

1. Fontainebleau (ExxonMobil segmentation, 288 x 288 x 300, 7.5 um): full morphoCT
   pipeline, compared with the published porosity, permeability and formation factor.
2. Berea (VSG segmentation, the benchmark 400^3 subvolume at 0.74 um): porosity,
   morphometrics, REV and network permeability.
3. Berea grey-scale slices 126 to 205 (same x, y window): morphoCT median + Otsu
   segmentation compared voxel by voxel with the published VSG segmentation.

Usage: python benchmarks/run_benchmarks.py [--data DIR] [--out DIR] [--quick]
"""

from __future__ import annotations

import argparse
import json
import time
import urllib.request
from pathlib import Path

import numpy as np

BASE = "https://github.com/fkrzikalla/drp-benchmarks/raw/master/"
PUBLISHED = {
    "fontainebleau": {
        "porosity": 0.147,
        "k_mD": {"EJ Stokes mean": 1913.7, "LBM mean": 1502.6, "Stokes z (Fraunhofer)": 2073.0, "LBM z (Kongju)": 1624.1},
        "F": {"EJ diffusion": 41.6, "FEM-C": 29.4, "FEM-F": 29.6},
        "source": "Andrä et al. 2013 Part II, simulations.ods sheets '2.2 Overview' and 'Details'",
    },
    "berea_sub400": {
        "porosity": 0.184,
        "k_mD": {"EJ Stokes x": 123.6, "LBM x": 113.4, "EJ Stokes mean": 105.6, "LBM mean": 95.5},
        "F": {"EJ diffusion": 22.5, "FEM-C": 22.0, "FEM-F": 20.2, "FDM": 29.7},
        "source": "Andrä et al. 2013 Part II, simulations.ods sheet 'Details' (VSG subvolume x,y 151-550, z 126-525)",
    },
}


def fetch(rel: str, data: Path) -> Path:
    p = data / Path(rel).name
    if not p.exists():
        data.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(BASE + rel, p)
    return p


def fontainebleau(data: Path, out: Path, quick: bool) -> dict:
    from morphoct.pipeline import run

    raw = fetch("images/fontainebleau/segmented-exxon.raw.gz", data)
    cfg = {
        "name": "fontainebleau_exxon", "output": str(out / "fontainebleau"),
        "input": {"path": str(raw), "format": "raw", "shape": [300, 288, 288], "dtype": "uint8", "voxel_size": 7.5},
        "binary": {"pore_value": 0},
        "rve": {"eps_rel": 0.05, "windows": [8, 12, 16, 24, 32, 48, 64, 96, 144]},
        "properties": {"subvolume": 96 if quick else 144, "axes": [0], "mirror": True, "tol": 1e-4},
        "network": {"axes": [0, 1, 2]},
    }
    return run(cfg)


def berea_seg(data: Path, out: Path) -> dict:
    """Segmentation check on 80 grey slices against the VSG segmentation."""
    import gzip

    import tifffile

    from morphoct import preprocess, segment
    from morphoct.io import Volume

    vsg = fetch("images/berea/segmented-vsg.raw.gz", data)
    lab = np.frombuffer(gzip.open(vsg).read(), np.uint8).reshape(1024, 1024, 1024)
    z0, z1 = 125, 205  # slices 126..205 (1-based file numbering)
    ref = (lab[z0:z1, 150:550, 150:550] == 0).astype(np.uint8)
    sub400 = np.ascontiguousarray(lab[125:525, 150:550, 150:550])
    del lab
    sl = []
    for i in range(z0 + 1, z1 + 1):
        f = fetch(f"images/berea/tif/B6_10_{i:04d}.rec.16bit.tif", data / "berea_tif")
        sl.append(tifffile.imread(f)[150:550, 150:550])
    g = Volume(np.stack(sl).astype(np.float32), 0.74)
    res = {}
    for name, kw in (("median3", {"method": "median", "size": 3}), ("nlm", {"method": "nlm", "h": 0.8})):
        t = time.time()
        d = preprocess.denoise(g, kw.pop("method"), **kw)
        m, thr = segment.otsu(d)
        res[name] = {"porosity": float(m.mean()), "dice_vs_vsg": segment.dice(m, ref), "threshold": thr,
                     "sensitivity": segment.threshold_sensitivity(d, thr), "seconds": round(time.time() - t, 1)}
        if name == "median3":
            segment.qc_figure(d.data, m, out / "berea_grey_qc.png", thresholds=[thr], title="Berea grey: median + Otsu")
    res["vsg_porosity_same_slices"] = float(ref.mean())
    return res, sub400


def berea_sub400(binary: np.ndarray, out: Path) -> dict:
    from morphoct import morphometrics, network, rve

    t = time.time()
    m = morphometrics.compute(binary == 0, 0.74, pore_sizes=False)
    r = rve.rev_porosity((binary == 0).astype(np.uint8), 0.74, windows=[10, 16, 25, 40, 64, 100, 133, 200])
    n = network.permeability((binary == 0).astype(np.uint8), 0.74, axes=(0, 1, 2))
    return {"morphometrics": m, "rve": r, "network": n, "seconds": round(time.time() - t, 1)}


def network_sensitivity(binary_pore: np.ndarray, voxel: float) -> list[dict]:
    """Network permeability (z) for every conduit geometry and diameter choice."""
    from morphoct import network

    snow = network.extract_dict(binary_pore, voxel)
    rows = []
    for pdm, tdm in (("inscribed", "inscribed"), ("equivalent", "inscribed"), ("equivalent", "equivalent")):
        for geo in network.GEOMETRIES:
            try:
                r = network.permeability(binary_pore, voxel, axes=(0,), snow_network=snow,
                                       geometry=geo, pore_diameter=pdm, throat_diameter=tdm)
                rows.append({"geometry": geo, "pore_d": pdm, "throat_d": tdm, "k_mD": r["axes"]["z"]["k_mD"]})
            except Exception as e:  # incompatible geometry/diameters
                rows.append({"geometry": geo, "pore_d": pdm, "throat_d": tdm, "error": str(e)[:80]})
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data")
    ap.add_argument("--out", default="benchmarks/out")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-berea", action="store_true")
    ap.add_argument("--sensitivity", action="store_true")
    a = ap.parse_args()
    data, out = Path(a.data), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    summary = {"published": PUBLISHED}
    f = fontainebleau(data, out, a.quick)
    summary["fontainebleau"] = {k: f[k] for k in ("morphometrics", "rve", "permeability", "formation_factor", "network",
                                                  "kozeny_carman_mD", "properties_subvolume", "timings_s", "total_seconds") if k in f}
    if a.sensitivity:
        b = np.load(out / "fontainebleau" / "binary.npz")["data"].astype(np.uint8)
        summary["fontainebleau_network_sensitivity"] = network_sensitivity(b, 7.5)
    if not a.skip_berea:
        seg, sub = berea_seg(data, out)
        summary["berea_grey_segmentation"] = seg
        summary["berea_sub400"] = berea_sub400(sub, out)
    summary["total_seconds"] = round(time.time() - t0, 1)
    from morphoct.report import write_json

    write_json(summary, out / "benchmark_summary.json")
    print(json.dumps({"total_seconds": summary["total_seconds"]}))


if __name__ == "__main__":
    main()
