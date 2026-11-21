"""The configured pipeline behind ``morphoct run config.yaml``.

Config keys (all but ``input`` optional)::

    name: fontainebleau
    output: out/fontainebleau
    input: {path: sample.raw.gz, format: raw, shape: [300, 288, 288], dtype: uint8,
            voxel_size: 7.5, crop: [[0, 300], [0, 288], [0, 288]]}
    binary: {pore_value: 0}          # input already segmented: which value is pore
    preprocess: {denoise: median, size: 3}
    segment: {method: otsu, pore: dark, min_size: 0}
    morphometrics: {pore_sizes: true}
    rve: {eps_rel: 0.05, windows: [8, 16, 24, 32, 48, 64, 96]}
    properties: {subvolume: 144, axes: [0], mirror: true, permeability: true,
                 formation_factor: true, max_steps: 30000}
    network: {axes: [0, 1, 2], geometry: spheres_and_cylinders}
"""

from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import yaml

from morphoct import io as mio
from morphoct import morphometrics, network, preprocess, properties, report, rve, segment


def _versions() -> dict:
    import morphogen
    import porespy

    import morphoct

    return {"morphoct": morphoct.__version__, "morphogen": morphogen.__version__, "porespy": porespy.__version__}


def run(config: dict | str | Path) -> dict:
    if not isinstance(config, dict):
        config = yaml.safe_load(Path(config).read_text())
    cfg = config
    name = cfg.get("name", "sample")
    out = Path(cfg.get("output", f"morphoct_out/{name}"))
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    timings = {}
    results: dict = {"name": name, "config": cfg, "versions": _versions(), "notes": []}
    figures = []

    inp = dict(cfg["input"])
    path = inp.pop("path")
    vol = mio.load_volume(path, **inp)
    vs = vol.isotropic_voxel
    timings["load"] = time.time() - t0

    gray = None
    if "binary" in cfg:
        pv = cfg["binary"].get("pore_value", 1)
        binary = (vol.data == pv).astype(np.uint8)
        results["segmentation"] = {"method": "provided", "pore_value": pv}
    else:
        pp = cfg.get("preprocess", {})
        if pp.get("crop"):
            vol = preprocess.crop(vol, **pp["crop"])
        dn = pp.get("denoise", "median")
        kw = {k: v for k, v in pp.items() if k not in ("denoise", "crop")}
        vol = preprocess.denoise(vol, dn, **kw)
        gray = vol.data
        sg = cfg.get("segment", {})
        method = sg.get("method", "otsu")
        if method == "otsu":
            binary, thr = segment.otsu(vol, pore=sg.get("pore", "dark"))
            thrs = [thr]
        elif method == "multi_otsu":
            labels, thrs = segment.multi_otsu(vol, sg.get("classes", 3))
            binary = (labels == 0).astype(np.uint8)
            thr = thrs[0]
        else:
            raise ValueError(f"segment.method {method!r}")
        binary = segment.clean(binary, sg.get("min_size", 0))
        results["segmentation"] = {"method": method, "threshold": thr, "thresholds": thrs,
                                   "sensitivity": segment.threshold_sensitivity(vol, thr, pore=sg.get("pore", "dark"))}
        figures.append(segment.qc_figure(gray, binary, out / "qc_segmentation.png", thresholds=thrs))
    timings["segment"] = time.time() - t0 - sum(timings.values())
    if "binary" in cfg:
        figures.append(segment.qc_figure(None, binary, out / "qc_segmentation.png", title="Provided segmentation"))
    mio.save_volume(out / "binary.npz", vol.with_data(binary))

    mcfg = cfg.get("morphometrics", {})
    m = morphometrics.compute(binary, vs, pore_sizes=mcfg.get("pore_sizes", True))
    results["morphometrics"] = m
    figures.append(report.s2_figure(m, out / "two_point_correlation.png"))
    timings["morphometrics"] = time.time() - t0 - sum(timings.values())

    if cfg.get("rve", {}) is not False:
        rc = cfg.get("rve", {}) or {}
        rv = rve.rev_porosity(binary, vs, windows=rc.get("windows"), eps_rel=rc.get("eps_rel", 0.05))
        results["rve"] = rv
        figures.append(rve.rev_figure(rv, out / "rev_porosity.png"))
    timings["rve"] = time.time() - t0 - sum(timings.values())

    pc = cfg.get("properties", {})
    if pc:
        sub = binary
        if pc.get("subvolume"):
            s = int(pc["subvolume"])
            size = [min(s, n) for n in binary.shape]
            st = [(n - z) // 2 for n, z in zip(binary.shape, size)]
            sub = binary[st[0]:st[0] + size[0], st[1]:st[1] + size[1], st[2]:st[2] + size[2]]
            results["properties_subvolume"] = {"start_zyx": st, "size_zyx": size, "porosity": float(sub.mean())}
            results["notes"].append(f"Voxel properties computed on a centred {size} subvolume (porosity {sub.mean():.4f}) to bound runtime.")
        axes = tuple(pc.get("axes", [0]))
        mirror = pc.get("mirror", True)
        if pc.get("permeability", True):
            results["permeability"] = properties.permeability(sub, vs, axes=axes, mirror_flow_axis=mirror,
                                                              max_steps=pc.get("max_steps"), tol=pc.get("tol", 1e-4))
            timings["permeability"] = time.time() - t0 - sum(timings.values())
        if pc.get("formation_factor", True):
            results["formation_factor"] = properties.formation_factor(sub, vs, axes=axes, mirror_flow_axis=mirror)
            timings["formation_factor"] = time.time() - t0 - sum(timings.values())
        if mirror:
            results["notes"].append("Periodic solvers: sample mirrored along the flow axis; transverse faces wrapped.")
        mk = m["minkowski"]["S_V_per_um"]
        results["kozeny_carman_mD"] = properties.kozeny_carman(m["porosity"], mk)

    nc = cfg.get("network")
    if nc:
        results["network"] = network.permeability(binary, vs, axes=tuple(nc.get("axes", [0, 1, 2])),
                                                  geometry=nc.get("geometry", "spheres_and_cylinders"),
                                                  pore_diameter=nc.get("pore_diameter", "inscribed"),
                                                  throat_diameter=nc.get("throat_diameter", "inscribed"))
        timings["network"] = time.time() - t0 - sum(timings.values())

    results["timings_s"] = {k: round(v, 1) for k, v in timings.items()}
    results["total_seconds"] = round(time.time() - t0, 1)
    report.write_json(results, out / "results.json")
    report.write_html(results, figures, out / "report.html")
    report.write_pdf(results, figures, out / "report.pdf")
    return results
