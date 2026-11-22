"""Command line: ``morphoct run config.yaml`` and ``morphoct foram``."""

from __future__ import annotations

import argparse
import json
import sys


def _foram_cfg_from_args(a) -> dict:
    import yaml

    cfg: dict = {}
    if a.config:
        cfg = yaml.safe_load(open(a.config)) or {}
    inp = dict(cfg.get("input") or {})
    if a.input:
        inp["path"] = a.input
    if a.shape:
        inp["shape"] = [int(x) for x in a.shape.replace("x", ",").split(",")]
    if a.dtype:
        inp["dtype"] = a.dtype
    if a.voxel_size is not None:
        inp["voxel_size"] = a.voxel_size
    if a.byteorder:
        inp["byteorder"] = a.byteorder
    if a.offset is not None:
        inp["offset"] = a.offset
    if not inp.get("path"):
        raise SystemExit("morphoct foram: pass a YAML config or --input")
    cfg["input"] = inp
    if a.output:
        cfg["output"] = a.output
    if a.name:
        cfg["name"] = a.name
    fcfg = dict(cfg.get("foram") or {})
    if a.close_radius is not None:
        fcfg["close_radius"] = a.close_radius
    if a.no_chambers:
        fcfg["chambers"] = False
    if a.no_thickness:
        fcfg["wall_thickness_map"] = False
    if a.h_ratio is not None:
        fcfg["h_ratio"] = a.h_ratio
    if a.smooth_sigma is not None:
        fcfg["smooth_sigma"] = a.smooth_sigma
    if a.min_chamber_voxels is not None:
        fcfg["min_chamber_voxels"] = a.min_chamber_voxels
    cal = dict(fcfg.get("calibration") or {})
    if a.calibration:
        extra = yaml.safe_load(open(a.calibration)) if a.calibration.endswith((".yml", ".yaml")) else json.loads(a.calibration)
        cal.update(extra or {})
    if a.air_grey is not None:
        cal["air_grey"] = a.air_grey
    if a.calcite_grey is not None:
        cal["calcite_grey"] = a.calcite_grey
    if a.calcite_density is not None:
        cal["calcite_density_g_cm3"] = a.calcite_density
    if a.density_method:
        cal["method"] = a.density_method
    if cal:
        fcfg["calibration"] = cal
    cfg["foram"] = fcfg
    if "output" not in cfg:
        cfg["output"] = f"morphoct_out/{cfg.get('name', 'foram')}"
    if "name" not in cfg:
        cfg["name"] = "foram"
    return cfg


def main(argv=None) -> int:
    p = argparse.ArgumentParser(prog="morphoct", description="microCT analysis of rocks, cores and foraminifera")
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run the pipeline described by a YAML config")
    r.add_argument("config")
    sub.add_parser("version", help="print versions")
    f = sub.add_parser("foram", help="segment a foraminifera test: chambers, wall thickness, density")
    f.add_argument("config", nargs="?", help="YAML config (input.path, foram options, output)")
    f.add_argument("--input", help="volume path (RAW, TIFF, NPY, ...)")
    f.add_argument("--shape", help="Z,Y,X for headerless RAW (x fastest on disk)")
    f.add_argument("--dtype", default=None, help="RAW dtype (default: from config or uint16)")
    f.add_argument("--voxel-size", type=float, dest="voxel_size", help="isotropic voxel size in um")
    f.add_argument("--byteorder", default=None, help="little or big")
    f.add_argument("--offset", type=int, default=None, help="RAW header offset in bytes")
    f.add_argument("--output", help="output directory")
    f.add_argument("--name", help="specimen name, used in QC titles")
    f.add_argument("--close-radius", type=int, dest="close_radius", help="aperture-closing radius in voxels")
    f.add_argument("--no-chambers", action="store_true", help="skip chamber labelling")
    f.add_argument("--no-thickness", action="store_true", help="skip the wall-thickness map")
    f.add_argument("--h-ratio", type=float, dest="h_ratio", help="EDT h-maxima ratio for chamber seeds")
    f.add_argument("--smooth-sigma", type=float, dest="smooth_sigma")
    f.add_argument("--min-chamber-voxels", type=int, dest="min_chamber_voxels")
    f.add_argument("--calibration", help="YAML or JSON file of density-calibration keys")
    f.add_argument("--air-grey", type=float, dest="air_grey")
    f.add_argument("--calcite-grey", type=float, dest="calcite_grey",
                   help="grey value of a calcite standard; required for a calibrated density")
    f.add_argument("--calcite-density", type=float, dest="calcite_density",
                   help="density of the calcite standard in g/cm3 (default 2.71)")
    f.add_argument("--density-method", dest="density_method", choices=("iwasaki_ct", "linear"))
    a = p.parse_args(argv)
    if a.cmd == "version":
        from morphoct.pipeline import _versions

        print(json.dumps(_versions()))
        return 0
    if a.cmd == "foram":
        from morphoct.foram import run_config

        if a.dtype is None and not a.config:
            a.dtype = "uint16"
        if a.byteorder is None and not a.config:
            a.byteorder = "little"
        cfg = _foram_cfg_from_args(a)
        payload = run_config(cfg)
        res = payload["foram"]
        cal = "calibrated" if res.get("density_calibrated") else "UNCALIBRATED"
        print(
            f"done: {payload['name']}: chambers {res.get('n_chambers')}; "
            f"shell {res['shell_volume_um3']:.4g} um3; wall median "
            f"{(res.get('wall_thickness_um') or {}).get('median', float('nan')):.3g} um; "
            f"density {cal}; report in {cfg['output']}"
        )
        return 0
    from morphoct.pipeline import run

    res = run(a.config)
    out = res["config"].get("output", "")
    print(f"done in {res['total_seconds']} s: porosity {res['morphometrics']['porosity']:.4f}; report in {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
