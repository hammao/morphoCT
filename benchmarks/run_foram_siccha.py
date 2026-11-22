"""Chamber count, wall thickness and density on Siccha et al. 2023 specimens."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
DATA = Path("/workspace/morphoCT_data/foram_siccha2023")
OUT = ROOT / "benchmarks" / "results" / "foram_siccha"

SPECIMENS = [
    {
        "name": "CA_GGL_1",
        "species": "Globigerinita glutinata",
        "path": DATA / "CA_GGL_1_207x205x187-1.2016_um-ushort_LE.raw",
        "shape": [187, 205, 207],  # Z,Y,X
        "dtype": "uint16",
        "byteorder": "little",
        "voxel_size": 1.2016,
        # Final-whorl counts from taxonomic diagnoses, not total chambers.
        "literature": {
            "final_whorl_chambers": "3-3.5",
            "source": "Pearson et al. 2018 / Mikrotax diagnosis of G. glutinata: 3 to 3.5 chambers in the final whorl. Total chamber number is not published for this specimen.",
        },
    },
    {
        "name": "GE_GTE_1",
        "species": "Globigerinoides tenellus",
        "path": DATA / "GE_GTE_1_227x220x241-1.2004_um-ushort_LE.raw",
        "shape": [241, 220, 227],
        "dtype": "uint16",
        "byteorder": "little",
        "voxel_size": 1.2004,
        "literature": {
            "final_whorl_chambers": "4",
            "source": "Mikrotax diagnosis of G. tenellus: 4 rapidly enlarging chambers in the final whorl. Total chamber number is not published for this specimen.",
        },
    },
    {
        "name": "CA_GUV_1",
        "species": "Globigerinita uvula",
        "path": DATA / "CA_GUV_1_152x180x138-1.2016_um-ushort_LE.raw",
        "shape": [138, 180, 152],
        "dtype": "uint16",
        "byteorder": "little",
        "voxel_size": 1.2016,
        "literature": {
            "final_whorl_chambers": "3-3.5",
            "source": "Mikrotax diagnosis of G. uvula: 3 to 3.5 chambers in the final whorl; high trochospire, so total chambers typically exceed the final-whorl count. Total chamber number is not published for this specimen.",
        },
    },
]


def main() -> int:
    sys.path.insert(0, str(ROOT / "src"))
    from morphoct.foram import run_config

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for spec in SPECIMENS:
        name = spec["name"]
        out = OUT / name
        print(f"=== {name} ({spec['species']}) ===", flush=True)
        t0 = time.time()
        cfg = {
            "name": name,
            "output": str(out),
            "input": {
                "path": str(spec["path"]),
                "shape": spec["shape"],
                "dtype": spec["dtype"],
                "byteorder": spec["byteorder"],
                "voxel_size": spec["voxel_size"],
            },
            "foram": {
                "close_radius": 8,
                "seal_method": "hull_enclosure",
                "min_size": 64,
                "largest_only": True,
                "chambers": True,
                "h_ratio": 0.2,
                "smooth_sigma": 1.0,
                "min_chamber_voxels": 80,
                "wall_thickness_map": True,
            },
        }
        (out).mkdir(parents=True, exist_ok=True)
        (out / "config.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))
        payload = run_config(cfg)
        elapsed = time.time() - t0
        foram = payload["foram"]
        row = {
            "name": name,
            "species": spec["species"],
            "elapsed_s": elapsed,
            "literature": spec["literature"],
            "n_chambers": foram.get("n_chambers"),
            "chamber_volumes_um3": foram.get("chamber_volumes_um3"),
            "shell_volume_um3": foram.get("shell_volume_um3"),
            "filled_volume_um3": foram.get("filled_volume_um3"),
            "test_porosity": foram.get("test_porosity"),
            "shell_to_hull_fraction": foram.get("shell_to_hull_fraction"),
            "surface_area_um2": foram.get("surface_area_um2"),
            "density_calibrated": foram.get("density_calibrated"),
            "density_method": foram.get("density_method"),
            "density_index": foram.get("density_index"),
            "air_grey": foram.get("air_grey"),
            "shell_density_g_cm3": foram.get("shell_density_g_cm3"),
            "shell_mass_ug": foram.get("shell_mass_ug"),
            "wall_thickness_um": foram.get("wall_thickness_um"),
            "chamber_wall_thickness_um": foram.get("chamber_wall_thickness_um"),
            "threshold": foram.get("threshold"),
            "notes": foram.get("notes"),
        }
        rows.append(row)
        print(
            f"  chambers={row['n_chambers']}  wall median="
            f"{(row['wall_thickness_um'] or {}).get('median')} um  "
            f"density={row['density_method']}  {elapsed:.1f}s",
            flush=True,
        )
    summary = {
        "dataset": "Siccha et al. 2023, PANGAEA 949585, CC-BY-4.0",
        "doi": "10.1594/PANGAEA.949585",
        "specimens": rows,
    }
    (OUT / "summary.json").write_text(json.dumps(summary, indent=2))
    print(f"wrote {OUT / 'summary.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
