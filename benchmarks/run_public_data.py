"""Public-data benchmark: Imperial College 2015 carbonates and Siccha et al. 2023 forams.

Carbonate volumes are the 1000^3 segmented images from the Imperial College
Pore-Scale Modelling 2015 release (Estaillades, Ketton), loaded from
/workspace/morphoCT_data (not redistributed here). Published transport values
for Estaillades are taken from the dataset .mhd header (DNS on the full image).

Foraminifera volumes are individual raw files from Siccha et al. 2023,
PANGAEA doi:10.1594/PANGAEA.949585 (CC-BY-4.0).

Usage:
    python benchmarks/run_public_data.py [--data DIR] [--out DIR] [--quick]
"""

from __future__ import annotations

import argparse
import gzip
import json
import time
from pathlib import Path

import numpy as np

# 1 mD = 9.869233e-16 m^2
M2_PER_MD = 9.869233e-16

PUBLISHED = {
    "estaillades_2015": {
        "source": (
            "Imperial College Pore-Scale Modelling group, 2015 release, "
            "Estaillades_1000c_3p31136um.mhd (DNS block Estaillades_1000c_3p31136um-1-X). "
            "Total image porosity is measured here (the header reports effective porosity only)."
        ),
        "shape": [1000, 1000, 1000],
        "voxel_size_um": 3.31136,
        "pore_value": 0,
        "effective_porosity": 0.109066,
        "K_x_m2": 2.18531e-13,
        "K_x_mD": 2.18531e-13 / M2_PER_MD,
        "FF_x": 166.697,
    },
    "ketton_2015": {
        "source": (
            "Imperial College Pore-Scale Modelling group, 2015 release, "
            "Ketton_1000c_3p00006um.mhd. The header has no bundled DNS; only the image "
            "geometry is published with the file."
        ),
        "shape": [1000, 1000, 1000],
        "voxel_size_um": 3.00006,
        "pore_value": 0,
    },
}


def stream_porosity(path: Path, pore_value=0, dtype=np.uint8, chunk=1000 * 1000 * 20) -> dict:
    """Histogram and porosity of a headerless raw.gz volume, without loading it all."""
    hist = np.zeros(256, np.int64)
    n = 0
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as f:
        while True:
            buf = f.read(chunk)
            if not buf:
                break
            a = np.frombuffer(buf, dtype)
            if a.size == 0:
                break
            hist += np.bincount(a.astype(np.uint8, copy=False) if a.dtype != np.uint8 else a, minlength=256)
            n += a.size
    unique = [int(i) for i in np.nonzero(hist)[0]]
    return {
        "n_voxels": int(n),
        "unique_values": unique,
        "porosity": float(hist[int(pore_value)] / n) if n else float("nan"),
        "value_counts": {str(i): int(hist[i]) for i in unique},
    }


def read_raw_gz_crop(path: Path, full_shape, crop_zyx, dtype=np.uint8) -> np.ndarray:
    """Read a z,y,x crop from a headerless gzip raw volume (x fastest). full_shape is z,y,x."""
    z, y, x = (int(s) for s in full_shape)
    (z0, z1), (y0, y1), (x0, x1) = crop_zyx
    out = np.empty((z1 - z0, y1 - y0, x1 - x0), dtype=dtype)
    item = np.dtype(dtype).itemsize
    row = x * item
    sl = y * row
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rb") as f:
        if z0:
            remain = z0 * sl
            while remain:
                got = f.read(remain)
                if not got:
                    raise EOFError(path)
                remain -= len(got)
        buf_size = sl
        for iz in range(z1 - z0):
            raw = b""
            while len(raw) < buf_size:
                chunk = f.read(buf_size - len(raw))
                if not chunk:
                    raise EOFError(path)
                raw += chunk
            slc = np.frombuffer(raw, dtype).reshape(y, x)
            out[iz] = slc[y0:y1, x0:x1]
    return np.ascontiguousarray(out)


def centre_crop(shape, edge: int):
    e = [min(edge, n) for n in shape]
    st = [(n - z) // 2 for n, z in zip(shape, e)]
    return tuple((s, s + z) for s, z in zip(st, e))


def percolates_x(pore) -> bool:
    from scipy import ndimage as ndi

    labels, n = ndi.label(pore.astype(bool))
    if n == 0:
        return False
    left = set(np.unique(labels[:, :, 0])) - {0}
    right = set(np.unique(labels[:, :, -1])) - {0}
    return bool(left & right)


def find_porosity_matched_window(pore, edge: int, target_phi: float, *, step: int = 32):
    """Return the ``edge``^3 window whose porosity is closest to ``target_phi`` and
    that percolates along x, plus a summary of all windows on the grid."""
    z, y, x = pore.shape
    edge = min(edge, z, y, x)
    rows = []
    for iz in range(0, z - edge + 1, step):
        for iy in range(0, y - edge + 1, step):
            for ix in range(0, x - edge + 1, step):
                sub = pore[iz:iz + edge, iy:iy + edge, ix:ix + edge]
                phi = float(sub.mean())
                rows.append((abs(phi - target_phi), phi, (iz, iy, ix)))
    rows.sort()
    picked = None
    checked = []
    for _d, phi, st in rows[:40]:
        sub = pore[st[0]:st[0] + edge, st[1]:st[1] + edge, st[2]:st[2] + edge]
        ok = percolates_x(sub)
        checked.append({"start_zyx": list(st), "porosity": phi, "percolates_x": ok})
        if ok and picked is None:
            picked = (st, phi)
    phis = np.array([r[1] for r in rows], float)
    stats = {
        "edge": edge, "step": step, "n_windows": len(rows),
        "porosity_min": float(phis.min()), "porosity_max": float(phis.max()),
        "porosity_mean": float(phis.mean()), "porosity_std": float(phis.std()),
        "target_porosity": float(target_phi),
        "nearest_checked": checked[:15],
    }
    return picked, stats


def run_carbonate(name: str, raw: Path, published: dict, out: Path, *, crop_edge: int,
                  prop_edge: int, quick: bool) -> dict:
    from morphoct import morphometrics, properties, rve, segment
    from morphoct.report import write_json

    t0 = time.time()
    vs = published["voxel_size_um"]
    shape = published["shape"]
    pv = published.get("pore_value", 0)
    full = stream_porosity(raw, pore_value=pv)
    crop = centre_crop(shape, 200 if quick else crop_edge)
    arr = read_raw_gz_crop(raw, shape, crop)
    pore = (arr == pv).astype(np.uint8)
    del arr
    mcfg_sizes = not quick
    morph = morphometrics.compute(pore, vs, pore_sizes=mcfg_sizes)
    windows = [8, 12, 16, 24, 32, 48, 64] if quick else [8, 12, 16, 24, 32, 48, 64, 96, 128, 160]
    windows = [w for w in windows if w <= min(pore.shape) // 2]
    rev = rve.rev_porosity(pore, vs, windows=windows, eps_rel=0.05)
    # Properties: a centred window can sit in a tight, low-porosity region of a
    # heterogeneous carbonate. Search a larger centre crop for a window whose
    # porosity matches the full image and that percolates along x (MHD DNS axis).
    pe = min(64 if quick else prop_edge, min(shape))
    search_edge = 256 if quick else 512
    search_crop = centre_crop(shape, search_edge)
    search = read_raw_gz_crop(raw, shape, search_crop)
    search_pore = (search == pv).astype(np.uint8)
    del search
    target = float(full["porosity"])
    picked, search_stats = find_porosity_matched_window(search_pore, pe, target, step=32 if not quick else 64)
    if picked is None:
        st = [(n - pe) // 2 for n in search_pore.shape]
        sub = search_pore[st[0]:st[0] + pe, st[1]:st[1] + pe, st[2]:st[2] + pe]
        selection = "fallback: centred window (no percolating match)"
    else:
        st, _phi = picked
        sub = search_pore[st[0]:st[0] + pe, st[1]:st[1] + pe, st[2]:st[2] + pe]
        selection = "window with porosity nearest the full-image mean that percolates along x"
    full_start = [search_crop[i][0] + st[i] for i in range(3)]
    # Compare with the MHD DNS along x (volume axis 2).
    perm = properties.permeability(sub, vs, axes=(2,), mirror_flow_axis=True, tol=1e-3 if quick else 1e-4)
    ff = properties.formation_factor(sub, vs, axes=(2,), mirror_flow_axis=True)
    dest = out / name
    dest.mkdir(parents=True, exist_ok=True)
    segment.qc_figure(None, pore, dest / "qc_segmentation.png", title=f"{name}: provided segmentation")
    rve.rev_figure(rev, dest / "rev_porosity.png")
    result = {
        "name": name,
        "published": published,
        "full_image": full,
        "crop_zyx": crop,
        "crop_shape": list(pore.shape),
        "crop_porosity": float(pore.mean()),
        "morphometrics": morph,
        "rve": rev,
        "subvolume_search": search_stats,
        "properties_subvolume": {
            "start_in_image_zyx": full_start,
            "start_in_search_crop_zyx": list(st),
            "size_zyx": [pe, pe, pe],
            "porosity": float(sub.mean()),
            "axis": "x",
            "selection": selection,
        },
        "permeability": perm,
        "formation_factor": ff,
        "seconds": round(time.time() - t0, 1),
    }
    write_json(result, dest / "results.json")
    return result


def run_foram(path: Path, shape, voxel_size: float, byteorder: str, out: Path, name: str) -> dict:
    from morphoct import foram
    from morphoct.io import read_raw
    from morphoct.report import write_json

    t0 = time.time()
    grey = read_raw(path, shape, dtype="uint16", byteorder=byteorder, order="zyx")
    res = foram.analyse(grey, voxel_size, min_size=64, largest_only=True)
    dest = out / name
    dest.mkdir(parents=True, exist_ok=True)
    # mid-slice QC
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from morphoct.foram import segment_test

    shell, _, _ = segment_test(grey, min_size=64, largest_only=True)
    z = grey.shape[0] // 2
    fig, ax = plt.subplots(1, 2, figsize=(8, 4))
    ax[0].imshow(grey[z], cmap="gray")
    ax[0].set_title("grey, mid z")
    ax[1].imshow(shell[z], cmap="gray")
    ax[1].set_title("Otsu shell mask")
    for a in ax:
        a.axis("off")
    fig.tight_layout()
    fig.savefig(dest / "qc_shell.png", dpi=120)
    plt.close(fig)
    outd = {"name": name, "file": path.name, "shape_zyx": list(shape), "voxel_size_um": voxel_size,
            "byteorder": byteorder, "foram": res.as_dict(), "seconds": round(time.time() - t0, 1)}
    write_json(outd, dest / "results.json")
    return outd


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="/workspace/morphoCT_data")
    ap.add_argument("--out", default="benchmarks/out_public")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--skip-carbonate", action="store_true")
    ap.add_argument("--skip-foram", action="store_true")
    a = ap.parse_args()
    data, out = Path(a.data), Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    summary: dict = {"published": PUBLISHED}
    t0 = time.time()
    from morphoct.report import write_json

    carb = data / "carbonate" / "ImagesICPSC2015"
    if not a.skip_carbonate:
        esta = carb / "Estaillades_1000c_3p31136um.raw.gz"
        if esta.exists():
            summary["estaillades"] = run_carbonate(
                "estaillades", esta, PUBLISHED["estaillades_2015"], out,
                crop_edge=400, prop_edge=128 if a.quick else 144, quick=a.quick,
            )
            write_json(summary, out / "public_summary.json")
        ket = carb / "Ketton_1000c_3p00006um.raw.gz"
        if ket.exists():
            # Ketton: full-image porosity plus a crop for morphometrics and REV.
            # Skip LBM/FFT unless Estaillades already demonstrated the solvers (saves time).
            from morphoct import morphometrics, rve, segment

            t1 = time.time()
            pub = PUBLISHED["ketton_2015"]
            full = stream_porosity(ket, pore_value=0)
            crop = centre_crop(pub["shape"], 200 if a.quick else 400)
            arr = read_raw_gz_crop(ket, pub["shape"], crop)
            pore = (arr == 0).astype(np.uint8)
            del arr
            morph = morphometrics.compute(pore, pub["voxel_size_um"], pore_sizes=not a.quick)
            windows = [8, 12, 16, 24, 32, 48, 64] if a.quick else [8, 12, 16, 24, 32, 48, 64, 96, 128, 160]
            windows = [w for w in windows if w <= min(pore.shape) // 2]
            rev = rve.rev_porosity(pore, pub["voxel_size_um"], windows=windows, eps_rel=0.05)
            dest = out / "ketton"
            dest.mkdir(parents=True, exist_ok=True)
            segment.qc_figure(None, pore, dest / "qc_segmentation.png", title="ketton: provided segmentation")
            rve.rev_figure(rev, dest / "rev_porosity.png")
            summary["ketton"] = {
                "name": "ketton", "published": pub, "full_image": full,
                "crop_zyx": crop, "crop_shape": list(pore.shape),
                "crop_porosity": float(pore.mean()), "morphometrics": morph, "rve": rev,
                "seconds": round(time.time() - t1, 1),
            }
            write_json(summary["ketton"], dest / "results.json")
            write_json(summary, out / "public_summary.json")

    if not a.skip_foram:
        foram_dir = data / "foram_siccha2023"
        specimens = [
            {
                "file": "CA_GGL_1_207x205x187-1.2016_um-ushort_LE.raw",
                "name": "CA_GGL_1_Globigerinita_glutinata",
                "shape": [187, 205, 207], "voxel": 1.2016, "byteorder": "little",
                "taxon": "Globigerinita glutinata",
            },
            {
                "file": "GE_GTE_1_227x220x241-1.2004_um-ushort_LE.raw",
                "name": "GE_GTE_1_Globigerinoides_tenellus",
                "shape": [241, 220, 227], "voxel": 1.2004, "byteorder": "little",
                "taxon": "Globigerinoides tenellus",
            },
        ]
        summary["forams"] = []
        for sp in specimens:
            p = foram_dir / sp["file"]
            if not p.exists():
                continue
            r = run_foram(p, sp["shape"], sp["voxel"], sp["byteorder"], out, sp["name"])
            r["taxon"] = sp["taxon"]
            summary["forams"].append(r)
            write_json(summary, out / "public_summary.json")

    summary["total_seconds"] = round(time.time() - t0, 1)
    write_json(summary, out / "public_summary.json")
    print(json.dumps({"total_seconds": summary["total_seconds"],
                      "keys": [k for k in summary if k not in ("published",)]}))


if __name__ == "__main__":
    main()
