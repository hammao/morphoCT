"""One-page HTML and PDF report plus a JSON of every result."""

from __future__ import annotations

import base64
import datetime as _dt
import html
import json
from pathlib import Path

import numpy as np


def _jsonable(o):
    if isinstance(o, dict):
        return {str(k): _jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [_jsonable(v) for v in o]
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return None if not np.isfinite(f) else f
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return _jsonable(o.tolist())
    return o


def write_json(results: dict, path) -> Path:
    path = Path(path)
    path.write_text(json.dumps(_jsonable(results), indent=2))
    return path


def _fmt(v, nd=3):
    if v is None:
        return "n/a"
    if isinstance(v, float):
        if not np.isfinite(v):
            return "n/a"
        return f"{v:.{nd}g}"
    return html.escape(str(v))


def summary_rows(r: dict) -> list[tuple[str, str]]:
    m = r.get("morphometrics", {})
    rows = [("Sample", r.get("name", "")), ("Voxel size (um)", _fmt(m.get("voxel_size_um"))),
            ("Shape (z, y, x)", str(m.get("shape_zyx", ""))), ("Porosity", _fmt(m.get("porosity"), 4))]
    if "connectivity" in m:
        rows.append(("Connected porosity", _fmt(m["connectivity"]["connected_porosity"], 4)))
    if "minkowski" in m:
        mk = m["minkowski"]
        rows += [("S_V (1/um)", _fmt(mk["S_V_per_um"])), ("M_V (1/um^2)", _fmt(mk["M_V_per_um2"])),
                 ("chi_V (1/um^3)", _fmt(mk["chi_V_per_um3"]))]
    if "mean_chord_length_pooled_um" in m:
        rows.append(("Mean pore chord (um)", _fmt(m["mean_chord_length_pooled_um"])))
    if "correlation_length_um" in m:
        rows.append(("Correlation length (um)", _fmt(m["correlation_length_um"])))
    if "pore_size_um" in m:
        p = m["pore_size_um"]
        rows.append(("Pore size d10/d50/d90 (um)", f"{_fmt(p['d10'])} / {_fmt(p['d50'])} / {_fmt(p['d90'])}"))
    rv = r.get("rve")
    if rv:
        ci = rv["rev_edge_um_ci95"]
        rows.append((f"REV edge, porosity +/-{rv['eps_rel']:.0%} (um)", f"{_fmt(rv['rev_edge_um'])} (95% CI {_fmt(ci[0])} to {_fmt(ci[1])})"))
        rows.append(("Sample larger than REV", str(rv["sample_is_rev"])))
    k = r.get("permeability")
    if k:
        for ax, v in k["axes"].items():
            rows.append((f"Permeability k_{ax}, LBM (mD)", _fmt(v["k_mD"], 4)))
    f = r.get("formation_factor")
    if f:
        for ax, v in f["axes"].items():
            rows.append((f"Formation factor F_{ax}", _fmt(v["F"], 4)))
        rows.append(("Archie m", _fmt(f["archie_m"])))
    n = r.get("network")
    if n:
        rows.append(("Network pores / throats", f"{n['n_pores']} / {n['n_throats']}"))
        for ax, v in n["axes"].items():
            rows.append((f"Network k_{ax} (mD)", _fmt(v["k_mD"], 4)))
    seg = r.get("segmentation")
    if seg and "threshold" in seg:
        rows.append(("Segmentation", f"{seg['method']} threshold {_fmt(seg['threshold'])}"))
    return rows


def write_html(results: dict, figures: list, path) -> Path:
    path = Path(path)
    rows = "".join(f"<tr><th>{html.escape(a)}</th><td>{b}</td></tr>" for a, b in summary_rows(results))
    imgs = ""
    for f in figures:
        f = Path(f)
        if f.exists():
            b64 = base64.b64encode(f.read_bytes()).decode()
            imgs += f'<figure><img src="data:image/png;base64,{b64}"/><figcaption>{html.escape(f.stem)}</figcaption></figure>'
    notes = "".join(f"<li>{html.escape(n)}</li>" for n in results.get("notes", []))
    doc = f"""<!doctype html><html><head><meta charset="utf-8"><title>morphoCT report: {html.escape(results.get('name',''))}</title>
<style>body{{font-family:system-ui,sans-serif;max-width:1100px;margin:24px auto;color:#222}}
table{{border-collapse:collapse}}th,td{{border-bottom:1px solid #ddd;padding:4px 10px;text-align:left}}
th{{font-weight:500;color:#555}}figure{{display:inline-block;margin:6px}}img{{max-width:520px}}</style></head><body>
<h1>morphoCT report: {html.escape(results.get('name',''))}</h1>
<p>morphoCT {html.escape(results.get('versions',{}).get('morphoct',''))}, morphoGen {html.escape(results.get('versions',{}).get('morphogen',''))},
PoreSpy {html.escape(results.get('versions',{}).get('porespy',''))}. Generated {_dt.datetime.now().isoformat(timespec='minutes')}.</p>
<table>{rows}</table><ul>{notes}</ul>{imgs}
<p>Full results: results.json. Licence of morphoCT: GPL-3.0-or-later.</p></body></html>"""
    path.write_text(doc)
    return path


def write_pdf(results: dict, figures: list, path) -> Path:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    rows = summary_rows(results)
    fig = plt.figure(figsize=(8.27, 11.69))
    fig.text(0.05, 0.965, f"morphoCT report: {results.get('name','')}", fontsize=14, weight="bold")
    ax = fig.add_axes([0.05, 0.55, 0.9, 0.39])
    ax.axis("off")
    tab = ax.table(cellText=[[a, b.replace("&#x27;", "'")] for a, b in rows], colWidths=[0.5, 0.5], loc="upper left", cellLoc="left")
    tab.auto_set_font_size(False)
    tab.set_fontsize(7.5)
    tab.scale(1, 1.15)
    figs = [Path(f) for f in figures if Path(f).exists()][:4]
    for i, f in enumerate(figs):
        a = fig.add_axes([0.03 + 0.48 * (i % 2), 0.28 - 0.26 * (i // 2), 0.46, 0.25])
        a.imshow(mpimg.imread(f))
        a.axis("off")
    fig.savefig(path)
    plt.close(fig)
    return Path(path)


def s2_figure(m: dict, path):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    t = m["two_point_correlation"]
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    ax.plot(t["r_um"], t["S2"], "-")
    ax.axhline(m["porosity"] ** 2, ls="--", color="0.5", label="phi^2")
    ax.set_xlabel("lag r (um)")
    ax.set_ylabel("S2(r)")
    ax.set_title("Two-point correlation of the pores")
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)
    return path
