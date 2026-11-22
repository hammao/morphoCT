"""FIB-SEM serial-section stacks: load, register, correct, then the morphoCT pipeline.

FIB-SEM volumes are typically anisotropic (xy pixel size differs from the mill
slice thickness), often in nanometres, and commonly show vertical curtaining,
charging, slow intensity drift, and in-plane slice-to-slice drift. This module
prepares a stack so the existing segmentation, morphometrics, REV and transport
solvers can be applied after optional resampling to isotropic voxels.

Methods (reimplemented from published descriptions, no vendor code copied):

* Slice registration by phase-correlation (Kuglin and Hines 1975; scikit-image
  ``phase_cross_correlation``).
* Curtaining (vertical stripe) suppression by a Fourier notch along the
  stripe axis, in the spirit of Muench, Trtik, Marone and Stampanoni 2009,
  doi:10.1364/OE.17.008567 (wavelet-Fourier stripe filter; we use a
  simpler FFT notch that targets the same vertical-line family).
* Shading / charging: subtract a strongly smoothed background (Gaussian or
  per-slice residual low-frequency field).
* Milled-face crop: bounding box of the high-variance (or high-intensity)
  mill window, which is usually a rectangular trench inside an unmilled frame.

Metadata: FEI/Thermo Atlas XML and FEI TIFF tags, generic ``.xml/.txt/.a3d``
sidecars, and ImageJ TIFF spacing. Lengths are converted to micrometres
(1 nm = 0.001 um).
"""

from __future__ import annotations

import re
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from scipy import ndimage as ndi

from morphoct.io import Volume

# FEI/Thermo and Zeiss tags often store metres or nanometres; we normalise to um
_UM = {"m": 1e6, "mm": 1e3, "um": 1.0, "µm": 1.0, "nm": 1e-3, "pm": 1e-6}


def _to_um(value: float, unit: str | None) -> float:
    if not unit:
        return float(value)
    u = unit.strip().lower().replace("μ", "u").replace("micro", "u")
    u = {"micrometer": "um", "micrometre": "um", "nanometer": "nm", "nanometre": "nm",
         "meter": "m", "metre": "m"}.get(u, u)
    if u not in _UM:
        raise ValueError(f"unknown length unit {unit!r}")
    return float(value) * _UM[u]


def _text(el) -> str:
    return (el.text or "").strip() if el is not None else ""


def _local(tag: str) -> str:
    return tag.split("}", 1)[-1] if tag else tag


def parse_metadata_xml(path) -> dict:
    """Parse an Atlas/FEI/Zeiss-style XML sidecar for voxel size and beam parameters.

    Returns a dict that may include ``voxel_size_um`` as (dz, dy, dx), plus any
    recognised milling/imaging keys (pixel size, slice thickness, HV, current,
    working distance, HFW). Missing keys are omitted, never invented.
    """
    path = Path(path)
    tree = ET.parse(path)
    root = tree.getroot()
    found: dict[str, tuple[float, str | None]] = {}
    texts: dict[str, str] = {}

    def walk(el, prefix=""):
        name = _local(el.tag)
        key = f"{prefix}/{name}" if prefix else name
        txt = _text(el)
        unit = None
        for a, v in el.attrib.items():
            an = _local(a).lower()
            if an in ("unit", "units"):
                unit = v
        if txt:
            texts[key] = txt
            try:
                found[key] = (float(txt), unit)
            except ValueError:
                pass
        for a, v in el.attrib.items():
            if _local(a).lower() in ("unit", "units"):
                continue
            try:
                found[f"{key}/@{_local(a)}"] = (float(v), unit)
            except ValueError:
                texts[f"{key}/@{_local(a)}"] = v
        for ch in el:
            walk(ch, key)

    walk(root)

    def first(*needles):
        for needle in needles:
            n = needle.lower()
            for k, (val, unit) in found.items():
                lk = k.lower().replace("_", "").replace("-", "")
                if n in lk.replace("/", ""):
                    return val, unit, k
        return None

    meta: dict = {"source": str(path), "raw_keys": {k: texts[k] for k in list(texts)[:80]}}
    pix = first("pixelsize", "xypixelsize", "pixelwidth", "xresolution", "pixel_size")
    dy = first("pixelsizey", "ypixelsize", "pixelheight")
    dx = first("pixelsizex", "xpixelsize")
    dz = first("slicethickness", "zspacing", "voxelsizez", "cutthickness", "millthickness",
               "slicespacing", "dz")
    if dx is None and pix is not None:
        dx = pix
    if dy is None and pix is not None:
        dy = pix
    vs = [None, None, None]
    if dz:
        vs[0] = _to_um(dz[0], dz[1] or "nm")
    if dy:
        vs[1] = _to_um(dy[0], dy[1] or (pix[1] if pix else "nm"))
    if dx:
        vs[2] = _to_um(dx[0], dx[1] or (pix[1] if pix else "nm"))
    if all(v is not None for v in vs):
        meta["voxel_size_um"] = tuple(vs)
    elif dx is not None and dy is not None and dz is None:
        meta["pixel_size_um"] = (_to_um(dy[0], dy[1] or "nm"), _to_um(dx[0], dx[1] or "nm"))
        meta["notes"] = ["slice thickness was not in the XML; set voxel_size dz explicitly"]
    for label, needles in {
        "high_voltage_v": ("hv", "accelerationvoltage", "beamvoltage", "highvoltage"),
        "beam_current_a": ("beamcurrent", "emissioncurrent", "probecurrent"),
        "working_distance_m": ("workingdistance", "wd"),
        "hfw_m": ("hfw", "horizontalfieldwidth", "fieldofview"),
        "dwell_time_s": ("dwelltime", "dwell"),
    }.items():
        hit = first(*needles)
        if hit:
            meta[label] = hit[0]
            if hit[1]:
                meta[label + "_unit"] = hit[1]
    return meta


def parse_metadata_text(path) -> dict:
    """Parse a simple key=value or 'key: value unit' sidecar (.txt)."""
    path = Path(path)
    meta: dict = {"source": str(path)}
    pairs = {}
    for line in path.read_text(errors="ignore").splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=" in line:
            k, v = line.split("=", 1)
        elif ":" in line:
            k, v = line.split(":", 1)
        else:
            continue
        pairs[k.strip().lower()] = v.strip()
    meta["raw_keys"] = pairs

    def take(keys, default_unit):
        for k in keys:
            if k not in pairs:
                continue
            m = re.match(r"^\s*([+-]?\d+(\.\d+)?([eE][+-]?\d+)?)\s*([a-zA-Zµμ]+)?\s*$", pairs[k])
            if not m:
                continue
            return _to_um(float(m.group(1)), m.group(4) or default_unit)
        return None

    dx = take(("pixel size x", "pixelsizex", "pixel size", "xy pixel size", "pixel_size"), "nm")
    dy = take(("pixel size y", "pixelsizey"), "nm") or dx
    dz = take(("slice thickness", "slicethickness", "z spacing", "voxel size z", "cut thickness"), "nm")
    if dx and dy and dz:
        meta["voxel_size_um"] = (dz, dy, dx)
    elif dx:
        meta["pixel_size_um"] = (dy, dx)
    return meta


def find_sidecar(stack_dir) -> Path | None:
    """Return the first metadata sidecar in or next to ``stack_dir``."""
    d = Path(stack_dir)
    if d.is_file():
        d = d.parent
    for pat in ("*.xml", "*.a3d", "*.txt", "*.hdr", "*metadata*"):
        hits = sorted(d.glob(pat)) + sorted(d.parent.glob(pat)) if d.parent else sorted(d.glob(pat))
        for h in hits:
            if h.is_file() and h.suffix.lower() in {".xml", ".a3d", ".txt", ".hdr"}:
                return h
    return None


def load_metadata(path) -> dict:
    """Load voxel-size metadata from a sidecar or a directory that contains one."""
    path = Path(path)
    if path.is_dir():
        side = find_sidecar(path)
        if side is None:
            return {"source": str(path), "notes": ["no sidecar found"]}
        path = side
    suf = path.suffix.lower()
    if suf in {".xml", ".a3d"}:
        return parse_metadata_xml(path)
    return parse_metadata_text(path)


def list_stack_frames(path) -> list[Path]:
    """Sorted frame paths of a TIFF/PNG/BMP directory, a glob, or a single multi-page TIFF."""
    p = Path(path)
    if p.is_file():
        return [p]
    if any(c in str(path) for c in "*?["):
        import glob as _glob
        return [Path(f) for f in sorted(_glob.glob(str(path)))]
    exts = {".tif", ".tiff", ".png", ".bmp", ".jpg", ".jpeg"}
    files = sorted(f for f in p.iterdir() if f.is_file() and f.suffix.lower() in exts)
    if not files:
        raise FileNotFoundError(f"no image frames in {path}")
    return files


def load_stack(path, *, voxel_size=None, metadata=None, max_frames=None) -> Volume:
    """Load a FIB-SEM frame stack as a ``[z, y, x]`` Volume.

    ``voxel_size`` is (dz, dy, dx) in micrometres, or a scalar. When omitted, a
    sidecar in the folder is parsed; if that still lacks a z-spacing the loader
    records a warning and uses 1 um placeholders only if nothing else is known.
    """
    import tifffile
    from PIL import Image

    frames = list_stack_frames(path)
    if max_frames:
        frames = frames[: int(max_frames)]
    meta = {"source": str(path), "n_frames": len(frames), "frame_names": [f.name for f in frames[:8]]}
    if metadata:
        meta.update(metadata)
    else:
        side = load_metadata(path)
        meta.update({k: v for k, v in side.items() if k != "raw_keys"})
        if "raw_keys" in side:
            meta["sidecar"] = side.get("source")
    planes = []
    for f in frames:
        if f.suffix.lower() in {".tif", ".tiff"}:
            arr = tifffile.imread(f)
        else:
            arr = np.asarray(Image.open(f))
        if arr.ndim == 3 and arr.shape[-1] in (3, 4):
            arr = np.mean(arr[..., :3], axis=-1)
        if arr.ndim == 3:
            # multi-page single file: take as a stack
            if len(frames) == 1:
                planes = [np.asarray(s) for s in (arr if arr.ndim == 3 else [arr])]
                break
            arr = arr[0]
        planes.append(np.asarray(arr))
    data = np.stack(planes, axis=0)
    vs = voxel_size or meta.get("voxel_size_um")
    if vs is None:
        pix = meta.get("pixel_size_um")
        if pix is not None:
            vs = (1.0, float(pix[0]), float(pix[1]))
            meta["voxel_size_warning"] = "slice thickness unknown; dz set to 1 um as a placeholder"
        else:
            vs = (1.0, 1.0, 1.0)
            meta["voxel_size_warning"] = "voxel size unknown, set voxel_size=(dz, dy, dx) in micrometres"
    return Volume(np.ascontiguousarray(data), vs, meta)


def register_slices(vol: Volume, *, upsample=10, ref="previous") -> tuple[Volume, np.ndarray]:
    """Align slices in-plane by phase correlation; return (aligned volume, shifts).

    ``shifts[z] = (dy, dx)`` applied to bring slice z onto the reference
    (the previous already-aligned slice, or the middle slice when
    ``ref='middle'``). Empty margins after the shift are filled with the
    slice median. A common overlap crop is applied so every slice is valid.
    """
    from skimage.registration import phase_cross_correlation
    from skimage.transform import warp, AffineTransform

    d = np.asarray(vol.data)
    nz, ny, nx = d.shape
    shifts = np.zeros((nz, 2), float)
    out = np.empty_like(d, dtype=np.float32)
    if ref == "middle":
        refs = d[nz // 2]
        out[nz // 2] = d[nz // 2]
        order = list(range(nz // 2 - 1, -1, -1)) + list(range(nz // 2 + 1, nz))
        for z in order:
            shift, _, _ = phase_cross_correlation(refs, d[z], upsample_factor=upsample)
            shifts[z] = -np.asarray(shift, float)
    else:
        out[0] = d[0]
        acc = np.zeros(2, float)
        for z in range(1, nz):
            # Incremental shift between adjacent original slices, then accumulate.
            inc, _, _ = phase_cross_correlation(d[z - 1], d[z], upsample_factor=upsample)
            acc = acc - np.asarray(inc, float)
            shifts[z] = acc
            tform = AffineTransform(translation=(acc[1], acc[0]))  # x, y
            out[z] = warp(d[z].astype(np.float32), tform, order=1,
                          mode="constant", cval=float(np.median(d[z])), preserve_range=True)
        # the first slice is already in out; later slices written in the loop
        # (z=0 is the reference)
        aligned = out
        # crop to common overlap
        # positive dy means the moving image was down relative to the ref, etc.
        y0 = int(np.ceil(max(0, shifts[:, 0].max())))
        y1 = int(ny + np.floor(min(0, shifts[:, 0].min())))
        x0 = int(np.ceil(max(0, shifts[:, 1].max())))
        x1 = int(nx + np.floor(min(0, shifts[:, 1].min())))
        if y1 - y0 < 8 or x1 - x0 < 8:
            y0, y1, x0, x1 = 0, ny, 0, nx
        cropped = aligned[:, y0:y1, x0:x1]
        meta = {**vol.meta, "registration": "phase_cross_correlation",
                "shifts_yx": shifts.tolist(), "overlap_crop": [y0, y1, x0, x1]}
        return Volume(np.ascontiguousarray(cropped), vol.voxel_size, meta), shifts

    # middle-ref path
    for z in range(nz):
        if z == nz // 2:
            continue
        tform = AffineTransform(translation=(shifts[z, 1], shifts[z, 0]))
        out[z] = warp(d[z].astype(np.float32), tform, order=1,
                      mode="constant", cval=float(np.median(d[z])), preserve_range=True)
    y0 = int(np.ceil(max(0, shifts[:, 0].max())))
    y1 = int(ny + np.floor(min(0, shifts[:, 0].min())))
    x0 = int(np.ceil(max(0, shifts[:, 1].max())))
    x1 = int(nx + np.floor(min(0, shifts[:, 1].min())))
    if y1 - y0 < 8 or x1 - x0 < 8:
        y0, y1, x0, x1 = 0, ny, 0, nx
    cropped = out[:, y0:y1, x0:x1]
    meta = {**vol.meta, "registration": "phase_cross_correlation",
            "shifts_yx": shifts.tolist(), "overlap_crop": [y0, y1, x0, x1]}
    return Volume(np.ascontiguousarray(cropped), vol.voxel_size, meta), shifts


def correct_curtaining(vol: Volume, *, sigma=0, width=3) -> Volume:
    """Attenuate vertical curtaining (mill stripes) with a Fourier notch.

    Vertical stripes concentrate energy at fx = 0 and small fx for all fy
    (Muench et al. 2009). A narrow Gaussian notch is applied around fx = 0,
    leaving the DC term intact so mean intensity is preserved. ``width`` is the
    half-width of the notch in frequency bins.
    """
    d = np.asarray(vol.data, dtype=np.float32)
    out = np.empty_like(d)
    ny, nx = d.shape[1], d.shape[2]
    # Vertical stripes (constant along y) live at fy ~ 0, fx != 0.
    fy_w = np.exp(-0.5 * (np.fft.fftfreq(ny) * ny / max(float(width), 0.5)) ** 2)
    notch = (1.0 - fy_w)[:, None] * np.ones((1, nx), np.float32)
    notch[0, 0] = 1.0
    for z in range(d.shape[0]):
        F = np.fft.fft2(d[z])
        out[z] = np.fft.ifft2(F * notch).real
    return vol.with_data(out, curtaining="fft_notch", curtaining_width=width)


def correct_shading(vol: Volume, *, sigma=32.0) -> Volume:
    """Remove a slowly varying shading / charging field on each slice.

    A Gaussian of ``sigma`` pixels is subtracted, then the slice mean is
    restored so global grey levels stay comparable for thresholding.
    """
    d = np.asarray(vol.data, dtype=np.float32)
    out = np.empty_like(d)
    for z in range(d.shape[0]):
        sl = d[z]
        bg = ndi.gaussian_filter(sl, sigma)
        adj = sl - bg
        out[z] = adj + (sl.mean() - adj.mean())
    return vol.with_data(out, shading="gaussian_residual", shading_sigma=sigma)


def crop_milled_region(vol: Volume, *, pad=4, min_frac=0.15) -> Volume:
    """Crop to the milled window, estimated from per-pixel variance along z.

    Unmilled frame is usually flatter across slices than the trench. The
    bounding box of voxels whose z-variance exceeds a high quantile is kept,
    padded by ``pad`` voxels. Falls back to the full field if the box would
    cover less than ``min_frac`` of the area.
    """
    d = np.asarray(vol.data, dtype=np.float32)
    if d.shape[0] < 3:
        return vol
    var = d.var(axis=0)
    try:
        from skimage.filters import threshold_otsu
        thr = float(threshold_otsu(var))
    except Exception:
        thr = float(np.quantile(var, 0.75))
    thr = max(thr, float(np.quantile(var, 0.7)))
    mask = var >= thr
    mask = ndi.binary_opening(mask, iterations=1)
    mask = ndi.binary_closing(mask, iterations=1)
    if not np.any(mask):
        return vol
    lab, n = ndi.label(mask)
    if n:
        counts = np.bincount(lab.ravel())
        counts[0] = 0
        mask = lab == int(np.argmax(counts))
    ys, xs = np.where(mask)
    y0, y1 = int(ys.min()) - pad, int(ys.max()) + 1 + pad
    x0, x1 = int(xs.min()) - pad, int(xs.max()) + 1 + pad
    y0, x0 = max(y0, 0), max(x0, 0)
    y1, x1 = min(y1, d.shape[1]), min(x1, d.shape[2])
    area = (y1 - y0) * (x1 - x0)
    if area < min_frac * d.shape[1] * d.shape[2]:
        return vol
    cropped = d[:, y0:y1, x0:x1]
    return vol.with_data(np.ascontiguousarray(cropped),
                         milled_crop=[y0, y1, x0, x1])


def resample_isotropic(vol: Volume, *, target_um: float | None = None) -> Volume:
    """Zoom the volume so dz = dy = dx = ``target_um`` (default: the smallest spacing)."""
    dz, dy, dx = vol.voxel_size
    target = float(target_um if target_um is not None else min(dz, dy, dx))
    zoom = (dz / target, dy / target, dx / target)
    if np.allclose(zoom, 1.0, rtol=1e-4):
        return vol.with_data(vol.data, resampled_isotropic_um=target)
    out = ndi.zoom(vol.data.astype(np.float32), zoom, order=1)
    return Volume(np.ascontiguousarray(out), (target, target, target),
                  {**vol.meta, "resampled_isotropic_um": target, "zoom": zoom})


def prepare(vol: Volume, *, register=True, curtaining=True, shading=True,
            crop_milled=True, isotropic=True, target_um=None,
            curtaining_width=3, shading_sigma=32.0) -> Volume:
    """Apply the usual FIB-SEM corrections in a fixed order and return a Volume."""
    if register and vol.shape[0] >= 2:
        vol, _ = register_slices(vol)
    if curtaining:
        vol = correct_curtaining(vol, width=curtaining_width)
    if shading:
        vol = correct_shading(vol, sigma=shading_sigma)
    if crop_milled and vol.shape[0] >= 3:
        vol = crop_milled_region(vol)
    if isotropic:
        vol = resample_isotropic(vol, target_um=target_um)
    return vol


def run_config(config: dict | str | Path) -> dict:
    """Prepare a FIB-SEM stack from YAML and run :func:`morphoct.pipeline.run`.

    Config keys (in addition to the usual pipeline keys)::

        input: {path: /path/to/frames, voxel_size: [0.01, 0.005, 0.005]}  # dz, dy, dx in um
        fibsem:
          register: true
          curtaining: true
          shading: true
          crop_milled: true
          isotropic: true
          max_frames: null
        # then segment / rve / properties as in morphoct run
    """
    import yaml

    from morphoct.pipeline import run

    if not isinstance(config, dict):
        config = yaml.safe_load(Path(config).read_text())
    cfg = dict(config)
    fcfg = dict(cfg.pop("fibsem", {}) or {})
    inp = dict(cfg["input"])
    path = inp.pop("path")
    max_frames = fcfg.pop("max_frames", None)
    vol = load_stack(path, voxel_size=inp.pop("voxel_size", None), max_frames=max_frames)
    # remaining input keys (crop etc.) are not used here
    vol = prepare(vol, **{k: v for k, v in fcfg.items() if k in
                          {"register", "curtaining", "shading", "crop_milled",
                           "isotropic", "target_um", "curtaining_width", "shading_sigma"}})
    out = Path(cfg.get("output", f"morphoct_out/{cfg.get('name', 'fibsem')}"))
    out.mkdir(parents=True, exist_ok=True)
    tmp = out / "prepared.npz"
    from morphoct.io import save_volume
    save_volume(tmp, vol)
    cfg = dict(cfg)
    cfg["input"] = {"path": str(tmp), "voxel_size": list(vol.voxel_size)}
    cfg.setdefault("output", str(out))
    # prepared volume is grey; keep preprocess/segment from config
    result = run(cfg)
    result["fibsem"] = {
        "shape_zyx": list(vol.shape),
        "voxel_size_um": list(vol.voxel_size),
        "meta": {k: vol.meta[k] for k in vol.meta if k not in ("frame_names", "raw_keys", "shifts_yx")},
    }
    return result
