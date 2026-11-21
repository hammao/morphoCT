import json

import numpy as np
import pytest
import yaml

from morphoct import foram
from morphoct.cli import main
from morphoct.io import Volume, save_volume


@pytest.mark.slow
def test_cli_end_to_end_grey(tmp_path, blobs):
    rng = np.random.default_rng(0)
    gray = (np.where(blobs == 1, 60, 180) + rng.normal(0, 20, blobs.shape)).clip(0, 255).astype(np.uint8)
    save_volume(tmp_path / "g.tif", Volume(gray, 5.0))
    cfg = {"name": "t", "output": str(tmp_path / "out"),
           "input": {"path": str(tmp_path / "g.tif"), "voxel_size": 5.0},
           "preprocess": {"denoise": "median", "size": 3}, "segment": {"method": "otsu"},
           "rve": {"windows": [4, 6, 8, 12, 16, 24]},
           "properties": {"subvolume": 24, "axes": [0], "max_steps": 3000, "tol": 1e-3},
           "network": {"axes": [0]}}
    (tmp_path / "c.yaml").write_text(yaml.safe_dump(cfg))
    assert main(["run", str(tmp_path / "c.yaml")]) == 0
    res = json.loads((tmp_path / "out" / "results.json").read_text())
    assert abs(res["morphometrics"]["porosity"] - blobs.mean()) < 0.03
    for f in ("report.html", "report.pdf", "qc_segmentation.png", "rev_porosity.png", "binary.npz"):
        assert (tmp_path / "out" / f).exists()


def test_foram_helpers_and_stubs():
    cal = foram.linear_density_calibration([10, 200], [0.0, 2.71])
    assert abs(cal(200) - 2.71) < 1e-9
    m = np.zeros((10, 10, 10), bool)
    m[2:8, 2:8, 2:8] = True
    assert foram.shell_volume(m, 2.0) == 216 * 8
    g = np.where(m, 150.0, 0.0)
    assert foram.density_index(g, m) == 150.0
    with pytest.raises(NotImplementedError):
        foram.segment_chambers(m)
