# morphoCT

Scientifically grounded microCT analysis of **rocks, core plugs and foraminifera tests**,
built on [morphoGen](https://github.com/jdirr/morphoGen) (Dirrenberger, GPL-3.0-or-later)
and [PoreSpy](https://github.com/PMEAL/porespy) / [OpenPNM](https://github.com/PMEAL/OpenPNM) (MIT).

Status: **0.1.0.dev0, Phase 1 MVP** (rock and core). The foraminifera module is an API stub.
Licence: **GPL-3.0-or-later** (see `LICENSE`), as required by the morphoGen dependency.

## Pipeline

```
import -> preprocess -> segment (+QC) -> morphometrics -> REV -> properties -> network -> report
```

| Module | What it does | Built on |
|---|---|---|
| `morphoct.io` | TIFF stack (multi-page or a folder of slices), RAW (shape, dtype, byte order, voxel size; `.gz` ok), NumPy `.npy/.npz`, NIfTI/NRRD/MHA via SimpleITK (optional extra `itk`) | tifffile, SimpleITK |
| `morphoct.preprocess` | crop, median/Gaussian/non-local-means denoise, block downsample, **mirroring** for periodic solvers | SciPy, scikit-image |
| `morphoct.segment` | Otsu, multi-Otsu, small-object cleaning, marker watershed of grains, threshold-sensitivity sweep, Dice, QC figure | scikit-image |
| `morphoct.morphometrics` | porosity, connected porosity and percolation, Minkowski densities (V_V, S_V, M_V, chi_V), two-point correlation and correlation length, chord lengths with stereological check, local-thickness pore sizes | morphoGen `analysis`, PoreSpy |
| `morphoct.rve` | window-method REV of porosity (Kanit et al. 2003), fit of D^2(V) = K V^-gamma, REV size with a 95 % CI | morphoGen `rve` |
| `morphoct.properties` | permeability by lattice Boltzmann (TRT D3Q19), formation factor and Archie m by FFT conduction, Kozeny-Carman check | morphoGen `homogenization` |
| `morphoct.network` | SNOW2 pore network, OpenPNM Stokes flow permeability | PoreSpy, OpenPNM |
| `morphoct.report` | `results.json`, one-page `report.html` (figures embedded) and `report.pdf` | matplotlib |
| `morphoct.foram` | planned API: test segmentation, interior separation, chamber watershed, wall thickness, density index (generic helpers implemented) | (Phase 3) |

## Install and run

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e ".[itk,test]"
morphoct run examples/fontainebleau.yaml
pytest -q              # fast tests
pytest -q -m slow      # end-to-end and network tests
```

A config (`examples/*.yaml`) names the input, whether it is already segmented
(`binary: {pore_value: 0}`) or grey (then `preprocess` and `segment` run), the REV
windows, the voxel-property subvolume and axes, and the network options. Outputs go to
`output:`.

### Units and axes
Volumes are indexed `[z, y, x]` (slice order); voxel sizes are micrometres. Permeability
is reported in mD (1 mD = 9.869233e-4 um^2), lengths in um, S_V per um.

## Benchmark (Phase 1)

`benchmarks/run_benchmarks.py` downloads the public Andrä et al. (2013) digital rock
physics benchmark images from https://github.com/fkrzikalla/drp-benchmarks and compares
morphoCT with the published values (Andrä et al. 2013, Part I, doi:10.1016/j.cageo.2012.09.005;
Part II, doi:10.1016/j.cageo.2012.09.008). Results: `benchmarks/RESULTS.md`.
The benchmark images are not redistributed here.

## Where a morphoGen fork patch would be needed

morphoCT pins `morphogen==1.0.2` from PyPI and does not vendor its code, so a fork is
optional. Patches worth upstreaming (or carrying in `hammao/morphoGen`) are:

1. **Non-periodic solvers**: `homogenization.flow` and the FFT solver are periodic only.
   morphoCT mirrors along the flow axis; a pressure-driven inlet/outlet LBM boundary
   (or buffer layers) in `flow.py` would remove the wrap of transverse faces.
2. **LBM convergence on large rock images**: `stokes_flow` defaults to
   `max_steps = 200 * n`; a multigrid start or a GPU path would cut runtime on 300^3+
   images.
3. **Readers**: `morphogen.io` lacks RAW/NIfTI/NRRD (morphoCT provides them).
4. **Model identification**: fitting morphoGen's random-set models to measured S2 and
   chord functions (planned morphoCT Phase 4) may belong upstream in `morphogen.analysis`.
5. **Units metadata**: `VoxelImage.metadata` could carry the voxel unit explicitly.

## Scientific basis (selected, DOIs verified)

* Kanit et al. 2003, RVE size: doi:10.1016/S0020-7683(03)00143-4
* Moulinec and Suquet 1998, FFT homogenisation: doi:10.1016/S0045-7825(97)00218-1
* Blunt et al. 2013, pore-scale imaging and modelling: doi:10.1016/j.advwatres.2012.03.003
* Vogel et al. 2010, Minkowski functions of soil structure: doi:10.1016/j.cageo.2010.03.007
* Gostick 2017, SNOW network extraction: doi:10.1103/PhysRevE.96.023307
* Gostick et al. 2019, PoreSpy: doi:10.21105/joss.01296
* Archie 1942, formation factor: doi:10.2118/942054-G

## Licences of what morphoCT uses

morphoGen GPL-3.0-or-later; PoreSpy MIT; OpenPNM MIT; scikit-image BSD-3; SimpleITK
Apache-2.0; tifffile BSD-3; NumPy/SciPy BSD-3; matplotlib PSF-style. All are compatible
with GPL-3.0-or-later. No code from ForamJ (no licence file) is included.
