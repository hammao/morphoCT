# Phase 1 benchmark: Andrä et al. (2013) digital rock physics benchmarks

Run on the box (8 CPU cores, 15 GB RAM), 21 Nov 2026 (UTC+3), morphoCT 0.1.0.dev0,
morphoGen 1.0.2, PoreSpy 3.3.0, OpenPNM 3.6.5. Data: https://github.com/fkrzikalla/drp-benchmarks.
Published values: `results/simulations.ods` of that repository (Andrä et al. 2013 Part II,
doi:10.1016/j.cageo.2012.09.008; segmentations from Part I, doi:10.1016/j.cageo.2012.09.005).
Raw outputs: `benchmarks/results/`.

## Fontainebleau sandstone (ExxonMobil segmentation, 288 x 288 x 300 voxels, 7.5 um)

| Quantity | morphoCT | Published | Note |
|---|---|---|---|
| Porosity | 0.1467 | 0.147 | full image |
| Permeability k_z (mD) | 1785 | LBM z 1624 (Kongju); Stokes z 2073 (Fraunhofer); means 1503 (LBM) and 1914 (Stokes) | morphoCT: LBM on the centred 144^3 subvolume (porosity 0.1481), mirrored along z; published: full image |
| Formation factor F_z | 27.0 | 28.1 (Kongju FEM), 38.7 (Fraunhofer EJ); means 29.4 to 41.6 | FFT on the same mirrored subvolume; 4 % below the lowest published z value |
| Archie m | 1.73 | | from F_z and subvolume porosity |
| Network k (mD), z / y / x | 339 / 286 / 323 | as above | SNOW2 + OpenPNM, spheres and cylinders, inscribed diameters: about 5x too low |
| Network k_z range over 12 conduit choices (mD) | 339 to 4295 | | the conduit model, not the image, dominates network permeability |
| Kozeny-Carman (c = 5) (mD) | 5975 | | sanity bound only |
| Connected porosity | 0.1463 | | 99.8 % of pores in spanning clusters |
| S_V (1/um) | 0.0103 | | Minkowski density, per bulk volume |
| Mean pore chord (um) | 54.5 (stereology 4 phi/S_V: 56.7) | | |
| REV edge for porosity +/-5 % (um) | 2090 (95 % CI 1364 to 3276) | | sample edge 2190 um: just at the REV |
| Runtime (s) | 1144 total: LBM 751, FFT 336, network 47, morphometrics + REV 8 | | a standalone rerun of LBM + FFT on the same subvolume took 534 + 75 s |

## Berea sandstone (VSG segmentation, benchmark subvolume x, y 151 to 550, z 126 to 525; 400^3 at 0.74 um)

| Quantity | morphoCT | Published | Note |
|---|---|---|---|
| Porosity | 0.1839 | 0.184 | same subvolume |
| Network k (mD), z / y / x | 9.7 / 10.0 / 13.1 | Stokes x 123.6, LBM x 113.4; means 105.6 and 95.5 | cones and cylinders (spheres failed: pores overlap); network about 10x too low |
| REV edge for porosity +/-5 % (um) | 2015 (95 % CI 819 to 5804) | | sample edge 296 um: far below the REV (gamma 0.57) |
| LBM permeability | not run | | 400^3 too slow for the CPU solver in this session |

## Segmentation from grey scale (Berea slices 126 to 205, same x, y window, 80 x 400 x 400)

| Method | Porosity | VSG porosity, same voxels | Dice vs VSG | Time (s) |
|---|---|---|---|---|
| median 3 + Otsu | 0.1773 | 0.1744 | 0.947 | 5 |
| non-local means + Otsu | 0.1812 | 0.1744 | 0.913 | 131 |

Threshold sensitivity (median + Otsu): shifting the threshold by -5, -2.5, +2.5, +5 % of
the grey range gives porosity 0.143, 0.168, 0.187, 0.199. Published porosities of three
segmentations of the same Berea image: 0.184 to 0.209 (Part I).

## Reading the results

* Voxel solvers (LBM, FFT) agree with the published voxel solvers within their own
  spread for Fontainebleau, after mirroring.
* Pore network permeability is not trustworthy without calibration of the conduit model;
  it is kept as a fast screening tool.
* Neither benchmark image is clearly larger than its porosity REV, so single-image
  property values carry a sampling uncertainty that the published table does not show.
