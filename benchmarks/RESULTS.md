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

## Public test data: Imperial College 2015 carbonates and Siccha et al. 2023 forams

Run on the box (8 CPU cores, 15 GB RAM), 21 Nov 2026 (Asia/Riyadh, UTC+3), morphoCT
0.1.0.dev0, morphoGen 1.0.2, PoreSpy 3.3.0. Volumes live in `/workspace/morphoCT_data`
(outside the repository; see that folder's README for provenance). Raw outputs:
`benchmarks/results/public/` and `benchmarks/out_public/`. The images are not
redistributed here.

### Datasets

| Dataset | Where | Licence | Citation |
|---|---|---|---|
| Estaillades and Ketton, 1000^3 uint8, 2015 release | Catalogue: https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/micro-ct-images-and-networks/ . Copy used: public Google Drive archive `ImagesICPSC2015.zip` (file id `1O_CCcoMxqdC6ebDOL-2DCM7jM__FHQEd`), the same four-rock bundle linked from the public `drd` loader; no account required (only the standard Drive virus-scan confirm). | No SPDX licence file in the archive. The Imperial College Pore-Scale Modelling group distributes the images as freely available research data. Related Imperial figshare deposits ask users to acknowledge the Qatar Carbonates and Carbon Storage Research Centre, Imperial College London. | Imperial College Pore-Scale Modelling group, Micro-CT images and networks, 2015 release (`Estaillades_1000c_3p31136um`, `Ketton_1000c_3p00006um`). |
| Planktic foraminifera microCT, Siccha et al. 2023 | Individual raw files (no login): https://download.pangaea.de/dataset/949585/files/<filename> . Table: https://doi.pangaea.de/10.1594/PANGAEA.949585?format=textfile . The all-files zip/tar requires a PANGAEA account and was not used. | CC-BY-4.0 | Siccha, M., Morard, R., Meilland, J., Iwasaki, S., Kucera, M. and Kimoto, K. (2023). Collection of X-ray micro computed tomography images of shells of planktic foraminifera with curated taxonomy. PANGAEA. https://doi.org/10.1594/PANGAEA.949585 . Paper: Siccha et al. (2023), Scientific Data 10:679, https://doi.org/10.1038/s41597-023-02498-0 |

Found but not downloaded (this host was blocked or the source needs a login): Shah, Crawshaw and Boek figshare carbonates (CC BY 4.0), doi:10.6084/m9.figshare.1230047 (Estaillades), 1230044 (Ketton), 1230046 (Portland), HTTP 403; Imperial 2009 C1/C2 on figshare, doi:10.6084/m9.figshare.1189257 and 1189258, also 403 (Dong 2007, Table 7.3, gives C1 porosity 23.3 %, mean K 1102 mD, mean FF 26.78 and C2 porosity 16.8 %, mean K 72.333 mD, mean FF 169.78); Digital Rocks Portal Estaillades, doi:10.17612/P73W2C, now redirects to digitalporousmedia.org and requires a login; Imperial Box shares from the catalogue pages returned HTTP 404 without a browser session.

### Estaillades carbonate (1000^3, 3.31136 um, 0 = pore, 1 = solid)

Published transport values are the DNS block `Estaillades_1000c_3p31136um-1-X` in the dataset `.mhd` header: effective porosity 0.109066, K_x = 2.18531e-13 m^2 (221.4 mD), FF_x = 166.697, on the full 1000^3 image. The header does not report total (including isolated) image porosity. morphoCT morphometrics and the porosity REV were computed on a centred 400^3 crop; permeability and formation factor on a 144^3 window, mirrored along x.

| Quantity | morphoCT | Published (MHD DNS, full image) | Note |
|---|---|---|---|
| Total porosity (full 1000^3) | 0.1273 | (not given; effective 0.109066) | streamed over the raw.gz; values are only 0 and 1 |
| Connected porosity (400^3 crop) | 0.1194 | effective 0.109066 | crop porosity 0.1310; 91.2 % of crop pores in spanning clusters, percolating in z, y and x |
| Permeability k_x (mD) | 152.4 | 221.4 | LBM on a 144^3 window at image index z,y,x = (404, 564, 308), porosity 0.1273, mirrored along x; 27250 steps, residual 9.7e-5, 921 s. 31 % below the full-image DNS |
| Formation factor F_x | 116.7 | 166.7 | FFT on the same mirrored window; 30 % below the DNS. Archie m = 2.31 |
| S_V (1/um) | 0.0104 | | Minkowski density, 400^3 crop |
| Mean pore chord (um) | 48.8 (stereology 4 phi/S_V: 50.5) | | crop |
| Pore diameter d10 / d50 / d90 (um) | 6.6 / 19.9 / 43.0 | | local thickness, crop |
| Correlation length (um) | 53.0 | | crop |
| REV edge for porosity +/-5 % (um) | 301666 (95 % CI 70270 to 1.75e6) | | fit on the 400^3 crop (edge 1325 um). gamma = 0.34, r^2 = 0.979. The crop, and the 3.31 mm full image, are far below this REV |
| Runtime (s) | 942 first pass + 921 LBM on the matched window | | |

A naïve centred 144^3 inside the 400^3 crop landed in a tight region (porosity 0.057). LBM there gave k_x = 1.54 mD and did not reach the residual; the FFT conductivity was non-positive, so F was not reported. A grid of 1728 windows of 144^3 (step 32) in a centred 512^3 crop has porosity 0.016 to 0.395 (mean 0.130, std 0.076). The window used above is the one whose porosity is nearest the full-image mean and that percolates along x. Even then, a 144^3 subvolume cannot be expected to reproduce full-image DNS on a carbonate that is far from a porosity REV; the similar ~30 % shortfall in both k_x and F_x is consistent with missing long-range conductive paths.

### Ketton carbonate (1000^3, 3.00006 um, 0 = pore, 1 = solid)

The Ketton `.mhd` header has no bundled DNS, so only image porosity and crop morphometrics are reported. LBM and FFT were not run (Estaillades already exercised the solvers).

| Quantity | morphoCT | Published with this file | Note |
|---|---|---|---|
| Total porosity (full 1000^3) | 0.1330 | (none in the header) | this is resolved (inter-oolith) image porosity; laboratory helium porosity of Ketton is commonly higher because of intra-oolith microporosity not resolved at 3 um |
| Connected porosity (400^3 crop) | 0.1404 | | crop porosity 0.1411; 99.5 % of crop pores in spanning clusters |
| S_V (1/um) | 0.00622 | | crop |
| Mean pore chord (um) | 83.1 (stereology 90.7) | | crop |
| Pore diameter d10 / d50 / d90 (um) | 15.0 / 42.0 / 78.0 | | local thickness, crop |
| Correlation length (um) | 75.0 | | crop |
| REV edge for porosity +/-5 % (um) | 17558 (95 % CI 1887 to 5.13e5) | | crop edge 1200 um, full image 3000 um; gamma = 0.52, r^2 = 0.916; the sample is not a porosity REV |
| Runtime (s) | 26 | | |

### Foraminifera (Siccha et al. 2023)

Baseline Otsu shell segmentation (bright calcite), small-object cleaning, largest connected component, then calcite volume, marching-cubes surface area, a morphological closing of radius 8 voxels followed by hole filling (to seal apertures), and the convex-hull packing fraction. Mean grey of an eroded shell is reported as an uncalibrated density index (no density standards were scanned with these specimens). Chamber labelling is not done. The dataset paper does not publish shell volumes for these specimen IDs, so there is no numerical comparison table.

| Specimen | Shell volume (um^3) | Filled test (um^3) | Test porosity | Shell / hull | Surface (um^2) | Density index | Threshold |
|---|---|---|---|---|---|---|---|
| CA_GGL_1, Globigerinita glutinata, 207 x 205 x 187 at 1.2016 um (LE) | 1.076e6 | 1.453e6 | 0.260 | 0.217 | 4.07e5 | 53346 | 32331 |
| GE_GTE_1, Globigerinoides tenellus, 227 x 220 x 241 at 1.2004 um (LE) | 2.011e6 | 2.593e6 | 0.224 | 0.278 | 6.06e5 | 48476 | 31634 |

Without the 8-voxel closing, hole filling alone left almost no enclosed void (test porosity 0.002 and 0.0008), because the apertures open the chambers to the exterior. Mid-slice QC figures are in `benchmarks/results/public/`.

### Failures and limits

* Figshare and some Imperial Box links were not downloadable from this host (HTTP 403 / 404); Digital Porous Media requires a login. The 2015 Imperial archive and the PANGAEA individual foram files were enough to run the tests.
* The first Estaillades 144^3 (image centre of the 400^3 crop) was unrepresentative (porosity 0.057); FFT formation factor failed there (non-positive effective conductivity). The porosity-matched percolating window succeeded.
* Ketton has no published k or F in the shipped header, so the solvers were not compared on that rock.
* Ambient-occlusion interior separation (the MFSE default of Titschack et al. 2018) is not implemented; morphoCT uses exterior erosion plus a convex-hull enclosure filter instead.
* A last chamber with a very wide primary aperture can remain connected to the umbilicus and be discarded by the enclosure-ratio filter (see CA_GGL_1 and CA_GUV_1 below).

## Foraminifera features: chambers, wall thickness and density

This section updates the Siccha et al. 2023 specimens after implementing chamber labelling, Hildebrand wall-thickness maps, and density calibration. The earlier table above used a morphological closing of radius 8 plus hole filling, which left the last chambers open and understated filled test volume. The figures below use `seal_method=hull_enclosure` (default): apertures narrower than 16 voxels (~19 um) are cut from the outside, then the void inside the convex hull is watershed-split and basins with a low shell-to-void contact ratio (the umbilicus and open bays) are discarded. Chamber labels are then re-ordered along a nearest-neighbour path that starts at the smallest chamber (a proloculus candidate). Methods were reimplemented from published descriptions; no ForamJ or MFSE source was copied.

Wall thickness is the Hildebrand and Ruegsegger (1997) diameter of the largest inscribed sphere, from PoreSpy `local_thickness` (which returns a radius in voxels; morphoCT reports `2 * radius * voxel_size` in um). Density uses the Iwasaki et al. (2015) relative CT number, CT = 1000 * (grey - air) / (calcite_standard - air), and an optional two-point linear fit. The Siccha 2023 scans do not include a calcite phantom, so all three runs are flagged UNCALIBRATED: `density_index` is the mean grey of an eroded shell, not g/cm^3, and shell mass is omitted. Background median grey is stored as an estimated `air_grey`. A synthetic three-chamber phantom with known count, wall diameter and density is in `tests/test_cli_foram.py`.

CLI: `morphoct foram examples/foram.yaml`, or flags (`--input`, `--shape` as Z,Y,X, `--voxel-size`, `--calcite-grey`, ...). Runner: `benchmarks/run_foram_siccha.py`. QC figures and `results.json` are in `benchmarks/results/foram_siccha/`.

Literature comparisons are to taxonomic diagnoses of the *final whorl*, not to a published total chamber count for these specimen IDs (the dataset paper does not report chamber counts). Pearson, Wade and Huber (2018), Atlas of Oligocene Planktonic Foraminifera, Cushman Foundation Special Publication 46, Chapter 16, p. 429-458, and the matching Mikrotax diagnoses: *Globigerinita glutinata* and *G. uvula* have 3 to 3.5 globular chambers in the final whorl; *G. uvula* is distinguished by a high to very high trochospire. Kennett and Srinivasan (1983), Neogene Planktonic Foraminifera, p. 80, and Mikrotax: *Globigerinoides tenellus* has four rapidly enlarging chambers in the final whorl.

| Specimen | Chambers (morphoCT) | Final-whorl (literature) | Shell (um^3) | Filled test (um^3) | Test porosity | Wall median / mean / SD (um) | Density |
|---|---|---|---|---|---|---|---|
| CA_GGL_1, G. glutinata, 207 x 205 x 187 at 1.2016 um | 5 | 3 to 3.5 | 1.076e6 | 2.752e6 | 0.609 | 7.21 / 6.35 / 1.65 | UNCALIBRATED, mean grey 53346 (air ~16013, threshold 32331) |
| GE_GTE_1, G. tenellus, 227 x 220 x 241 at 1.2004 um | 5 | 4 | 2.011e6 | 4.364e6 | 0.539 | 12.00 / 10.43 / 3.97 | UNCALIBRATED, mean grey 48476 (air ~17288, threshold 31634) |
| CA_GUV_1, G. uvula, 152 x 180 x 138 at 1.2016 um | 5 | 3 to 3.5 (high spire) | 4.325e5 | 1.602e6 | 0.730 | 4.81 / 4.71 / 1.35 | UNCALIBRATED, mean grey 53655 (air ~15811, threshold 31365) |

Per-chamber volumes (um^3), growth-path order starting at the smallest chamber:

* CA_GGL_1: 6.52e4, 7.01e4, 2.06e5, 4.04e5, 9.31e5
* GE_GTE_1: 1.05e5, 1.33e5, 2.03e5, 5.80e5, 1.33e6
* CA_GUV_1: 1.37e5, 3.34e5, 1.85e5, 1.97e5, 3.17e5

### Comparison and QC notes

* CA_GGL_1: five labelled chambers, of which the mid-slice QC shows three last-whorl chambers filled and one last chamber with a wide aperture left unlabelled (it stays connected to the umbilicus and is dropped by the enclosure filter). Two smaller earlier chambers are included. The count is a total, not a final-whorl count, and is a lower bound because of the missed last chamber. Wall thickness is 5 to 8 um on the outer wall, with locally thicker septa; the histogram is stepped at 2 * voxel_size (~2.4 um), as expected for discrete local thickness. Runtime 23.5 s.
* GE_GTE_1: mid-slice QC shows four final-whorl chambers plus one smaller earlier chamber, in agreement with the four-chamber final-whorl diagnosis. The large last chamber is recovered. Outer walls are thicker than in the two Globigerinita specimens (median 12.0 um), consistent with a macroperforate cancellate wall; local maxima reach ~22 um. Runtime 35.4 s.
* CA_GUV_1: five labelled chambers. Mid-slice QC shows three chambers in the lower spire; the last, hook-shaped chamber remains open to the exterior and is missed. The lower heart-shaped region may fuse two chambers. Because of the high trochospire, the true total should exceed the 3 to 3.5 final-whorl figure; 5 is a lower bound. The wall is thin (median 4.81 um), as expected for this microperforate species. Runtime 9.8 s.

### Unit tests

`pytest tests/test_cli_foram.py -q -k "not slow"`: 8 passed. These include a three-sphere phantom with known chamber count (3), a 4-voxel wall whose Hildebrand diameter is recovered to within 1 um, Iwasaki and linear calibrations that recover 2.71 g/cm^3 and the expected mass, an explicit UNCALIBRATED flag when no standard is given, aperture closing, hollow-cube volumes, and `morphoct foram` on the phantom (QC figures and `results.json`). Ambient-occlusion interior separation still raises NotImplementedError.

### Limits (feature work)

* Last chambers with a primary aperture wider than the enclosure criterion can be missed (GGL, GUV). Speijer et al. 2008 already noted that closing foramina can err by about 15 %.
* Growth-path order is a nearest-neighbour walk from the smallest chamber, not a verified ontogenetic sequence (kummerform or bullate last chambers can break volume order).
* Discrete local thickness is quantized at 2 * voxel_size; thin walls (a few voxels) have a coarse histogram.
* No calcite standard was scanned with these specimens, so g/cm^3 and shell mass are not reported here. The calibration path is covered by the phantom tests.
* Full ambient occlusion (Titschack et al. 2018 / MFSE) is still a stub.
