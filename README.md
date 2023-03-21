# Huoguo: Improved Photometry and Stellar Mass Distributions of Massive Galaxies in the Hyper Suprime-Cam Subaru Strategic Program

---- 

- Song Huang 

- We have recently repurposed this repo to carry out photometric analysis & stellar mass distribution estimates of 0.1<z<0.6 massive galaxies in the Hyper Suprime-Cam Subaru Strategic Program (HSC-SSP).
- We will use this repo to keep track of the progress of the analysis. We will store the `Jupyter` notebooks, `Python` scripts, and necessary `Markdown` notes in this repo.
## Introduction

- The public facing part of this project will use the [HSC PDR3 data release](https://hsc-release.mtk.nao.ac.jp/doc/index.php/data-access__pdr3/).

### General Workflow 

1. Select bright & extended objects from the HSC database as candidates of massive galaxies.
2. Estimate the photometric redshifts of the candidates & compile the spectroscopic redshifts available for these candidates. 
3. Estimate the mass-to-light ratio (M/L) of the candidate based on the HSC photometry & the best-available redshift estimates. 
4. Perform careful photometric analysis to estimate the stellar mass distribution of the candidate.
5. Summarize the analysis into 1-D stellar mass density profiles and the curve-of-growth (COG) of the candidate.

### Sample Selection 

- We will make an initial selection based on the magnitude & optical colors of the objects from the HSC database. We will use our previous sample from the `S16A` data release as a starting point. And we will also use the [`COSMOS2020` catalog](https://astroweaver.github.io/project/cosmos2020-galaxy-catalog/) to check the completeness of our sample.
- The sample selection will take the HSC full-depth full-color (FDFC) & bright star masks into account.

#### Galaxy Clusters 

- In addition to the general sample of massive galaxies, we will pay particular attention to the available catalogs of galaxy clusters as well. We will carry out the same photometric analysis of the massive central & satellite galaxies in these clusters regardless of whether they have made into our sample. 
- Currently, we are considering the cluster catalogs based on the [`redMaPPer`](https://github.com/erykoff/redmapper) & [`CAMIRA`](https://github.com/oguri/cluster_catalogs/tree/main/hsc_s20a_camira) algorithms. Both richness-based catalogs are available for the HSC footprint. We will also consider the following catalogs of galaxy clusters & groups: 
  - [Wen & Han 2021 for HSC PDR2](https://ui.adsabs.harvard.edu/abs/2021MNRAS.500.1003W/abstract). The catalogs are [available here](shttp://zmtt.bao.ac.cn/galaxy_clusters/catalogs.html)
  - [Zou et al. 2022 for DESI DR9 & HSC PDR3](https://ui.adsabs.harvard.edu/abs/2022RAA....22f5001Z/abstract). The catalogs are [available here](https://www.scidb.cn/en/detail?dataSetId=7797b553a23846a187b7746d8fc555a5)
  - [Yang et al. 2021 for DESI DR8](https://ui.adsabs.harvard.edu/abs/2021ApJ...909..143Y/abstract). There is also a version for [DESI DR9](https://gax.sjtu.edu.cn/data/DESI.html)

### Photometric Redshift  

- The HSC database has provided photometric redshifts estimated by the `Mizuki`, `DEmP`, and `DNNz` algorithms.
- We will also run [`Frankenz`](https://github.com/joshspeagle/frankenz) on the HSC PDR3 data release to estimate the photometric redshifts of the candidates.

### SED Fitting & M/L Estimation 

- We will use the [`Bagpipes`](https://github.com/ACCarnall/bagpipes) code as the main code to perform the SED fitting & M/L estimation.
- We will also use the [`Prospector`](https://github.com/bd-j/prospector) and [`dense-basis](https://github.com/kartheikiyer/dense_basis) codes to test the robustness of our estimates. 
- We will consider the impact of different choices of the Initial Mass Function (IMF), the dust attenuation law, and the stellar population synthesis models on the M/L estimates.

### Photometric Analysis

- The centerpiece of the photometric analysis is the extraction of the surface brightness profiles of the massive galaxies. For this purpose, we will build a non-parametric 2-D model of the galaxy or perform isophotal analysis of the galaxy.
  - We will use the `Python` code [`AutoProf-Legacy`](https://github.com/ConnorStoneAstro/AutoProf-Legacy) to perform the isophotal analysis. And we will also explore the application of the new [`AutoProf`](https://github.com/ConnorStoneAstro/AutoProf) code to build non-parametric models. 
- We will consider the impact of background subtraction, contaminating objects, blending scenarios, and the PSF on the surface brightness profiles.
- We will also consider the forward modeling approach to build the 2-D models of the galaxies.

#### Downloading Imaging Data 

- We will update the [`unagi`](https://github.com/dr-guangtou/unagi) code to download the imaging data from the HSC database.

----

## Naming 

- Huoguo (火锅; Hotpot) is a Chinese dish that is a combination of many ingredients.