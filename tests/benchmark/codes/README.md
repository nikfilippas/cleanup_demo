# Guide to benchmark codes

The source codes included in this directory can be used to reproduce all benchmark files used for the CCL unit tests. Several of these tests are based on the same software (e.g. the CosmoMAD library), but the output of that software was compared against other codes (i.e. CCL has not just been compared to CosmoMAD or any one particular package, but to several of them).

* **distances_bm.py** : this script generates the comoving distances used in tests distances_cosmomad_lowz:model_1-5. To run it (and several of the other scripts below) you need to download and install the cosmomad library and python wrapper (https://github.com/damonge/CosmoMAD). The version of CosmoMAD used to produce the benchmarks is the release under tag `ccl_validation`. CosmoMAD is publicly available under the GSL3.0 license.
* **distances_hiz_bm.py** : generates the high-redshift distances used in tests distances_cosmomad_hiz:model_1-3
* **growth_lowz_bm.py** and **growth_hiz_bm.py** : same as distances_bm and distances_hiz_bm for the growth functions used in tests growth_lowz:model_1-5 and growth_hiz:model_1-3
* **create_CLASS_distance_benchmarks.ipynb** : this notebook generates the distance benchmarks using `CLASS` for the tests distances_class:model1-11, where models 7-11 include massive neutrinos. The benchmarks are created for the low redshifts (`lowz`), high redshifts (`hiz`), and the whole redshift range (`allz`, using 10 logarithmically spaced redshifts between 0.01 and 1000)
* **multiple_neutrino_distances.ipynb** : this notebook uses `astropy` to generate the distance benchmarks for the tests distances_astropy_mnu_lowz:model1-5 and distances_astropy_mnu_hiz:model1-5
* **growth_allz.py** : generates the growth factors for the `allz` redshift range. Used in the tests growth_allz:model1-5
* **bbks_bm.py** and **ehpk_bm.py** : generates the benchmark BBKS and Eisenstein&Hu power spectra used in bbks:model_1-3 and eh:model_1.
* **sigmaM_bm.py** : generates the benchmark sigma(M) values used in sigmam:model_1-3.
* **mfunc_bm.py** : generates the mass function predictions used in massfunc:model_1.
* **cl_corr_bm** : this folder contains the software needed to reproduce the benchmark files for angular power spectra and correlation functions used in tests cls:histo, cls:analytic, corrs:analytic_bessel and corrs:analytic_fftlog. The C code within that folder (limberjack) should first be compiled before running the script run_all.py.
* **cl_cmbl_bm.py** : generates the benchmark file for the CMB lensing angular power spectrum used in test cls:cmblens.
* **bcm_bm.c** : this file shows an excerpt from a modified version of CLASS that was used to generate a benchmark power spectrum containing baryon corrections (used for test bcm:model_1). The file is just an excerpt, and therefore won't compile or run on its own.
* **3dcorr_benchmark.ipynb** : this python notebook can be used to generate the files used in tests corrs_3d:model_1-3.

Notice that for the cosmic emulator checks, we benchmark against simulated power spectra directly provided by the authors of the cosmic emulator work (Lawrence et al., 2017, arXiv:1705,03388). In addition, VARRIC tests discussed in the CCL paper assess the accuracy of CCL splines against direct CLASS outputs. Similarly for Figure 3 of the CCL paper. These checks are currently not automated. Finally, predictions of angular power spectra delivered by Angpow are benchmarked against a slower, brute force integration method that is available within CCL.