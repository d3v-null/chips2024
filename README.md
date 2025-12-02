CHIPS is a simple visibility-based power spectrum estimator for 21cm data. It is designed to work with data from the MWA, but has been adapted for SKA-Low data. It uses time-interleaving to produce power spectra that are not noise power-biased, and visibility weights to estimate the thermal noise uncertainty.

CHIPS is composed of three separate codes, all written in ANSI-C:
1. grid_vis_PB_chips.c --> gridvisdiff: Reads UVFITS calibrated data files and grids onto the uv-plane (u,v,nu). This code can be run over multiple UVFITS files to grid onto the same (u,v,nu) plane.
2. prepare_cube_chips.c --> prepare_diff: Reads gridded visibility files + noise files + weights files, folds onto a half-uv plane and rearranges the data structure.
3. fft_stripped.c --> lssa_fg_simple: performs spectral estimation and cylindrical averaging of the cross power spectra.


