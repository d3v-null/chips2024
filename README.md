CHIPS is a simple visibility-based power spectrum estimator for 21cm data. It is designed to work with data from the MWA, but has been adapted for SKA-Low data. It uses time-interleaving to produce power spectra that are not noise power-biased, and visibility weights to estimate the thermal noise uncertainty.

CHIPS is composed of three separate codes, all written in ANSI-C:
1. grid_vis_PB_chips.c --> gridvisdiff
