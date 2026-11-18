CHIPS is a simple visibility-based power spectrum estimator for 21cm data. It is designed to work with data from the MWA, but has been adapted for SKA-Low data. It uses time-interleaving to produce power spectra that are not noise power-biased, and visibility weights to estimate the thermal noise uncertainty.

CHIPS grids channel-based visibilities onto a uv-plane using a Blackman-Harris gridding kernel, matched to the Fourier primary beam size at the centre of the band. w-terms are not accounted for, and instead a cut is made for w > 30 wavelengths.

Spectral analysis uses a simple weighted Fourier Transform, and a Blackman-Nuttall window. Earlier versions of CHIPS used LSSA with an inverse covariance weigthing comprising thermal noise and a foreground model, and (optionally) kriging (to interpolate missing channels). This current version of CHIPS uses inverse variance weighting only using the number of samples in the dataset, and a regular FFT - no foreground covariance or kriging are used.

CHIPS software:

    ANSI C, compiled with gcc, using OpenMP parallelisation
    dependencies required: CFITSIO, LAPACK, BLAS, OpenMP, LIBSLA (astronomy)
    file input: UVFITS
    file outputs: complex floats, doubles (binary files, little endian)

CHIPS is composed of three separate codes, all written in ANSI-C:
1. grid_vis_PB_chips.c --> gridvisdiff: Reads UVFITS calibrated data files and grids onto the uv-plane (u,v,nu). This code can be run over multiple UVFITS files to grid onto the same (u,v,nu) plane.
2. prepare_cube_chips.c --> prepare_diff: Reads gridded visibility files + noise files + weights files, folds onto a half-uv plane and rearranges the data structure.
3. fft_krig_stripped.c --> lssa_fg_simple: performs spectral estimation and cylindrical averaging of the cross power spectra.

Usage of each program can be obtained by running the binary without any command line arguments. For example:

Usage: ./gridvisdifffine <options> uvfits_filename obs_id output_extension band 
	 -p period (seconds)
	 -ew flag 14m EW baselines
	 -c chanwidth (Hz)
	 -n bottom frequency (Hz)
	 -field fieldnum (0=EoR0, 1=EoR1, >2 anything)
	 -u umax

Usage: <options> ./prepare_diff input_ext Nchan freq_index_start pol output_extension band
	 -p period (seconds)
	 -c chanwidth
	 -u umax
	 -n bottom frequency (Hz)

Usage: <options> ./lssa_fg_simple input_extension Nchan nbins pol maxu output_extension bias_mode (0/10/11/12/13) band 
	 -p period
	 -c chanwidth
	 -f flag_dc

Here: band = 1 (EoR High-band); band = 0 (EoR low-band)
bias_mode sets the sub-band of the power spectrum estimation:
0 = full band
10 = top half of band
11 = bottom half of band
12 = middle half of band

The code also uses environment variables to define the location of the intermediate files and output files:
CODEDIR=/location of binaries
OUTPUTDIR=/location of intermediate files
DATADIR=/location of uvfits data
INPUTDIR=/location of intermediate files

export CODEDIR
export DATADIR
export OUTPUTDIR
export INPUTDIR


Typical usage for an EoR high-band observation from MWA might be:

./gridvisdifffine file.uvfits 145000000 testdata 1
./prepare_diff testdata 384 0 'yy' testdata 1
./lssa_fg_simple testdata 384 80 'yy' 300. testdata 0 1



## Measuring the signal normalisation (`gridvisdiff -S`)

Two visibilities see the same sky mode only if they are close in the uv-plane: for the MWA their sky signal is 50% correlated at a separation of about 1.75 wavelengths. The gridding kernel is 10 wavelengths wide, so a cell averages visibilities whose sky signal is only partly the same, and part of the sky power is lost. The more densely the uv-plane is filled (more baselines, more observations at different hour angles), the larger the loss. A fixed normalisation constant cannot follow this.

`gridvisdiff -S <seed>` measures it. It grids the observation exactly as usual (same positions, weights, flags and kernel) but replaces the value of every visibility by a simulated random sky:

- the sky is a field of unit-variance complex Gaussian numbers on a uv grid of 0.5 wavelengths, Hermitian, independent between channels;
- a visibility is that field seen through a Gaussian aperture footprint (`-W`, sigma in metres, default 1.7, which gives the 1.75 wavelengths above for the MWA tile at 185 MHz);
- the sky depends only on the seed, the channel frequency and (u,v), so every observation gridded with the same seed sees the same sky. w and the pointing of the primary beam are ignored.

The outputs have the usual names and go through `prepare_diff`, `combine_data` and `lssa` unchanged, so use a separate extension. `sky_reference.py` then writes the grid that would result with no loss (the same variance in every cell that has data, same weights). The ratio of the two power spectra is the correction for the data that were gridded and combined the same way:

    ./gridvisdiff file.uvfits 145000000 sky 1 -S 1
    ./prepare_diff sky 384 0 'yy' sky 1
    python3 sky_reference.py yy sky skyref 384
    ./lssa_fg_simple sky 384 80 'yy' 300. sky 0 1
    ./lssa_fg_simple skyref 384 80 'yy' 300. skyref 0 1
    # correction(k) = crosspower(skyref) / crosspower(sky)

Without `-S` the outputs of `gridvisdiff` are unchanged. `make test_sky` checks the simulated sky itself (variance, correlation against separation, Hermitian symmetry); `tests/` has a small synthetic uvfits generator and a check of the gridded simulated sky against direct evaluation.
