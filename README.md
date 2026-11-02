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


