# Makefile for corr2uvfits.

CFLAGS=-g -O -D_FILE_OFFSET_BITS=64 -L.
CFITSIO_INCS=$(shell pkg-config --silence-errors --cflags cfitsio)
CFITSIO_LIBS=$(shell pkg-config --silence-errors --libs cfitsio)
INCS=$(shell python -c "if len('${INCLUDE}')>0:print ' '.join(['-I ' + s for s in '${INCLUDE}'.split(':')])") -L${CFITSLIB} -I${CFITSINC}
#CHOLMOD_INCS=-I/usr/include/suitesparse/

CFITSIO_INCS=-I/usr/local/cfitsio-3.49/include/
CFITSIO_LIBS=-L/usr/local/cfitsio-3.49/lib/

BLAS_INCS=-I/usr/include/x86_64-linux-gnu/
BLAS_LIBS=-L/usr/lib/x86_64-linux-gnu/

SRCDIR=CONV2UVFITS

TARGETS=grid_bn fft_thermal prepare_diff
DRIPSTAR=dripsall prepare_lssa_drips lssa_drips_thermal
SKA=grid_ska prepare_ska lssa_ska

all: $(TARGETS)
ska: $(SKA)
drip: $(DRIPSTAR)

grid_vis: grid_vis.c uvfits.c primary_beam.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o gridvis grid_vis.c uvfits.c primary_beam.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp

grid_diff: grid_visDIFF.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	gcc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisdiff grid_visDIFF.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -pg
	
#grid_gauss: grid_visGAUSS.c uvfits.c primary_beamGAUSS.c cspline.c libsla.a
#	ranlib libsla.a
#	gcc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisgauss grid_visGAUSS.c uvfits.c primary_beamGAUSS.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -pg

grid_ska: chipska/grid_visDIFFBNUV.c chipska/uvfitska.c chipska/primary_beamBNanalytic.c cspline.c
#	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  $(BLAS_INCS)  -o gridvisska chipska/grid_visDIFFBNUV.c chipska/uvfitska.c chipska/primary_beamBNanalytic.c cspline.c $(CFITSIO_LIBS) $(BLAS_LIBS) -lstarlink_pal -lcfitsio  -lm -fopenmp -pg

#grid_sims: grid_visSIMS.c uvfits.c primary_beamSIMS.c cspline.c libsla.a
#	ranlib libsla.a
#	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(CHOLMOD_INCS) -o gridvissims grid_visSIMS.c uvfits.c primary_beamSIMS.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg


delay: obs_delayspectrum.c uvfits.c kde_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o obs_delayspectrum obs_delayspectrum.c uvfits.c kde_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg
	
extract: obs_extractvis.c uvfits.c kde_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o obs_extractvis obs_extractvis.c uvfits.c kde_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg


extractrts: obs_extractvis_rts.c uvfits.c kde_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o obs_extractvis_rts obs_extractvis_rts.c uvfits.c kde_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg


grid_sim: grid_visDIFFsimsUV.c uvfits.c primary_beamDIFFsimsUV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisdiffsims grid_visDIFFsimsUV.c uvfits.c primary_beamDIFFsimsUV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg

grid_fine: grid_visDIFFfine.c uvfits.c primary_beamDIFFfine.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisdifffine grid_visDIFFfine.c uvfits.c primary_beamDIFFfine.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg

grid_uv: grid_visDIFFfineUV.c uvfits.c primary_beamDIFFfineUV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisdifffineuv grid_visDIFFfineUV.c uvfits.c primary_beamDIFFfineUV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg

grid_bn: grid_visDIFFBNUV.c uvfits.c primary_beamBNanalytic.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS)  -o gridvisdifffine grid_visDIFFBNUV.c uvfits.c primary_beamBNanalytic.c cspline.c $(CFITSIO_LIBS) $(BLAS_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg
	
grid_all: grid_vis_PB_chips.c uvfits.c primary_beamBNanalytic.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS)  -o gridvisdifffine grid_vis_PB_chips.c uvfits.c primary_beamBNanalytic.c cspline.c $(CFITSIO_LIBS) $(BLAS_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg
	
kde_all: grid_kde_memory.c uvfits.c primary_beamBNanalytic.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS)  -o grid_kde_memory grid_kde_memory.c uvfits.c primary_beamBNanalytic.c cspline.c $(CFITSIO_LIBS) $(BLAS_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg

grid_stack: grid_visDIFFBNUVstack.c uvfits.c primary_beamBNanalyticWstack.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisdifffine grid_visDIFFBNUVstack.c uvfits.c primary_beamBNanalyticWstack.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg


grid_bnsim: grid_visDIFFBNsim.c uvfits.c primary_beamBNsimstack.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisdifffine grid_visDIFFBNsim.c uvfits.c primary_beamBNsimstack.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg


grid_gauss: grid_visGAUSSanalytic.c uvfits.c primary_beamGAUSSanalytic.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisdifffine grid_visGAUSSanalytic.c uvfits.c primary_beamGAUSSanalytic.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg

grid_ultra: grid_vis_ultra.c uvfits.c primary_beam_ultra.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisultra grid_vis_ultra.c uvfits.c primary_beam_ultra.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg

grid_check: grid_visCHECK.c uvfits.c primary_beamCHECK.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvischeck grid_visCHECK.c uvfits.c primary_beamCHECK.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg


lyman: lyman.c uvfits.c primary_beam_lyman.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lyman lyman.c uvfits.c primary_beam_lyman.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg

prepare_lyman: prepare_lyman.c uvfits.c primary_beam_lyman.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_lyman prepare_lyman.c uvfits.c primary_beam_lyman.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_diff: prepare_diff.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_diff prepare_diff.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_cube: prepare_cube_chips.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_diff prepare_cube_chips.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_full: prepare_diff_full.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_diff prepare_diff_full.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_add: prepare_add.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_add prepare_add.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_ska: chipska/prepare_diff.c chipska/uvfitska.c chipska/primary_beamBNanalytic.c cspline.c
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_ska chipska/prepare_diff.c chipska/uvfitska.c chipska/primary_beamBNanalytic.c cspline.c $(CFITSIO_LIBS) -lstarlink_pal -lfftw3 -lcfitsio   -lm -fopenmp

prepare_ska_woden: chipska/prepare_diff_woden.c uvfits.c chipska/primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_ska_woden chipska/prepare_diff_woden.c uvfits.c chipska/primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_ultra: prepare_ultra.c uvfits.c primary_beam_ultra.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_ultra prepare_ultra.c uvfits.c primary_beam_ultra.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

#kde:	kdemain.c uvfits.c kde_build3D.c cspline.c libsla.a
#	ranlib libsla.a
#	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o kdemain kdemain.c uvfits.c kde_build3D.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

kde:	kdemain.c uvfits.c kde_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o kdemain kdemain.c uvfits.c kde_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg


kdec:	kdemain_cube.c uvfits.c kde_build3D_freq_cube.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o kdemain kdemain_cube.c uvfits.c kde_build3D_freq_cube.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

vis: output_vis.c uvfits.c kde_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o output_vis output_vis.c uvfits.c kde_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg

checksum: checksum.c uvfits.c kde_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o check_metrics checksum.c uvfits.c kde_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg

ps_metric: ps_power_metrics.c uvfits.c kde_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o ps_metrics ps_power_metrics.c uvfits.c kde_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg
	
ps_power: ps_power.c uvfits.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o ps_metrics ps_power.c uvfits.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg

mmode: output_vis_mmode.c uvfits.c kde_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o output_mmode output_vis_mmode.c uvfits.c kde_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg

redundant: obs_redundant.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o obs_redundant obs_redundant.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg

bispec3d: obs_bispectrum.c uvfits.c bispectrum_function.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o obs_bispectrum obs_bispectrum.c uvfits.c bispectrum_function.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg

bidelay: obs_delaybispectrum.c uvfits.c delaybispectrum_function.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o obs_delaybispectrum obs_delaybispectrum.c uvfits.c delaybispectrum_function.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -pg

kdecs:	kdemain_cube_stokes.c uvfits.c kde_build3D_freq_cube_stokes.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o kdemain kdemain_cube_stokes.c uvfits.c kde_build3D_freq_cube_stokes.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

wpmain:	winpowmain.c uvfits.c winpow_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o winpowmain winpowmain.c uvfits.c winpow_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

clos:	closuremain.c uvfits.c closure_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o closure closuremain.c uvfits.c closure_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

closf:	closuremain.c uvfits.c closure_function.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o closure closuremain.c uvfits.c closure_function.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

bispectrum:	bispectrum_main.c uvfits.c closure_function.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o bispectrum bispectrum_main.c uvfits.c closure_function.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg


winpow: compute_delayspectrum_windowpower.c uvfits.c winpow_build3D_freq.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o winpow compute_delayspectrum_windowpower.c uvfits.c winpow_build3D_freq.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp


kdeext:	kdemainext.c uvfits.c kde_extract.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o kdemain kdemainext.c uvfits.c kde_extract.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

dripsall:	../drips/grid_visDIFF_drips.c uvfits.c ../drips/primary_beamDIFF_drips.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvisdrips ../drips/grid_visDIFF_drips.c uvfits.c ../drips/primary_beamDIFF_drips.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

dev: grid_visDEV.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvis grid_visDEV.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

noise: compute_noise_power.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o compute_noise_power compute_noise_power.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

noisefrb: compute_noise_power_frb.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o compute_noise_power compute_noise_power_frb.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

noiseqa: compute_noise_power_qa.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o compute_noise_power compute_noise_power_qa.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

sim: grid_vis_sim.c uvfits.c calc_sims.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvis_sim grid_vis_sim.c uvfits.c calc_sims.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -pg

#dev: grid_visDEV.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c libsla.a
#	ranlib libsla.a
#	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvis grid_visDEV.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

lssa: lssa_fg.full.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa lssa_fg.full.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

makecov: make_ptsrc_covariance.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o make_ptsrc_covariance make_ptsrc_covariance.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

lssa_fg: lssa_fg.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg lssa_fg.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio  -lm -fopenmp -llapack -lblas

lssa_fg_test: lssa_fg.test.iter.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg lssa_fg.test.iter.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

#lssa_fg_thermal: lssa_fg.krig.c uvfits.c primary_beamDEV.c cspline.c libsla.a
#	ranlib libsla.a
#	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg lssa_fg.krig.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

lssa_fg_thermal: lssa_fg.krig.billabong.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg lssa_fg.krig.billabong.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

fft_thermal: fft_krig.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg fft_krig.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -lgsl -lgslcblas
	
fft_mdss: fft_mdss.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_mdss fft_mdss.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -lgsl -lgslcblas
	
fft_nfft: fft_nfft.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_nfft fft_nfft.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -lgsl -lgslcblas
	

fft_simple: fft_krig_stripped.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_simple fft_krig_stripped.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -lgsl -lgslcblas
	
	
fft_test: fft_test.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_test fft_test.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -lgsl -lgslcblas
	
fft_ska: chipska/fft_krig.c chipska/uvfitska.c chipska/primary_beamBNanalytic.c cspline.c
#	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_ska chipska/fft_krig.c chipska/uvfitska.c chipska/primary_beamBNanalytic.c cspline.c $(CFITSIO_LIBS) -lstarlink_pal -lcfitsio -lm -fopenmp -lgsl -lgslcblas

fft_thermalsub: fft_krig_subtract.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg fft_krig_subtract.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -lgsl -lgslcblas

fft_mean: fft_krig_submean.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_sub fft_krig_submean.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

skew: skew_spectrum.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o skew_spec skew_spectrum.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

dft_kde: dft_kde.krig.c uvfits.c kde_build.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o dft_kde dft_kde.krig.c uvfits.c kde_build.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp

lssa_fg_gpr: lssa_fg.krig_gpr.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg lssa_fg.krig_gpr.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas


lssa_ska: chipska/lssa_fg.SKA.c uvfits.c chipska/primary_beamSKA.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_ska chipska/lssa_fg.SKA.c uvfits.c chipska/primary_beamSKA.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

lssa_fg_bp: lssa_fg.thermal_bpweights.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg lssa_fg.thermal_bpweights.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

lssa_drips: ../drips/lssa_fg.test.iter_drips.c uvfits.c ../drips/primary_beamDEV_drips.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_drips ../drips/lssa_fg.test.iter_drips.c uvfits.c ../drips/primary_beamDEV_drips.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

lssa_drips_thermal: ../drips/lssa_fg.drips_thermal.c uvfits.c ../drips/primary_beamDIFF_drips.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_drips ../drips/lssa_fg.drips_thermal.c uvfits.c ../drips/primary_beamDIFF_drips.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

lssa_fg_invert: lssa_fg.invert.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_invert lssa_fg.invert.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

lssa_fg_eigen: lssa_fg.test.eigen.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_eigen lssa_fg.test.eigen.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -llapack -lblas

make_weights: make_weights.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o make_weights make_weights.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -pg

sum_w: sum_w.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o sum_w sum_w.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -pg

make_weights_choose: make_weights_choose.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o make_weights make_weights_choose.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -pg

make_w: combine_w.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o make_w combine_w.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio -lm -fopenmp -pg

gen_beams: gen_beams.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gen_beams gen_beams.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio -lm -fopenmp

gen_w: gen_w.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gen_w gen_w.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio -lm -fopenmp


combine: cholmod_combine.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o combine cholmod_combine.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp

w_convolve: cholmod_w_convolve.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o w_convolve cholmod_w_convolve.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp

combine_bb: cholmod_combine_bb.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o combine_bb cholmod_combine_bb.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp


combine_w: cholmod_combine_w.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o combine_w cholmod_combine_w.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

combine_w_divide: cholmod_combine_w_divide.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o combine_w cholmod_combine_w_divide.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

combine_w_choose: cholmod_combine_w_choose.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o combine_w cholmod_combine_w_choose.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

combine_wtest: cholmod_combine_w_test.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o combine_w cholmod_combine_w_test.c uvfits.c primary_beamDEV.c circ_shift.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_lssa: prepare_lssa.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_lssa prepare_lssa.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_kde: prepare_kde_freq3D_grid.c uvfits.c kde_build.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_kde prepare_kde_freq3D_grid.c uvfits.c kde_build.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_kden: prepare_kde_freq3D_cube.c uvfits.c kde_build.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_kde prepare_kde_freq3D_cube.c uvfits.c kde_build.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp

prepare_kdecg: prepare_kde_freq3D_cubegrid.c uvfits.c kde_build.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_kde prepare_kde_freq3D_cubegrid.c uvfits.c kde_build.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp

prepare_ds: prepare_delayspectrum.c uvfits.c kde_build.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_ds prepare_delayspectrum.c uvfits.c kde_build.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_lssa_drips: ../drips/prepare_diff_drips.c uvfits.c ../drips/primary_beamDIFF_drips.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_lssa_drips ../drips/prepare_diff_drips.c uvfits.c ../drips/primary_beamDIFF_drips.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_clean_lssa: prepare_clean_lssa.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_clean_lssa prepare_clean_lssa.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

prepare_poly: prep_poly.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prep_poly prep_poly.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp -lgsl -lgslcblas

combine_data: combine_data.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o combine_data combine_data.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

combine_stack: combine_data_stack.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o combine_stack combine_data_stack.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp

make_beams_grid: make_beams_grid.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o make_beams_grid make_beams_grid.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lsla -lcfitsio   -lm -fopenmp


dev_sparse: grid_visDEV.c uvfits.c primary_beam.dev_sparse.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o gridvis grid_visDEV.c uvfits.c primary_beam.dev_sparse.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp -pg

cholmod_construct_projection: cholmod_construct_projection.c primary_beamDEV.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o construct cholmod_construct_projection.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp

cholmod_build_projection: cholmod_build_projection.c primary_beamDEV.c cholmod_extras.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o build cholmod_build_projection.c uvfits.c primary_beamDEV.c cholmod_extras.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp

get_diag: get_diag.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o get_diag get_diag.c uvfits.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm -fopenmp

calc_grid: calc_grid.c uvfits.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o calc_grid calc_grid.c uvfits.c cspline.c $(CFITSIO_LIBS) -lsla -lcfitsio   -lm 

test_readuvfits: test_readuvfits.c uvfits.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) -o $@ $^ $(CFITSIO_LIBS) -lcfitsio -lsla -lm

libsla.a:
	cd SLALIB_C ; make
	rm -f SLALIB_C/*.o
clean:
	rm -f *.o $(TARGETS) libsla.a SLALIB_C/*.o 
	rm -rf corr2uvfits.dSYM

cleanfiles:
	rm /data/gridded_vis/*xx*.cst
	rm /data/gridded_vis/*xx*.dat
	rm /data/gridded_vis/*yy*.cst
	rm /data/gridded_vis/*yy*.dat

