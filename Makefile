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

TARGETS=grid_all fft_simple prepare_cube


all: $(TARGETS)



grid_all: grid_vis_PB_chips.c uvfits.c primary_beamBNanalytic.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS)  -o gridvisdifffine grid_vis_PB_chips.c uvfits.c primary_beamBNanalytic.c cspline.c $(CFITSIO_LIBS) $(BLAS_LIBS) -lcfitsio  -lm -fopenmp -pg
	
prepare_cube: prepare_cube_chips.c uvfits.c primary_beamDIFF.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o prepare_diff prepare_cube_chips.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3 -lcfitsio   -lm -fopenmp


fft_simple: fft_krig_stripped.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_simple fft_krig_stripped.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lcfitsio -lm -fopenmp -lgsl -lgslcblas
	
fft_general: fft_krig_simple.c uvfits.c primary_beamDEV.c cspline.c libsla.a
	ranlib libsla.a
	cc $(CFLAGS) $(INCS) $(CFITSIO_INCS)  -o lssa_fg_general fft_krig_simple.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) -lcfitsio -lm -fopenmp -lgsl -lgslcblas
	

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

