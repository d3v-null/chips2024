# Makefile for corr2uvfits.

# on setonix:
# module load cfitsio/4.3.0 spack/default fftw/3.3.10
# spack install gsl pal
# spack load spack load gsl@2.7.1/th42khq pal

# on macos
# brew install llvm libomp icecube/icecube/pal openblas
# export PKG_CONFIG_PATH="$(brew --prefix)/opt/openblas/lib/pkgconfig"
# export LDFLAGS="-L$(brew --prefix)/opt/llvm/lib -L$(brew --prefix)/lib"
# export CFLAGS="-g -O -D_FILE_OFFSET_BITS=64 -L. -I$(brew --prefix)/opt/llvm/include -I$(brew --prefix)/include"
# export CC=$(brew --prefix)/opt/llvm/bin/clang
# export CXX=$(brew --prefix)/opt/llvm/bin/clang++

CFLAGS ?= -g -O -D_FILE_OFFSET_BITS=64 -L.
BLAS_INCS=$(shell pkg-config --silence-errors --cflags openblas)
BLAS_LIBS=$(shell pkg-config --silence-errors --libs openblas)
GSL_INCS=$(shell pkg-config --silence-errors --cflags gsl)
GSL_LIBS=$(shell pkg-config --silence-errors --libs gsl)
CFITSIO_INCS=$(shell pkg-config --silence-errors --cflags cfitsio)
CFITSIO_LIBS=$(shell pkg-config --silence-errors --libs cfitsio)

TARGETS=gridvisdiff lssa_fg_simple prepare_diff

# Prefix for all installed files
PREFIX ?= /usr/local

all: $(TARGETS)

gridvisdiff: grid_vis_PB_chips.c uvfits.c primary_beamBNanalytic.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS) $(GSL_INCS) $(LDFLAGS) -o gridvisdiff grid_vis_PB_chips.c uvfits.c primary_beamBNanalytic.c cspline.c $(CFITSIO_LIBS) $(BLAS_LIBS)  -lm -fopenmp -pg -lpal

prepare_diff: prepare_cube_chips.c uvfits.c primary_beamDIFF.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(GSL_INCS) $(LDFLAGS) -o prepare_diff prepare_cube_chips.c uvfits.c primary_beamDIFF.c cspline.c $(CFITSIO_LIBS) -lfftw3   -lm -fopenmp -lpal

fft_simple: fft_krig_stripped.c uvfits.c primary_beamDEV.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS) $(GSL_INCS) $(LDFLAGS) -o lssa_fg_simple fft_krig_stripped.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) $(GSL_LIBS) -lm -fopenmp -lgsl -lgslcblas -lpal

fft_general: fft_krig_simple.c uvfits.c primary_beamDEV.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS) $(GSL_INCS) $(LDFLAGS) -o lssa_fg_general fft_krig_simple.c uvfits.c primary_beamDEV.c cspline.c $(CFITSIO_LIBS) $(GSL_LIBS) -lm -fopenmp -lgsl -lgslcblas -lpal

install: all
	mkdir -p $(PREFIX)/bin
	mv $(TARGETS) $(PREFIX)/bin

clean:
	rm -f *.o $(TARGETS)
