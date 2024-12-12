# Makefile for corr2uvfits.

# on setonix:
# module load cfitsio/4.3.0 spack/default fftw/3.3.10
# spack install gsl pal
# spack load spack load gsl@2.7.1/th42khq pal

# on macos
# brew install llvm libomp icecube/icecube/pal openblas
# export PKG_CONFIG_PATH="$(brew --prefix)/opt/openblas/lib/pkgconfig"
# export LDFLAGS="-L$(brew --prefix)/opt/llvm/lib -L$(brew --prefix)/lib"
# export INCS="-I$(brew --prefix)/opt/llvm/include -I$(brew --prefix)/include"
# export CC=$(brew --prefix)/opt/llvm/bin/clang
# export CXX=$(brew --prefix)/opt/llvm/bin/clang++

# in the docker container
# export INCS="-I/usr/include/pal/"
# export LDFLAGS="-L/usr/lib"

CFLAGS ?= -g -O -D_FILE_OFFSET_BITS=64 -L.
CFITSIO_INCS=$(shell pkg-config --cflags cfitsio)
CFITSIO_LIBS=$(shell pkg-config --libs cfitsio)
PAL_LIBS=-lpal
BLAS_INCS=$(shell pkg-config --cflags openblas)
BLAS_LIBS=$(shell pkg-config --libs openblas)
LAPACK_INCS=$(shell pkg-config --cflags lapack)
LAPACK_LIBS=$(shell pkg-config --libs lapack)
FFTW3_INCS=$(shell pkg-config --cflags fftw3)
FFTW3_LIBS=$(shell pkg-config --libs fftw3)
GSL_INCS=$(shell pkg-config --cflags gsl)
GSL_LIBS=$(shell pkg-config --libs gsl)
COMMON_LIBS=-lm -fopenmp

TARGETS=gridvisdiff lssa_fg_simple lssa_fg_general prepare_diff ps_metrics combine_data

# Prefix for all installed files
PREFIX ?= /usr/local

all: $(TARGETS)

ps_metrics: ps_power.c uvfits.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS) $(LDFLAGS) \
		-o ps_metrics ps_power.c uvfits.c cspline.c \
		$(CFITSIO_LIBS) $(BLAS_LIBS) ${COMMON_LIBS} $(PAL_LIBS)

gridvisdiff: grid_vis_PB_chips.c uvfits.c primary_beamBNanalytic.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(GSL_INCS) $(BLAS_INCS) $(LDFLAGS) \
		-o gridvisdiff grid_vis_PB_chips.c uvfits.c primary_beamBNanalytic.c cspline.c \
		$(CFITSIO_LIBS) $(GSL_LIBS) $(BLAS_LIBS) $(COMMON_LIBS) -pg $(PAL_LIBS)

prepare_diff: prepare_cube_chips.c uvfits.c primary_beamDIFF.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(GSL_INCS) $(FFTW3_INCS) $(LDFLAGS) \
		-o prepare_diff prepare_cube_chips.c uvfits.c primary_beamDIFF.c cspline.c \
		$(CFITSIO_LIBS) $(GSL_LIBS) ${FFTW3_LIBS} $(COMMON_LIBS) $(PAL_LIBS)

combine_data: combine_data.c uvfits.c primary_beamDEV.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(BLAS_INCS) ${FFTW3_INCS} $(LDFLAGS) \
		-o combine_data combine_data.c uvfits.c primary_beamDIFF.c cspline.c \
		$(CFITSIO_LIBS) $(BLAS_LIBS) ${FFTW3_LIBS} ${LAPACK_LIBS} ${COMMON_LIBS} $(PAL_LIBS)

# TODO: it's confusing that the stripped compiles to simple, and simple compiles to general
lssa_fg_simple: fft_krig_stripped.c uvfits.c primary_beamDEV.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(GSL_INCS) $(BLAS_INCS) $(LDFLAGS) \
		-o lssa_fg_simple fft_krig_stripped.c uvfits.c primary_beamDEV.c cspline.c \
		$(CFITSIO_LIBS) $(GSL_LIBS) $(BLAS_LIBS) $(COMMON_LIBS) $(PAL_LIBS)

lssa_fg_general: fft_krig_simple.c uvfits.c primary_beamDEV.c cspline.c
	$(CC) $(CFLAGS) $(INCS) $(CFITSIO_INCS) $(GSL_INCS) $(BLAS_INCS) $(LDFLAGS) \
		-o lssa_fg_general fft_krig_simple.c uvfits.c primary_beamDEV.c cspline.c \
		$(CFITSIO_LIBS) $(GSL_LIBS) $(BLAS_LIBS) $(COMMON_LIBS) $(PAL_LIBS)

install: all
	mkdir -p $(PREFIX)/bin
	mv $(TARGETS) $(PREFIX)/bin

clean:
	rm -f *.o $(TARGETS)
