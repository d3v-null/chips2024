# syntax=docker/dockerfile:1
# cross-platform, cpu-only dockerfile for CHIPS
# on amd64, arm64
# ref: https://docs.docker.com/build/building/multi-platform/
FROM python:3.11-slim-bookworm

# Suppress perl locale errors
ENV LC_ALL=C
RUN apt-get update -y \
    && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    tzdata \
    build-essential \
    pkg-config \
    cmake \
    curl \
    git \
    autoconf \
    libtool \
    unzip \
    wget \
    zip \
    gfortran \
    zlib1g-dev \
    liberfa-dev \
    libstarlink-pal-dev \
    libgsl-dev \
    libopenblas-dev \
    libfftw3-dev \
    && apt-get clean all \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    && apt-get -y autoremove

# install pal
# WORKDIR /pal
# RUN wget https://github.com/Starlink/pal/releases/download/v0.9.8/pal-0.9.8.tar.gz \
#     && tar -xvf pal-0.9.8.tar.gz \
#     && cd pal-0.9.8 \
#     && ./configure --prefix=/usr/local --without-starlink \
#     && make \
#     && make install \
#     && ldconfig \
#     && rm -rf /pal

# install cfitsio
WORKDIR /cfitsio
RUN wget http://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.3.0.tar.gz \
    && tar -zxvf cfitsio-4.3.0.tar.gz \
    && cd cfitsio-4.3.0/ \
    && CFLAGS="-O3" ./configure --prefix=/usr/local --enable-reentrant --enable-ssse3 --enable-sse2 --disable-curl \
    && make -j $(nproc) \
    && make install \
    && ldconfig \
    && rm -rf /cfitsio

# provide source code files from the host
ADD . /chips
WORKDIR /chips
RUN make install PAL_LIBS="-lstarlink_pal" PREFIX=/usr/local

# docker build . -t d3vnull0/chips2024:latest
# docker run -it --rm d3vnull0/chips2024:latest /bin/bash
# docker push d3vnull0/chips2024:latest