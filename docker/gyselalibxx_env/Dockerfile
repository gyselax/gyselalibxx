# SPDX-License-Identifier: MIT

FROM ubuntu:noble

COPY bash_run /bin/
ENV BASH_ENV=/etc/profile
ENV PATH=/opt/cmake-3.27/bin${PATH:+:${PATH}}
SHELL ["/bin/bash", "-c"]

# Layers are ordered from least to most frequently updated so that a version bump only
# rebuilds the layers that follow it. Each layer removes its own temporary files as files
# deleted in a later layer still take up space in the image.
# Libraries built for a given Kokkos backend are installed in /opt/<backend>/<Library>.
# Dependencies are found by passing CMAKE_PREFIX_PATH=/opt/<backend>.

# System packages
RUN chmod +x /bin/bash_run \
 && export DEBIAN_FRONTEND=noninteractive \
 && echo tzdata tzdata/Areas string Etc | debconf-set-selections \
 && echo tzdata tzdata/Zones/Etc string UTC | debconf-set-selections \
 && apt-get update -y \
 && apt-get upgrade -y \
 && apt-get install -y --no-install-recommends \
    apt-transport-https \
    apt-utils \
    ca-certificates \
    wget \
 && echo "deb [ arch=amd64 ] https://repo.pdi.dev/ubuntu noble main" > /etc/apt/sources.list.d/pdi.list \
 && wget -O /etc/apt/trusted.gpg.d/pdidev-archive-keyring.gpg https://repo.pdi.dev/ubuntu/pdidev-archive-keyring.gpg \
 && chmod a+r /etc/apt/trusted.gpg.d/pdidev-archive-keyring.gpg \
 && apt-get update -y \
 && apt-get install -y --no-install-recommends \
    build-essential \
    clang-18 \
    libomp-dev \
    clang-format-14 \
    doxygen \
    pkg-config \
    git \
    libfftw3-dev \
    liblapacke-dev \
    libopenmpi-dev \
    openssh-client \
    pdidev-archive-keyring \
    hdf5-tools \
    pdiplugin-mpi-openmpi \
    pdiplugin-decl-hdf5-openmpi \
    pdiplugin-set-value \
    python3-pip \
    python3-venv \
    libpdi-dev \
    cmake \
 && ln -sf /usr/bin/clang-18 /usr/bin/clang \
 && ln -sf /usr/bin/clang++-18 /usr/bin/clang++ \
 && apt-get purge -y \
    apt-transport-https \
    apt-utils \
    wget \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && rm -rf /var/lib/apt/lists/*

# CI user
RUN useradd -d /ci -m -U ci \
 && mkdir /data \
 && chown -R ci:ci /ci /data

# Python environment
RUN python3 -m venv /opt/venv \
 && /opt/venv/bin/python -m pip install -U --no-cache-dir \
    PyYAML \
    dask \
    h5py \
    imageio \
    matplotlib \
    numpy \
    pip \
    scipy \
    xarray \
    'pylint==3.3.9' \
 && chown -R ci:ci /opt/venv

# Googletest
RUN git clone --branch v1.17.0 --depth 1 https://github.com/google/googletest.git \
 && cmake -DCMAKE_BUILD_TYPE=Release -B build -S googletest \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/googletest \
 && rm -rf build googletest

# Cppcheck
RUN git clone --branch 2.22.0 --depth 1 https://github.com/cppcheck-opensource/cppcheck.git \
 && cd cppcheck \
 && cmake -S . -B build -DUSE_MATCHCOMPILER=ON -DCMAKE_BUILD_TYPE=Release \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build \
 && cd .. \
 && rm -rf cppcheck

# Ginkgo
RUN git clone --branch v1.11.0 --depth 1 https://github.com/ginkgo-project/ginkgo.git \
 && cd ginkgo \
 && cmake -DCMAKE_BUILD_TYPE=Release -DGINKGO_BUILD_TESTS=OFF -DGINKGO_BUILD_EXAMPLES=OFF -DGINKGO_BUILD_BENCHMARKS=OFF -B build \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build \
 && cd .. \
 && rm -rf ginkgo

# Kokkos-family dependencies (serial backend)
RUN git clone --branch 5.1.1 --depth 1 https://github.com/kokkos/kokkos.git \
 && git clone --branch v1.1.0 --depth 1 https://github.com/kokkos/kokkos-fft.git \
 && git clone --branch 5.1.1 --depth 1 https://github.com/kokkos/kokkos-kernels.git \
 && git clone --branch develop --depth 1 https://github.com/kokkos/kokkos-tools.git \
 && cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_DEPRECATED_CODE_4=OFF -DKokkos_ENABLE_SERIAL=ON -B build -S kokkos \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/serial/Kokkos \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/serial -B build -S kokkos-fft \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/serial/KokkosFFT \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/serial -B build -S kokkos-kernels \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/serial/KokkosKernels \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/serial -B build -S kokkos-tools \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/serial/KokkosTools \
 && rm -rf build \
 && rm -rf kokkos kokkos-fft kokkos-kernels kokkos-tools

# Kokkos-family dependencies (OpenMP backend)
RUN git clone --branch 5.1.1 --depth 1 https://github.com/kokkos/kokkos.git \
 && git clone --branch v1.1.0 --depth 1 https://github.com/kokkos/kokkos-fft.git \
 && git clone --branch 5.1.1 --depth 1 https://github.com/kokkos/kokkos-kernels.git \
 && git clone --branch develop --depth 1 https://github.com/kokkos/kokkos-tools.git \
 && cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_DEPRECATED_CODE_4=OFF -DKokkos_ENABLE_OPENMP=ON -B build -S kokkos \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/openmp/Kokkos \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/openmp -B build -S kokkos-fft \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/openmp/KokkosFFT \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/openmp -B build -S kokkos-kernels \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/openmp/KokkosKernels \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/openmp -B build -S kokkos-tools \
 && cmake --build build --parallel $(nproc) \
 && cmake --install build --prefix /opt/openmp/KokkosTools \
 && rm -rf build \
 && rm -rf kokkos kokkos-fft kokkos-kernels kokkos-tools

# GMGPolar
RUN git clone --branch v2.4.1 --depth 1 https://github.com/SciCompMod/GMGPolar.git \
 && for backend in serial openmp; do \
      cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/${backend} -B build -S GMGPolar \
   && cmake --build build --parallel $(nproc) \
   && cmake --install build --prefix /opt/${backend}/GMGPolar \
   && rm -rf build \
   || exit 1; \
    done \
 && rm -rf GMGPolar

# Koliop
RUN git clone --branch v0.2.0 --depth 1 https://gitlab.com/cines/code.gysela/libkoliop \
 && for backend in serial openmp; do \
      cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/${backend} -B build -S libkoliop \
   && cmake --build build --parallel $(nproc) \
   && cmake --install build --prefix /opt/${backend}/koliop \
   && rm -rf build \
   || exit 1; \
    done \
 && rm -rf libkoliop

# DDC
RUN git clone --branch v0.16.0 --depth 1 https://github.com/CExA-project/ddc.git \
 && for backend in serial openmp; do \
      cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/${backend} -DDDC_BUILD_EXAMPLES=OFF -DDDC_BUILD_TESTS=OFF -B build -S ddc \
   && cmake --build build --parallel $(nproc) \
   && cmake --install build --prefix /opt/${backend}/DDC \
   && rm -rf build \
   || exit 1; \
    done \
 && rm -rf ddc

USER ci:ci
WORKDIR /ci
ENV HOME=/ci
ENV PATH="/opt/venv/bin:$PATH"
ENV CMAKE_PREFIX_PATH="/opt/googletest"

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
