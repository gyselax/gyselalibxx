# SPDX-License-Identifier: MIT

FROM ubuntu:jammy

COPY bash_run /bin/
ENV BASH_ENV=/etc/profile
ENV PATH=/opt/cmake-3.27/bin:${PATH:+:${PATH}}
SHELL ["/bin/bash", "-c"]


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
 && echo "deb [ arch=amd64 ] https://repo.pdi.dev/ubuntu jammy main" > /etc/apt/sources.list.d/pdi.list \
 && wget -O /etc/apt/trusted.gpg.d/pdidev-archive-keyring.gpg https://repo.pdi.dev/ubuntu/pdidev-archive-keyring.gpg \
 && chmod a+r /etc/apt/trusted.gpg.d/pdidev-archive-keyring.gpg \
 && apt-get update -y \
 && apt-get install -y --no-install-recommends \
    build-essential \
    clang-format \
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
    libpdi-dev \
 && wget https://cmake.org/files/v3.27/cmake-3.27.9-linux-x86_64.tar.gz \
 && tar --extract --file cmake-3.27.9-linux-x86_64.tar.gz \
 && mv cmake-3.27.9-linux-x86_64 /opt/cmake-3.27 \
 && git clone -b v1.8.0 https://github.com/ginkgo-project/ginkgo.git \
 && cd ginkgo \
 && cmake -DCMAKE_BUILD_TYPE=Release -DGINKGO_BUILD_TESTS=OFF -DGINKGO_BUILD_EXAMPLES=OFF -DGINKGO_BUILD_BENCHMARKS=OFF -B build \
 && cmake --build build \
 && cmake --install build \
 && cd .. \
 && rm -rf ginkgo \
 && git clone --branch 4.4.01 --depth 1 https://github.com/kokkos/kokkos.git \
 && git clone --branch develop-spline-kernels-v2 --depth 1 https://github.com/yasahi-hpc/kokkos-kernels \
 && git clone --branch develop --depth 1 https://github.com/kokkos/kokkos-tools.git \
 && git clone --branch v0.2.1 --depth 1 https://github.com/kokkos/kokkos-fft.git \
 && export Kokkos_ROOT=/opt/serial/kokkos \
 && export KokkosFFT_ROOT=/opt/serial/kokkos-fft \
 && export KokkosKernels_ROOT=/opt/serial/kokkos-kernels \
 && export KokkosTools_ROOT=/opt/serial/kokkos-tools \
 && cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_DEPRECATED_CODE_4=OFF -DKokkos_ENABLE_SERIAL=ON -B build -S kokkos \
 && cmake --build build \
 && cmake --install build --prefix $Kokkos_ROOT \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -B build -S kokkos-fft \
 && cmake --build build \
 && cmake --install build --prefix $KokkosFFT_ROOT \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -B build -S kokkos-kernels \
 && cmake --build build \
 && cmake --install build --prefix $KokkosKernels_ROOT \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -B build -S kokkos-tools \
 && cmake --build build \
 && cmake --install build --prefix $KokkosTools_ROOT \
 && rm -rf build \
 && export Kokkos_ROOT=/opt/openmp/kokkos \
 && export KokkosFFT_ROOT=/opt/openmp/kokkos-fft \
 && export KokkosKernels_ROOT=/opt/openmp/kokkos-kernels \
 && export KokkosTools_ROOT=/opt/openmp/kokkos-tools \
 && cmake -DCMAKE_BUILD_TYPE=Release -DKokkos_ENABLE_DEPRECATED_CODE_4=OFF -DKokkos_ENABLE_OPENMP=ON -B build -S kokkos \
 && cmake --build build \
 && cmake --install build --prefix $Kokkos_ROOT \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -B build -S kokkos-fft \
 && cmake --build build \
 && cmake --install build --prefix $KokkosFFT_ROOT \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -B build -S kokkos-kernels \
 && cmake --build build \
 && cmake --install build --prefix $KokkosKernels_ROOT \
 && rm -rf build \
 && cmake -DCMAKE_BUILD_TYPE=Release -B build -S kokkos-tools \
 && cmake --build build \
 && cmake --install build --prefix $KokkosTools_ROOT \
 && rm -rf build kokkos kokkos-fft kokkos-kernels kokkos-tools \
 && git clone -b 2.14.2 https://github.com/danmar/cppcheck.git \
 && cd cppcheck \
 && cmake -S . -B build \
 && cmake --build build \
 && cmake --install build \
 && cd .. \
 && rm -rf cppcheck \
 && apt-get purge -y \
    apt-transport-https \
    apt-utils \
    wget \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && python3 -m pip install -U --no-cache-dir \
    PyYAML \
    dask \
    h5py \
    imageio \
    matplotlib \
    numpy \
    pip \
    scipy \
    xarray \
    'pylint==2.17.4' \
 && rm -rf /var/lib/apt/lists/* \
 && useradd -d /ci -m -U ci \
 && mkdir /data \
 && chown ci:ci /ci /data

USER ci:ci
WORKDIR /ci
ENV HOME=/ci

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
