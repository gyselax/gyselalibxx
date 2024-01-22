# SPDX-License-Identifier: MIT

FROM ubuntu:jammy

COPY bash_run /bin/
ENV BASH_ENV=/etc/profile
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
 && echo "deb [ arch=amd64 ] https://raw.githubusercontent.com/pdidev/repo/pdi-master.ubuntu jammy main" > /etc/apt/sources.list.d/pdi.list \
 && wget -O /etc/apt/trusted.gpg.d/pdidev-archive-keyring.gpg https://raw.githubusercontent.com/pdidev/repo/pdi-master.ubuntu/pdidev-archive-keyring.gpg \
 && chmod a+r /etc/apt/trusted.gpg.d/pdidev-archive-keyring.gpg \
 && apt-get update -y \
 && apt-get install -y --no-install-recommends \
    build-essential \
    clang-format \
    cmake \
    doxygen \
    pkg-config \
    cmake-data \
    git \
    libfftw3-dev \
    liblapack-dev \
    openssh-client \
    pdidev-archive-keyring \
    hdf5-tools \
    pdiplugin-decl-hdf5-serial \
    pdiplugin-set-value \
    python3-pip \
    libpdi-dev \
 && git clone -b v1.7.0 https://github.com/ginkgo-project/ginkgo.git \
 && cd ginkgo \
 && cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DGINKGO_BUILD_TESTS=OFF -DGINKGO_BUILD_EXAMPLES=OFF -DGINKGO_BUILD_BENCHMARKS=OFF \
 && cmake --build build \
 && cmake --install build \
 && cd .. \
 && rm -rf ginkgo \
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
