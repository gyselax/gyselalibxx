# SPDX-License-Identifier: MIT

FROM ubuntu:focal

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
 && echo "deb [ arch=amd64 ] https://raw.githubusercontent.com/pdidev/repo/pdi-master.ubuntu focal main" > /etc/apt/sources.list.d/pdi.list \
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
    pdiplugin-decl-hdf5-serial \
    pdiplugin-set-value \
    python3-pip \
    libpdi-dev \
 && apt-get purge -y \
    apt-transport-https \
    apt-utils \
    wget \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && python3 -m pip install --system -U --no-cache-dir \
    PyYAML \
    dask \
    h5py \
    imageio \
    matplotlib \
    numpy \
    pip \
    scipy \
    xarray \
    pylint \
 && rm -rf /var/lib/apt/lists/* \
 && useradd -d /ci -m -U ci \
 && mkdir /data \
 && chown ci:ci /ci /data

USER ci:ci
WORKDIR /ci
ENV HOME=/ci

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
