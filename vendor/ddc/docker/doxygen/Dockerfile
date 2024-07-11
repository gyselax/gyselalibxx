FROM ubuntu:noble

# Copyright (C) The DDC development team, see COPYRIGHT.md file
#
# SPDX-License-Identifier: MIT

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
 && mkdir --parents --mode=0755 /etc/apt/keyrings \
 && echo "deb [arch=amd64 signed-by=/etc/apt/keyrings/pdidev-archive-keyring.gpg] https://raw.githubusercontent.com/pdidev/repo/ubuntu noble main" > /etc/apt/sources.list.d/pdi.list \
 && wget -q -O /etc/apt/keyrings/pdidev-archive-keyring.gpg https://raw.githubusercontent.com/pdidev/repo/ubuntu/pdidev-archive-keyring.gpg \
 && apt-get update -y \
 && apt-get install -y --no-install-recommends \
    build-essential \
    clang \
    cm-super-minimal \
    cmake \
    cmake-data \
    libfftw3-dev \
    liblapacke-dev \
    libpdi-dev \
    ghostscript \
    git \
    graphviz \
    pkg-config \
    doxygen-latex \
 && git clone -b v1.7.0 https://github.com/ginkgo-project/ginkgo.git \
 && cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DGINKGO_BUILD_OMP=OFF -DGINKGO_BUILD_TESTS=OFF -DGINKGO_BUILD_EXAMPLES=OFF -DGINKGO_BUILD_BENCHMARKS=OFF -B build-ginkgo -S ginkgo \
 && cmake --build build-ginkgo \
 && cmake --install build-ginkgo \
 && rm -rf ginkgo \
 && apt-get purge -y \
    apt-transport-https \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && rm -rf /var/lib/apt/lists/* \
 && useradd -d /data -m -U ci

USER ci:ci
WORKDIR /data

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
