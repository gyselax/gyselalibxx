# SPDX-License-Identifier: MIT

FROM ubuntu:focal

LABEL "org.opencontainers.image.source"="https://github.com/Maison-de-la-Simulation/ddc"

ARG BACKEND

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
 && if [ "xhip" = "x${BACKEND}" ] \
 ;  then echo 'deb [arch=amd64] https://repo.radeon.com/amdgpu/latest/ubuntu focal main' > /etc/apt/sources.list.d/amdgpu.list \
      && echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/5.7.3 focal main' > /etc/apt/sources.list.d/rocm.list \
      && wget -q -O /etc/apt/trusted.gpg.d/rocm.asc https://repo.radeon.com/rocm/rocm.gpg.key \
      && echo -e 'Package: *\nPin: release o=repo.radeon.com\nPin-Priority: 600' > /etc/apt/preferences.d/rocm-pin-600 \
 ;  fi \
 && echo "deb [arch=amd64] https://raw.githubusercontent.com/pdidev/repo/ubuntu focal main" > /etc/apt/sources.list.d/pdi.list \
 && wget -q -O /etc/apt/trusted.gpg.d/pdidev-archive-keyring.gpg https://raw.githubusercontent.com/pdidev/repo/ubuntu/pdidev-archive-keyring.gpg \
 && chmod -R a+r /etc/apt/trusted.gpg.d/ \
 && apt-get update -y \
 && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    pkg-config \
    cmake-data \
    git \
    libelf-dev \
    libfftw3-dev \
    libhwloc-dev \
    libnuma1 \
    libpdi-dev \
    pdidev-archive-keyring \
 && case "${BACKEND}" in \
    "cpu") \
      apt-get install -y --no-install-recommends \
      clang \
 ;; "hip") \
      apt-get install -y --no-install-recommends \
      rocm-hip-sdk \
 ;; esac \
 && git clone -b v1.6.0 https://github.com/ginkgo-project/ginkgo.git \
 && cd ginkgo \
 && case "${BACKEND}" in \
    "cpu") \
      cmake -S . -B build -DCMAKE_BUILD_TYPE=RelWithDebInfo -DGINKGO_BUILD_OMP=ON -DGINKGO_BUILD_TESTS=OFF -DGINKGO_BUILD_EXAMPLES=OFF -DGINKGO_BUILD_BENCHMARKS=OFF \
 ;; "hip") \
      cmake -S . -B build -DCMAKE_PREFIX_PATH=/opt/rocm -DCMAKE_BUILD_TYPE=RelWithDebInfo -DCMAKE_HIP_ARCHITECTURES=gfx90a -DGINKGO_BUILD_HIP=ON -DGINKGO_BUILD_TESTS=OFF -DGINKGO_BUILD_EXAMPLES=OFF -DGINKGO_BUILD_BENCHMARKS=OFF \
 ;; esac \
 && cmake --build build \
 && cmake --install build \
 && cd .. \
 && rm -rf ginkgo \
 && apt-get purge -y \
    apt-transport-https \
    apt-utils \
    ca-certificates \
    wget \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && rm -rf /var/lib/apt/lists/* \
 && useradd -d /data -m -U ci \
 && if [ "xhip" = "x${BACKEND}" ] \
 ;  then echo 'export LD_LIBRARY_PATH="/opt/rocm/lib"' > /etc/profile.d/10-rocm.sh \
 ;  fi

USER ci:ci
WORKDIR /data

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
