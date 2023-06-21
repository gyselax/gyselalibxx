FROM ubuntu:jammy
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
 && echo 'deb [arch=amd64] https://repo.radeon.com/amdgpu/latest/ubuntu jammy main' > /etc/apt/sources.list.d/amdgpu.list \
 && echo 'deb [arch=amd64] https://repo.radeon.com/rocm/apt/latest jammy main' > /etc/apt/sources.list.d/rocm.list \
 && echo "deb [arch=amd64] https://raw.githubusercontent.com/pdidev/repo/ubuntu jammy main" > /etc/apt/sources.list.d/pdi.list \
 && wget -q -O /etc/apt/trusted.gpg.d/rocm.asc https://repo.radeon.com/rocm/rocm.gpg.key \
 && wget -q -O /etc/apt/trusted.gpg.d/pdidev-archive-keyring.gpg https://raw.githubusercontent.com/pdidev/repo/ubuntu/pdidev-archive-keyring.gpg \
 && chmod -R a+r /etc/apt/trusted.gpg.d/ \
 && echo -e 'Package: *\nPin: release o=repo.radeon.com\nPin-Priority: 600' > /etc/apt/preferences.d/rocm-pin-600 \
 && apt-get update -y \
 && apt-get install -y --no-install-recommends \
    build-essential \
    clang \
    clang-format \
    cmake \
    gcc-10 \
    g++-10 \
    gfortran-10 \
    git \
    libfftw3-dev \
    libpdi-dev \
    llvm-14-tools \
    nvidia-cuda-toolkit \
    pdidev-archive-keyring \
    pkg-config \
    rocm-hip-sdk \
 && apt-get purge -y \
    apt-transport-https \
    apt-utils \
    ca-certificates \
    wget \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && rm -rf /var/lib/apt/lists/* \
 && useradd -d /data -m -U ci

USER ci:ci
WORKDIR /data

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
