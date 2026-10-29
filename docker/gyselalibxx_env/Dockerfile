# SPDX-License-Identifier: MIT

FROM ubuntu:noble

COPY --chmod=0755 bash_run /bin/bash_run
ENV BASH_ENV=/etc/profile
ENV DEBIAN_FRONTEND=noninteractive
SHELL ["/bin/bash", "-c"]


RUN echo tzdata tzdata/Areas string Etc | debconf-set-selections \
 && echo tzdata tzdata/Zones/Etc string UTC | debconf-set-selections \
 && apt-get update -y \
 && apt-get upgrade -y \
 && apt-get install -y --no-install-recommends \
    apt-transport-https \
    apt-utils \
    ca-certificates \
    wget \
 && apt-get update -y \
 && apt-get install -y --no-install-recommends \
    build-essential \
    clang-18 \
    libomp-dev \
    clang-format-14 \
    git \
    openssh-client \
    python3-pip \
    python3-venv \
    cmake \
 && ln -sf /usr/bin/clang-18 /usr/bin/clang \
 && ln -sf /usr/bin/clang++-18 /usr/bin/clang++ \
 && git clone --branch 2.22.0 --depth 1 https://github.com/cppcheck-opensource/cppcheck.git \
 && cd cppcheck \
 && cmake -S . -B build -DUSE_MATCHCOMPILER=ON -DCMAKE_BUILD_TYPE=Release \
 && cmake --build build \
 && cmake --install build \
 && cd .. \
 && rm -rf cppcheck \
 && apt-get purge -y \
    apt-transport-https \
    apt-utils \
    cmake \
    wget \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && python3 -m venv /opt/venv \
 && /opt/venv/bin/python -m pip install -U --no-cache-dir \
    pip \
    'pylint==3.3.9' \
 && rm -rf /var/lib/apt/lists/* \
 && useradd -d /ci -m -U ci \
 && mkdir /data \
 && chown -R ci:ci /ci /data /opt/venv

# Install what is necessary for Spack
# See https://spack.readthedocs.io/en/latest/installing_prerequisites.html#linux
RUN apt-get update --yes \
 && apt-get install --yes --no-install-recommends \
    file \
    bzip2 \
    ca-certificates \
    g++ \
    gcc \
    gfortran \
    git \
    gzip \
    lsb-release \
    patch \
    python3 \
    tar \
    unzip \
    xz-utils \
    zstd

USER ci:ci
WORKDIR /ci
ENV HOME=/ci
ENV PATH="/opt/venv/bin:$PATH"

COPY spack-env /data/spack-env
RUN --mount=type=secret,id=github_user_id,env=GITHUB_USER \
    --mount=type=secret,id=github_token_id,env=GITHUB_TOKEN \
    git clone --branch v1.2.2 --depth 1 https://github.com/spack/spack.git /ci/spack \
    && . /ci/spack/share/spack/setup-env.sh \
    && for backend in serial openmp; do \
        spack env create "gyselalibxx-${backend}" "/data/spack-env/spack-${backend}.yaml" \
        && spack --env "gyselalibxx-${backend}" config add config:install_tree:padded_length:128 \
        && spack --env "gyselalibxx-${backend}" mirror set --autopush github-buildcache \
        && spack --env "gyselalibxx-${backend}" repo update \
        && spack --env "gyselalibxx-${backend}" install --fail-fast; \
    done

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
