FROM ubuntu:focal

COPY bash_run /bin/
ENV BASH_ENV=/etc/profile
SHELL ["/bin/bash", "-c"]

RUN chmod +x /bin/bash_run \
 && apt-get update -y \
 && apt-get upgrade -y \
 && echo -e '13\n33\n' \
  | apt-get install -qy --no-install-recommends \
    tzdata \
 && apt-get install -y --no-install-recommends \
    build-essential \
    clang-format \
    cmake \
    git \
    libfftw3-dev \
    liblapack-dev \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && rm -rf /var/lib/apt/lists/* \
 && useradd -m -U default \
 && mkdir /data

USER 1001
WORKDIR /data

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
