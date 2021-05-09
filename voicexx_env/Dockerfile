FROM ubuntu:focal as prebuild

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
    liblapack-dev \
 && apt-get autoremove -y \
 && apt-get clean -y \
 && apt-get autoclean -y \
 && rm -rf /var/lib/apt/lists/*

ENV HOME /home/default
RUN useradd -d ${HOME} -m -u 1001 -U default
USER 1001
WORKDIR ${HOME}

ENTRYPOINT ["/bin/bash_run"]
CMD ["/bin/bash", "-li"]
