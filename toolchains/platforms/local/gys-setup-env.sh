#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

# Personal workstation / laptop: no modules, no Slurm, no shared install.
# Everything can be overridden from the environment, e.g.
#   GYSELA_SPACK_BASE=/opt/gysela GYSELA_SPACK_JOBS=4 ./install.sh --platform local ...
: "${GYSELA_SPACK_BASE:=${XDG_DATA_HOME:-${HOME}/.local/share}/gysela}"

GYSELA_SPACK_GROUP="$(id -gn)"
export GYSELA_SPACK_GROUP
export GYSELA_SPACK_JOBS="${GYSELA_SPACK_JOBS:-$(nproc 2>/dev/null || echo 4)}"
GYSELA_SPACK_VERSION="$(tr -d '[:space:]' < "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")/../../spack-version.txt")"
export GYSELA_SPACK_VERSION
export GYSELA_SPACK_PREFIX="${GYSELA_SPACK_BASE}/gysela-spack-${GYSELA_SPACK_VERSION}"

# There is no site-specific externals.yaml: let Spack detect these system packages.
export GYSELA_DETECT_EXTERNALS="gcc"

export SPACK_DISABLE_LOCAL_CONFIG=true

if [ -f "${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh" ]; then
    . "${GYSELA_SPACK_PREFIX}"/share/spack/setup-env.sh
fi
