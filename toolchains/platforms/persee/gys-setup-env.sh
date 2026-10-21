#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

: "${GYSELA_SPACK_BASE:=/data/gyselarunner}"

export GYSELA_SPACK_GROUP="gysela"
export GYSELA_SPACK_JOBS=32
GYSELA_SPACK_VERSION="$(tr -d '[:space:]' < "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")/../../spack-version.txt")"
export GYSELA_SPACK_VERSION
export GYSELA_SPACK_PREFIX="${GYSELA_SPACK_BASE}/gysela-spack-${GYSELA_SPACK_VERSION}"

export SPACK_DISABLE_LOCAL_CONFIG=true

if [ -f "${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh" ]; then
    . "${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh"
fi
