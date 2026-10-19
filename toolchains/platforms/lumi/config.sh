#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

# `module` may be missing in non-interactive shells (e.g. under srun)
if command -v module >/dev/null 2>&1; then module purge; fi

export GYSELA_SPACK_GROUP="project_465002992"
export GYSELA_SPACK_JOBS=128
GYSELA_SPACK_VERSION="$(tr -d '[:space:]' < "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")/../../spack-version.txt")"
export GYSELA_SPACK_VERSION
export GYSELA_SPACK_PREFIX=/scratch/project_465002992/gysela-spack-${GYSELA_SPACK_VERSION}

export PYTHONPYCACHEPREFIX="/flash/project_465002992/pycache"
export SPACK_DISABLE_LOCAL_CONFIG=true

if [ -f "${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh" ]; then
    . "${GYSELA_SPACK_PREFIX}"/share/spack/setup-env.sh
fi
