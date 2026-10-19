#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

module purge

export GYSELA_SPACK_GROUP="project_465002992"
export GYSELA_SPACK_JOBS=128
export GYSELA_SPACK_VERSION="1.2.2"
export GYSELA_SPACK_PREFIX=/scratch/project_465002992/gysela-spack-${GYSELA_SPACK_VERSION}

export PYTHONPYCACHEPREFIX="/flash/project_465002992/pycache"
export SPACK_DISABLE_LOCAL_CONFIG=true

if [ -f "${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh" ]; then
    . ${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh
fi
