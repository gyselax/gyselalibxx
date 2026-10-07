#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

if [[ $# -lt 1 ]]; then
    echo "Usage: . environment.sh <environment_name>" >&2
    return 1
fi
environment_name="$1"

module purge

SPACK_VERSION="1.2.2"
export SPACK_PREFIX=${ALL_CCFRSCRATCH}/gysela-spack-${SPACK_VERSION}
export SPACK_DISABLE_LOCAL_CONFIG=true

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=$ALL_CCFRSCRATCH/pycache

. ${SPACK_PREFIX}/share/spack/setup-env.sh
spack env activate ${environment_name}

# Add Kokkos Tools to the `LD_LIBRARY_PATH`
export LD_LIBRARY_PATH="$(spack location -i kokkos-tools)/lib64:$LD_LIBRARY_PATH"
