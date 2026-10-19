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

# Defines GYSELA_SPACK_PREFIX and loads Spack
. "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")/config.sh"

if [[ -n "${SPACK_ENV:-}" ]]; then
    spack env deactivate
fi

spack env activate ${environment_name}

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=16

# Add Kokkos Tools to the `LD_LIBRARY_PATH`
kokkos_tools_dir="$(spack location --install-dir kokkos-tools)"
export LD_LIBRARY_PATH="${kokkos_tools_dir}/lib64${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
