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

if command -v spack >/dev/null 2>&1
then
    spack env deactivate
else
    . /data/gyselarunner/gysela-spack-1.2.2/share/spack/setup-env.sh
fi

spack env activate ${environment_name}

export OMP_PROC_BIND=spread
export OMP_PLACES=threads
export OMP_NUM_THREADS=16

# Add Kokkos Tools to the `LD_LIBRARY_PATH`
export LD_LIBRARY_PATH="$(spack --env ${environment_name} location --install-dir kokkos-tools)/lib64:$LD_LIBRARY_PATH"
