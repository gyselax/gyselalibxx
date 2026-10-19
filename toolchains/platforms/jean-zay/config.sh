#!/bin/bash

if [ "${BASH_SOURCE[0]}" -ef "$0" ]
then
    echo "This script must be sourced not executed."
    echo ". $0"
    exit 1
fi

# `module` may be missing in non-interactive shells (e.g. under srun)
if command -v module >/dev/null 2>&1; then module purge; fi

export GYSELA_SPACK_GROUP="${IDRPROJ}"
export GYSELA_SPACK_JOBS=128
GYSELA_SPACK_VERSION="$(tr -d '[:space:]' < "$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]}")")/../../spack-version.txt")"
export GYSELA_SPACK_VERSION
export GYSELA_SPACK_PREFIX=${ALL_CCFRSCRATCH}/gysela-spack-${GYSELA_SPACK_VERSION}

export PYTHONPYCACHEPREFIX="${ALL_CCFRSCRATCH}/pycache"
export SPACK_DISABLE_LOCAL_CONFIG=true

if [ -f "${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh" ]; then
    . "${GYSELA_SPACK_PREFIX}"/share/spack/setup-env.sh
fi

# Build on GPU compute nodes (no internet there: sources are mirrored from the login node)
export GYSELA_USE_SLURM=1
export GYSELA_SLURM_CPUS_PER_GPU=24
export GYSELA_SLURM_QOS="qos_gpu_h100-dev"
export GYSELA_SLURM_TIME="02:00:00"
export GYSELA_SLURM_GPUS_CONCRETIZE=1
export GYSELA_SLURM_GPUS_INSTALL=3
# GYSELA_SLURM_ACCOUNT must be provided (env var or --slurm-account)
