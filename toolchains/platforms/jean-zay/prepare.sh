#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

# Gets the slurm account to use for compiling on a compute node
if [ $# -ne 1 ]; then
    echo usage: ${0} h100_slurm_account
    exit 1
fi
slurm_account=${1}

set -eu

module purge

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"
CPUS_PER_GPU=24

GYSELA_SPACK_GROUP="ksw"
GYSELA_SPACK_JOBS=128
GYSELA_SPACK_VERSION="1.2.2"
GYSELA_SPACK_PREFIX=${ALL_CCFRSCRATCH}/gysela-spack-${GYSELA_SPACK_VERSION}
export SPACK_DISABLE_LOCAL_CONFIG=true
export PYTHONDONTWRITEBYTECODE=True

if [ ! -d "${GYSELA_SPACK_PREFIX}" ]; then
    mkdir --parents "${GYSELA_SPACK_PREFIX}"
    chmod g+s "${GYSELA_SPACK_PREFIX}"
    git clone --branch v${GYSELA_SPACK_VERSION} --depth 1 https://github.com/spack/spack.git "${GYSELA_SPACK_PREFIX}"
fi

. ${GYSELA_SPACK_PREFIX}/share/spack/setup-env.sh

# Bootstrap must happen before getting on a compute node
spack bootstrap now

for arch in h100; do
    env="gyselalibxx-${arch}"
    env_file="${TOOLCHAIN_ROOT_DIRECTORY}/${arch}/gyselalibxx-spack-environment.yaml"

    echo "Preparing the Spack environment ${env}"

    spack env remove --yes-to-all "${env}"
    spack env create "${env}" "${env_file}"

    spack env activate "${env}"
    spack repo update
    # Concretize on the compute node
    GPUS=1
    CPUS_PER_TASK=$((${GPUS} * ${CPUS_PER_GPU}))
    srun --account=${slurm_account} --constraint=h100 --ntasks=1 --gres=gpu:${GPUS} --cpus-per-task=${CPUS_PER_TASK} --time=02:00:00 --qos=qos_gpu_h100-dev spack concretize --quiet

    # Get sources on the login node
    spack mirror create --all --directory ${GYSELA_SPACK_PREFIX}/var/spack/mirror
    spack mirror add local_filesystem ${GYSELA_SPACK_PREFIX}/var/spack/mirror

    spack spec --install-status --namespaces

    # Install on the compute node
    GPUS=3
    CPUS_PER_TASK=$((${GPUS} * ${CPUS_PER_GPU}))
    srun --account=${slurm_account} --constraint=h100 --ntasks=1 --gres=gpu:${GPUS} --cpus-per-task=${CPUS_PER_TASK} --time=02:00:00 --qos=qos_gpu_h100-dev spack install --jobs ${CPUS_PER_TASK}

    spack env deactivate
done

find "${GYSELA_SPACK_PREFIX}" -user "${USER}" \
    -exec chgrp "${GYSELA_SPACK_GROUP}" {} + \
    -exec chmod g+rwX {} +
