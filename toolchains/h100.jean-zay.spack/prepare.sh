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

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"
CPUS_PER_GPU=24

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

export SPACK_USER_PREFIX=$ALL_CCFRWORK/spack-user-install
export SPACK_USER_CONFIG_PATH=$SPACK_USER_PREFIX/configuration
export SPACK_USER_CACHE_PATH=$SPACK_USER_PREFIX/cache

# Avoid too many temporary files in the Spack installation tree
export PYTHONPYCACHEPREFIX=$ALL_CCFRSCRATCH/pycache

module purge

git clone --branch v1.1.0 --depth 1 https://github.com/spack/spack.git $ALL_CCFRWORK/spack || true
. $ALL_CCFRWORK/spack/share/spack/setup-env.sh

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

# Bootstrap must happen before getting on a compute node
spack bootstrap now

echo "Preparing the Spack environment..."

# Concretize on the compute node
GPUS=1
CPUS_PER_TASK=$((${GPUS} * ${CPUS_PER_GPU}))
srun --account=${slurm_account} --constraint=h100 --ntasks=1 --gres=gpu:${GPUS} --cpus-per-task=${CPUS_PER_TASK} --time=02:00:00 --qos=qos_gpu_h100-dev spack --env gyselalibxx-spack-environment concretize --force --jobs ${CPUS_PER_TASK}

# Get sources on the login node
spack --env gyselalibxx-spack-environment mirror create --all --directory ${SPACK_USER_CONFIG_PATH}/mirror
spack --env gyselalibxx-spack-environment mirror add local_filesystem ${SPACK_USER_CONFIG_PATH}/mirror

# Install on the compute node
GPUS=3
CPUS_PER_TASK=$((${GPUS} * ${CPUS_PER_GPU}))
srun --account=${slurm_account} --constraint=h100 --ntasks=1 --gres=gpu:${GPUS} --cpus-per-task=${CPUS_PER_TASK} --time=02:00:00 --qos=qos_gpu_h100-dev spack --env gyselalibxx-spack-environment install --concurrent-packages ${GPUS} --jobs ${CPUS_PER_GPU}
