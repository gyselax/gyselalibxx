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

git clone --branch v0.23.1 --depth 1 https://github.com/spack/spack.git $ALL_CCFRWORK/spack || true
. $ALL_CCFRWORK/spack/share/spack/setup-env.sh

mkdir -p -- "${SPACK_USER_CONFIG_PATH}/external-repositories"

# Inject recent PDI recipes into our repository.
git clone https://github.com/pdidev/spack "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi" || true
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi" fetch
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi" checkout 5483cfea7d2d39d654c2962114248f597b3ecf46

git clone https://github.com/gyselax/spack "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx" || true
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx" fetch
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx" checkout releases/v0.23
git -C "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx" pull

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-spack-environment repo add "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi"
spack --env gyselalibxx-spack-environment repo add "${SPACK_USER_CONFIG_PATH}/external-repositories/gyselalibxx"

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
srun --account=${slurm_account} --constraint=h100 --ntasks=1 --gres=gpu:${GPUS} --cpus-per-task=${CPUS_PER_TASK} --time=02:00:00 --qos=qos_gpu_h100-dev spack --env gyselalibxx-spack-environment install --jobs ${CPUS_PER_TASK}
