#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

if [[ "$(id -gn)" != "gysela" ]]; then
    echo "Primary group must be 'gysela'!" >&2
    exit 1
fi

set -eu

module purge

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

export SPACK_VERSION=1.1.1
export SPACK_PATH=/data/gyselarunner/spack

if [ ! -d "$SPACK_PATH" ]; then
    git clone --branch "v${SPACK_VERSION}" --depth 1 https://github.com/spack/spack.git $SPACK_PATH
fi
. $SPACK_PATH/share/spack/setup-env.sh

spack config --scope site add 'config:install_tree:projections:all:"{compiler.name}-{compiler.version}/{name}-{version}-{hash}"'

spack config --scope site add 'packages:all:providers:blas:[openblas]'
spack config --scope site add 'packages:all:providers:lapack:[openblas]'
spack config --scope site add 'packages:git:version:[":2.46"]'

module load gcc/13
spack compiler find --scope site ${GCC_HOME}
module purge

spack env remove --yes-to-all gyselalibxx-env-omp-cuda
spack env create gyselalibxx-env-omp-cuda "${TOOLCHAIN_ROOT_DIRECTORY}/v100/gyselalibxx-spack-environment.yaml"

echo "Preparing the Spack environment..."

spack --env gyselalibxx-env-omp-cuda external find cuda
spack --env gyselalibxx-env-omp-cuda repo update
spack --env gyselalibxx-env-omp-cuda concretize --force
spack --env gyselalibxx-env-omp-cuda install --concurrent-packages 2 --jobs 16

spack env remove --yes-to-all gyselalibxx-env-omp
spack env create gyselalibxx-env-omp "${TOOLCHAIN_ROOT_DIRECTORY}/xeon/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-env-omp repo update
spack --env gyselalibxx-env-omp concretize --force
spack --env gyselalibxx-env-omp install --concurrent-packages 2 --jobs 16
