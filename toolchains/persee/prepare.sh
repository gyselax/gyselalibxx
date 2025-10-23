#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

set -eu

module purge

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

SPACK_VERSION="0.23.1"

export SPACK_PREFIX=/data/gyselarunner/spack-${SPACK_VERSION}

cd /tmp
wget https://github.com/spack/spack/releases/download/v${SPACK_VERSION}/spack-${SPACK_VERSION}.tar.gz
tar -xvf spack-${SPACK_VERSION}.tar.gz
rm spack-${SPACK_VERSION}.tar.gz
mv /tmp/spack-${SPACK_VERSION} ${SPACK_PREFIX}

. ${SPACK_PREFIX}/share/spack/setup-env.sh

spack config --scope site add 'config:install_tree:projections:all:"{compiler.name}-{compiler.version}/{name}-{version}-{hash}"'
spack config --scope site add 'config:connect_timeout:60'

spack config --scope site add 'packages:all:permissions:read:world'
spack config --scope site add 'packages:all:permissions:write:group'
spack config --scope site add 'packages:all:permissions:group:gysela'
spack config --scope site add 'packages:all:providers:blas:[openblas]'
spack config --scope site add 'packages:all:providers:lapack:[openblas]'
spack config --scope site add 'packages:git:version:[":2.46"]'

spack compiler find --scope site

mkdir -p -- "${SPACK_PREFIX}/external-repositories"

# Inject recent PDI recipes into our repository.
git clone https://github.com/pdidev/spack "${SPACK_PREFIX}/external-repositories/pdi" || true
git -C "${SPACK_PREFIX}/external-repositories/pdi" fetch
git -C "${SPACK_PREFIX}/external-repositories/pdi" checkout 5483cfea7d2d39d654c2962114248f597b3ecf46

git clone https://github.com/gyselax/spack "${SPACK_PREFIX}/external-repositories/gyselalibxx" || true
git -C "${SPACK_PREFIX}/external-repositories/gyselalibxx" fetch
git -C "${SPACK_PREFIX}/external-repositories/gyselalibxx" checkout releases/v0.23
git -C "${SPACK_PREFIX}/external-repositories/gyselalibxx" pull

spack repo add --scope site "${SPACK_PREFIX}/external-repositories/pdi"
spack repo add --scope site "${SPACK_PREFIX}/external-repositories/gyselalibxx"

spack env remove --yes-to-all gyselalibxx-env-omp-cuda
spack env create gyselalibxx-env-omp-cuda "${TOOLCHAIN_ROOT_DIRECTORY}/v100/gyselalibxx-spack-environment.yaml"

echo "Preparing the Spack environment..."

spack --env gyselalibxx-env-omp-cuda external find cuda
spack --env gyselalibxx-env-omp-cuda concretize --fresh --force
spack --env gyselalibxx-env-omp-cuda install --jobs 16

spack env remove --yes-to-all gyselalibxx-env-omp
spack env create gyselalibxx-env-omp "${TOOLCHAIN_ROOT_DIRECTORY}/xeon/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-env-omp concretize --fresh --force
spack --env gyselalibxx-env-omp install --jobs 16
