#!/bin/bash

TOOLCHAIN_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

module purge

#set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

#export SPACK_USER_PREFIX="${WORK}/gysela/spack-install-A100"
export SPACK_USER_PREFIX="${HOME}/gysela/spack-install-A100"

mkdir -p ${SPACK_USER_PREFIX}/config_user_spack/local-repo
mkdir -p ${SPACK_USER_PREFIX}/config_user_spack/local-repo/packages

module load spack
module load gcc/12.2.0

# Inject PDI recipes into our local repo.
git clone  https://github.com/pdidev/spack pdi.spack || true
cd pdi.spack; git checkout a318a62631b054adca836f8267ca8905f3aa6d73; cd ..

# NOTE: A sparse checkout would be great.
git clone --depth=1 --branch=v0.22.2 https://github.com/spack/spack spack_0.22.2 || true

# FIXME: Add Ginkgo 1.8.0 to its package.py
sed -i '/version("master", branch="master")/a\    version("1.8.0", commit="586b1754058d7a32d4bd1b650f9603484c2a8927")  # v1.8.0' spack_0.22.2/var/spack/repos/builtin/packages/ginkgo/package.py

. spack_0.22.2/share/spack/setup-env.sh
spack repo add pdi.spack

echo "Preparing the Spack environment..."

spack env create -d ${HOME}/gyselalib_env ${TOOLCHAIN_ROOT_DIRECTORY}/spack.yaml
cd ${HOME}
spack -e gyselalib_env install
