#!/bin/bash

# Ensures the script is not being sourced
if [[ "${BASH_SOURCE[0]}" != "${0}" ]]; then
    echo "This script must be executed, not sourced!" >&2
    return 1
fi

TOOLCHAIN_ROOT_DIRECTORY="$(dirname -- "$(readlink -f -- "${BASH_SOURCE[0]:-${0}}")")"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-MI250/Configuration.${SPACK_USER_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/cache"

module purge
module load "${SPACK_USER_VERSION}"
module list

which spack
spack debug report

mkdir -p -- "${SPACK_USER_CONFIG_PATH}/external-repositories"

# Inject recent PDI recipes into our repository.
git clone https://github.com/pdidev/spack "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi.spack" || true
cd -- "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi.spack" && git fetch && git checkout 5483cfea7d2d39d654c2962114248f597b3ecf46 && cd -

git clone https://github.com/spack/spack spack.spack || true
cd spack.spack && git fetch && git checkout 8e7489b && cd ..
cp -rf -- spack.spack/var/spack/repos/builtin/packages/ginkgo "${SPACK_USER_CONFIG_PATH}/local-repo/packages"

sed -i '/args.append("-DGINKGO_HIP_AMDGPU={0}".format(arch_str))/a\                args.append("-DCMAKE_HIP_ARCHITECTURES={0}".format(arch_str))' "${SPACK_USER_CONFIG_PATH}/local-repo/packages/ginkgo/package.py"

spack env remove --yes-to-all gyselalibxx-spack-environment
spack env create gyselalibxx-spack-environment "${TOOLCHAIN_ROOT_DIRECTORY}/gyselalibxx-spack-environment.yaml"

spack --env gyselalibxx-spack-environment repo add "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi.spack"

echo "Preparing the Spack environment..."

spack --env gyselalibxx-spack-environment concretize --force

for ((i = 0; i < 8; ++i)); do
    spack --env gyselalibxx-spack-environment install --jobs 48 &
done

wait

# spack --env gyselalibxx-spack-environment module tcl refresh --delete-tree --yes-to-all
