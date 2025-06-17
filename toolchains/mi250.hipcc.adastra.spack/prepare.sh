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

spack config --scope user update --yes-to-all packages
spack config --scope user update --yes-to-all config
spack config --scope user add 'packages:all:permissions:read:world'
spack config --scope user add 'packages:all:permissions:write:group'
spack config --scope user add 'concretizer:unify:true'
# spack config --scope user add 'modules:default:tcl:hash_length:4'
spack debug report

# Inject PDI recipes into our repos.
mkdir -p -- "${SPACK_USER_CONFIG_PATH}/external-repositories"
git clone https://github.com/pdidev/spack "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi.spack" || true
cd -- "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi.spack" && git fetch && git checkout 5483cfea7d2d39d654c2962114248f597b3ecf46 && cd -
spack repo add --scope user "${SPACK_USER_CONFIG_PATH}/external-repositories/pdi.spack" || true

# NOTE: A sparse checkout would be great.
git clone https://github.com/spack/spack spack.spack || true
cd spack.spack && git fetch && git checkout 8e7489b && cd ..
# NOTE: We may be overriding some CINES modified recipes.
cp -rf -- spack.spack/var/spack/repos/builtin/packages/ginkgo "${SPACK_USER_CONFIG_PATH}/local-repo/packages"

# FIXME: Add a recent Ginkgo to its package.py
sed -i '/args.append("-DGINKGO_HIP_AMDGPU={0}".format(arch_str))/a\                args.append("-DCMAKE_HIP_ARCHITECTURES={0}".format(arch_str))' "${SPACK_USER_CONFIG_PATH}/local-repo/packages/ginkgo/package.py"

echo "Preparing the Spack environment..."

# # Either do not use the source mirrors or use --no-check-signature. There is an
# # ongoing issue a CINES to fix the issue.
# echo "# Disabled" >"${SPACK_USER_CONFIG_PATH}/mirrors.yaml"

# We use GCC as a base compiler (c/c++/fortran) and implicitly, ROCm's hipcc when the +rocm variant is specified.
PRODUCT_SPEC_LIST="
ninja%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
libyaml%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
paraconf%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
pdi%gcc@13.2.1.mi250+python arch=linux-rhel8-zen3
pdiplugin-decl-hdf5%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
pdiplugin-set-value%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
pdiplugin-trace%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
pdiplugin-mpi%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
ginkgo%gcc@13.2.1.mi250+openmp+rocm~shared amdgpu_target=gfx90a arch=linux-rhel8-zen3
py-matplotlib%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
py-xarray%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
py-h5py%gcc@13.2.1.mi250 arch=linux-rhel8-zen3
"
# openblas@0.3.26%gcc@12.1.generic~bignuma~consistent_fpcsr+dynamic_dispatch+fortran~ilp64+locking+pic+shared build_system=makefile symbol_suffix=none threads=none arch=linux-rhel8-zen3

which spack

# If we start preparing a new environment, ensure we wont get name clashes by
# uninstalling previous products.
# NOTE: This may fail if the product are not already installed. IMO this is a
# bug, a rm of something that does not exists is a success not a failure.
echo "Removing old packages (errors are expected)."
spack uninstall --dependents --all --yes-to-all ${PRODUCT_SPEC_LIST} || true
echo "Installing new packages (errors are NOT expected)."
spack spec -lt ${PRODUCT_SPEC_LIST}
for ((i = 0; i < 2; ++i)); do
    spack install --no-check-signature --fresh ${PRODUCT_SPEC_LIST}
done

# Ensure we expose modules for every installed software.
spack module tcl refresh --delete-tree --yes-to-all ${PRODUCT_SPEC_LIST}
