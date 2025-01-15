#!/bin/bash

TOOLCHAIN_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

module purge

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-MI250"

# FIXME: This loads unneeded modules. It should not interfere with our build.
module load develop
module load spack-user-4.0.0

# Inject PDI recipes into our local repo.
git clone https://github.com/pdidev/spack pdi.spack || true
cd pdi.spack && git checkout ac5b78d && cd ..
cp -rf -- pdi.spack/packages "${SPACK_USER_PREFIX}/config_user_spack/local-repo"

# NOTE: A sparse checkout would be great.
git clone https://github.com/spack/spack spack.spack || true
cd spack.spack && git fetch && git checkout ba6cb62 && cd ..
# NOTE: We may be overriding some CINES modified recipes.
cp -rf -- spack.spack/var/spack/repos/builtin/packages/ginkgo "${SPACK_USER_PREFIX}/config_user_spack/local-repo/packages"

# FIXME: Add a recent Ginkgo to its package.py
sed -i '/version("master", branch="master")/a\    version("1.8.0", commit="586b1754058d7a32d4bd1b650f9603484c2a8927")  # v1.8.0' ${SPACK_USER_PREFIX}/config_user_spack/local-repo/packages/ginkgo/package.py

echo "Preparing the Spack environment..."

# # Either do not use the source mirrors or use --no-check-signature. There is an
# # ongoing issue a CINES to fix the issue.
# echo "# Disabled" >"${SPACK_USER_PREFIX}/config_user_spack/mirrors.yaml"

# We use GCC as a base compiler (c/c++/fortran) and implicitly, ROCm's hipcc when the +rocm variant is specified.
PRODUCT_SPEC_LIST="
ninja%gcc+re2c build_system=generic arch=linux-rhel8-zen3
libyaml%gcc build_system=autotools arch=linux-rhel8-zen3
paraconf%gcc~fortran~ipo~shared~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3
pdi%gcc~benchs~docs+fortran~ipo+python~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3
pdiplugin-decl-hdf5%gcc~benchs~fortran~ipo~mpi~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3
pdiplugin-set-value%gcc~ipo~tests build_system=cmake build_type=Release generator=ninja arch=linux-rhel8-zen3
pdiplugin-trace%gcc~ipo~tests build_system=cmake build_type=Release generator=ninja arch=linux-rhel8-zen3
ginkgo%gcc~cuda~develtools~full_optimizations~hwloc~ipo~mpi+openmp+rocm~sde~shared~sycl amdgpu_target=gfx90a build_system=cmake build_type=Release generator=ninja arch=linux-rhel8-zen3
eigen%gcc~ipo build_system=cmake build_type=Release generator=ninja arch=linux-rhel8-zen3
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
spack spec --reuse-deps ${PRODUCT_SPEC_LIST}
for ((i = 0; i < 2; ++i)); do
    spack install --no-check-signature --fresh ${PRODUCT_SPEC_LIST}
done

# Ensure we expose modules for every installed software.
spack module tcl refresh --delete-tree --yes-to-all ${PRODUCT_SPEC_LIST}
