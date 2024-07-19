#!/bin/bash

TOOLCHAIN_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

module purge

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/spack-install-MI250"

# FIXME: This loads unneeded modules. It should not interfere with our build.
module load spack-MI250-3.1.0

# Inject PDI recipes into our local repo.
git clone https://github.com/pdidev/spack pdi.spack || true
cd pdi.spack && git pull && cd ..
cp -rf -- pdi.spack/packages "${SPACK_USER_PREFIX}/config_user_spack/local-repo"

# NOTE: A sparse checkout would be great.
git clone https://github.com/spack/spack spack.spack || true
cd spack.spack && git pull && cd ..
# NOTE: We may be overriding some CINES modified recipes.
cp -rf -- spack.spack/var/spack/repos/builtin/packages/ginkgo "${SPACK_USER_PREFIX}/config_user_spack/local-repo/packages"

# FIXME: Add Ginkgo 1.8.0 to its package.py
sed -i '/version("master", branch="master")/a\    version("1.8.0", commit="93432abadfd5be0ba8c931c220be9cd4797a5aca")  # v1.8.0' ${SPACK_USER_PREFIX}/config_user_spack/local-repo/packages/ginkgo/package.py

echo "Preparing the Spack environment..."

# # Either do not use the source mirrors or use --no-check-signature. There is an
# # ongoing issue a CINES to fix the issue.
# echo "# Disabled" >"${SPACK_USER_PREFIX}/config_user_spack/mirrors.yaml"

PRODUCT_SPEC_LIST="
libyaml@0.2.5%gcc@12.1 build_system=autotools arch=linux-rhel8-zen3
paraconf@1.0.0%gcc@12.1~fortran~ipo~shared~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3
pdi@1.6.0%gcc@12.1~benchs~docs+fortran~ipo+python~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3
pdiplugin-decl-hdf5@1.6.0%gcc@12.1~benchs~fortran~ipo~mpi~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3
pdiplugin-set-value@1.6.0%gcc@12.1~ipo~tests build_system=cmake build_type=Release generator=ninja arch=linux-rhel8-zen3
pdiplugin-trace@1.6.0%gcc@12.1~ipo~tests build_system=cmake build_type=Release generator=ninja arch=linux-rhel8-zen3
ginkgo@1.8.0%gcc@12.1~cuda~develtools~full_optimizations~hwloc~ipo~mpi+openmp+rocm~sde~shared~sycl amdgpu_target=gfx90a build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3
eigen@3.4.0%gcc@12.1~ipo build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3
ninja@1.11.1%gcc@12.1+re2c build_system=generic arch=linux-rhel8-zen3
"

# If we start preparing a new environment, ensure we wont get name clashes by
# uninstalling previous products.
# NOTE: This may fail if the product are not already installed. IMO this is a
# bug, a rm of something that does not exists is a success not a failure.
spack uninstall --dependents --all --yes-to-all \
    ${PRODUCT_SPEC_LIST} || true

# GCC based with hipcc for the ROCm variant:
spack install --no-check-signature --reuse-deps \
    ${PRODUCT_SPEC_LIST}

# Ensure we expose modules for every installed software.
spack module tcl refresh --delete-tree --yes-to-all
