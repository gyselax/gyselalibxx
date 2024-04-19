#!/bin/bash

TOOLCHAIN_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

module purge

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/spack-install-MI250"

# FIXME: Currently in development.
module load develop
# FIXME: This loads unneeded modules. It should not interfere with out build.
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

echo "Preparing the Spack environment..."

# # Either do not use the source mirrors or use --no-check-signature. There is an
# # ongoing issue a CiINES to fix the issue.
# echo "# Disabled" >"${SPACK_USER_PREFIX}/config_user_spack/mirrors.yaml"

# NOTE: We do two passes because Spack is bugged and even though we specified
# ninja/cmake as explicit target it does NOT register it that way, thus failing
# to generate modules for it.
for ((i = 0; i < 2; ++i)); do
    # GCC based with hipcc for the ROCm variant:
    spack install --no-check-signature --reuse-deps \
        libyaml@0.2.5%gcc@12.1 build_system=autotools arch=linux-rhel8-zen3 \
        paraconf@1.0.0%gcc@12.1~fortran~ipo~shared~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3 \
        pdi@1.6.0%gcc@12.1~benchs~docs+fortran~ipo+python~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3 \
        pdiplugin-decl-hdf5@1.6.0%gcc@12.1~benchs~fortran~ipo~mpi~tests build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3 \
        ginkgo@1.7.0%gcc@12.1~cuda~develtools~full_optimizations~hwloc~ipo~mpi+openmp+rocm~sde~shared~sycl amdgpu_target=gfx90a build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3 \
        eigen@3.4.0%gcc@12.1~ipo build_system=cmake build_type=Release generator==ninja arch=linux-rhel8-zen3 \
        cmake@3.27.7%gcc@12.1~doc+ncurses+ownlibs build_system=generic build_type=Release arch=linux-rhel8-zen3 \
        ninja@1.11.1%gcc@12.1+re2c build_system=generic arch=linux-rhel8-zen3

    # Ensure we expose modules for every installed software.
    spack module tcl refresh --delete-tree --yes-to-all
done
