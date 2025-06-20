#!/bin/bash

SPACK_USER_VERSION="spack-user-4.0.0"

export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-MI250/Configuration.${SPACK_USER_VERSION}"
export SPACK_USER_CACHE_PATH="${SPACK_USER_PREFIX}/cache"

module purge

module load "${SPACK_USER_VERSION}"
which spack
spack debug report
# Spack must work in a clean, purged environment so it can load modules without
# having to purge itself or clearing environment variables (which it does not
# do..). When we spack env activate, the same constraint applies.
# Use spack load instead of an environment activation as it should limit the
# inode produced by the environment's view.
# eval -- "$(spack env activate --prompt --sh gyselalibxx-spack-environment)"
# unalias despacktivate
# unset despacktivate
# function despacktivate() {
#     eval "$(spack env deactivate --sh)"
# }

eval -- "$(
    spack \
        --env gyselalibxx-spack-environment \
        load --sh \
        cmake \
        fftw \
        ginkgo \
        googletest \
        ninja \
        paraconf \
        pdi \
        pdiplugin-decl-hdf5 \
        pdiplugin-decl-netcdf \
        pdiplugin-mpi \
        pdiplugin-set-value \
        pdiplugin-trace \
        python \
        py-dask \
        py-h5py \
        py-imageio \
        py-matplotlib \
        py-netcdf4 \
        py-numpy \
        py-scipy \
        py-sympy \
        py-xarray \
        py-pyyaml
)"

# Due to https://github.com/gyselax/gyselalibxx/pull/198#issuecomment-2943081411
# we have a different spack to build dependencies and to build gyselalib. It is
# fine except we must unset some environment variables defined by spack's HIP
# packages so they do not interfere with HIPCC/clang in the 6.3.3 version.
# Notably:
unset ROCM_PATH                  # /opt/rocm-6.1.2
unset HIP_CLANG_PATH             # /opt/rocm-6.1.2/llvm/bin
unset HSA_PATH                   # /opt/rocm-6.1.2
unset ROCMINFO_PATH              # /opt/rocm-6.1.2
unset DEVICE_LIB_PATH            # /opt/rocm-6.1.2/amdgcn/bitcode
unset HIP_DEVICE_LIB_PATH        # /opt/rocm-6.1.2/amdgcn/bitcode
unset LLVM_PATH                  # /opt/rocm-6.1.2/llvm
unset COMGR_PATH                 # /opt/rocm-6.1.2
unset HIPCC_COMPILE_FLAGS_APPEND # --rocm-path=/opt/rocm-6.1.2 --gcc-toolchain=/opt/rh/gcc-toolset-13/root/usr --rocm-path=/opt/rocm-6.1.2 --gcc-toolchain=/opt/rh/gcc-toolset-13/root/usr
unset HIPFLAGS                   # --gcc-toolchain=/opt/rh/gcc-toolset-13/root/usr --gcc-toolchain=/opt/rh/gcc-toolset-13/root/usr

module load cpe/24.07
module load craype-x86-trento craype-accel-amd-gfx90a
# NOTE: Force 6.3.3 due to startup failures (https://github.com/gyselax/gyselalibxx/pull/198#issuecomment-2943081411)
module load PrgEnv-gnu-amd amd-mixed/6.3.3

module list
