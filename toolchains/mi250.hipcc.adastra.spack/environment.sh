ENVIRONMENT_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

module purge

module load cpe/23.12
# FIXME:
# craype-accel-amd-gfx90a Error with the cray wrappers clang: error: unsupported option '-fopenmp-targets=' for language mode 'HIP'
module load craype-x86-trento
module load PrgEnv-cray-amd amd-mixed/5.7.1

module load rocm/5.7.1
module load cray-fftw
module load cray-hdf5-parallel
module load cray-python

# FIXME: SPACK_USER_PREFIX path too long
# export SPACK_USER_PREFIX="${ENVIRONMENT_ROOT_DIRECTORY}/spack-install-MI250"
export SPACK_USER_PREFIX="${SHAREDWORKDIR}/spack-install-MI250"
module use "${SPACK_USER_PREFIX}/modules/tcl/linux-rhel8-zen3"

module load \
    gcc/12.1/zen3/libyaml/0.2.5 \
    gcc/12.1/zen3/paraconf/1.0.0 \
    gcc/12.1/zen3/pdi/1.6.0 \
    gcc/12.1/zen3/pdiplugin-decl-hdf5/1.6.0 \
    gcc/12.1/zen3/ginkgo/1.7.0 \
    gcc/12.1/zen3/eigen/3.4.0 \
    gcc/12.1/zen3/cmake/3.27.7 \
    gcc/12.1/zen3/ninja/1.11.1

module list
