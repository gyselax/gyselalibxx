ENVIRONMENT_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

module purge

module load cpe/24.07
# FIXME:
# craype-accel-amd-gfx90a Error with the cray wrappers clang: error: unsupported option '-fopenmp-targets=' for language mode 'HIP'
module load craype-x86-trento
module load PrgEnv-cray-amd amd-mixed/6.1.2

module load rocm/6.1.2
module load cray-fftw
module load cray-hdf5-parallel
module load cray-python

# FIXME: SPACK_USER_PREFIX path too long
# export SPACK_USER_PREFIX="${ENVIRONMENT_ROOT_DIRECTORY}/gyselalibxx-spack-install-MI250"
export SPACK_USER_PREFIX="${SHAREDWORKDIR}/gyselalibxx-spack-install-MI250"
module use "${SPACK_USER_PREFIX}/modules/tcl/linux-rhel8-zen3"

module load \
    cmake \
    gcc/13.2.mi250/zen3/libyaml \
    gcc/13.2.mi250/zen3/paraconf \
    gcc/13.2.mi250/zen3/pdi \
    gcc/13.2.mi250/zen3/pdiplugin-decl-hdf5 \
    gcc/13.2.mi250/zen3/pdiplugin-set-value \
    gcc/13.2.mi250/zen3/pdiplugin-trace \
    gcc/13.2.mi250/zen3/ginkgo \
    gcc/13.2.mi250/zen3/eigen \
    gcc/13.2.mi250/zen3/ninja

module list
