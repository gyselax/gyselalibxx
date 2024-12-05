module purge
# necessary dependancies for build
module load \
    gcc/11.2.0/gcc-4.8.5 \
    cmake/3.28.3/gcc-11.2.0 \
    cuda/12.2.1/gcc-11.2.0 \
    pdi/1.6.0/gcc-11.2.0 \
    libyaml/0.2.5/gcc-11.2.0 \
    paraconf/1.0.0/gcc-11.2.0 \
    fftw/3.3.10/gcc-11.2.0-openmpi \
    ginkgo/1.8.0/gcc-11.2.0 \
    openblas/0.3.8/gcc-9.2.0 \
    spack/0.21.1/gcc-11.2.0

# additional module for tests and simulations
module load \
    pdiplugin-decl-hdf5/1.6.0/gcc-11.2.0-openmpi \
    pdiplugin-set-value/1.6.0/gcc-11.2.0 \
    pdiplugin-trace/1.6.0/gcc-11.2.0 \
    python/3.9.10/gcc-11.2.0
