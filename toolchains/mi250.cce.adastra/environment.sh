module purge

module load cpe/23.12
module load craype-accel-amd-gfx90a craype-x86-trento
module load PrgEnv-cray

module load amd-mixed/5.7.1
module load rocm/5.7.1
module load cray-fftw
module load cray-hdf5-parallel
module load cray-python

pip3 install --user --upgrade ninja cmake

module list
