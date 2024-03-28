#!/bin/bash

TOOLCHAIN_ROOT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd -P)"

set -eu

cd -- "${TOOLCHAIN_ROOT_DIRECTORY}"

source ./environment.sh

cd dependencies

rm -rf InstallDirectory

BUILD_SHARED_LIBS=ON
_CMAKE_INSTALL_PREFIX="$(pwd)/InstallDirectory"
export CMAKE_PREFIX_PATH="${_CMAKE_INSTALL_PREFIX}:${CMAKE_PREFIX_PATH}"

git clone https://github.com/yaml/libyaml -b 0.2.5 || true

rm -rf libyaml/build
cmake -B libyaml/build -S libyaml -DCMAKE_CXX_COMPILER=CC \
    -DCMAKE_INSTALL_PREFIX="${_CMAKE_INSTALL_PREFIX}" \
    -DBUILD_SHARED_LIBS="${BUILD_SHARED_LIBS}"
cmake --build libyaml/build --parallel --target install

git clone https://github.com/pdidev/paraconf.git -b v1.0 || true

rm -rf paraconf/build
cmake -B paraconf/build -S paraconf -DCMAKE_CXX_COMPILER=CC \
    -DCMAKE_INSTALL_PREFIX="${_CMAKE_INSTALL_PREFIX}" \
    -DBUILD_SHARED_LIBS="${BUILD_SHARED_LIBS}" \
    -DUSE_yaml="SYSTEM"
cmake --build paraconf/build --parallel --target install

git clone https://gitlab.maisondelasimulation.fr/pdidev/pdi.git -b 1.6.0 || true

cd pdi
git apply ../pdi.patch || true
cd ..

rm -rf pdi/build
cmake -B pdi/build -S pdi -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn \
    -DCMAKE_INSTALL_PREFIX="${_CMAKE_INSTALL_PREFIX}" \
    -DBUILD_SHARED_LIBS="${BUILD_SHARED_LIBS}" \
    -DUSE_yaml="SYSTEM" \
    -DUSE_HDF5="SYSTEM" \
    -DUSE_paraconf="SYSTEM" \
    -DBUILD_BENCHMARKING="OFF" \
    -DBUILD_DECL_HDF5_PLUGIN="ON" \
    -DBUILD_FORTRAN="OFF" \
    -DBUILD_NETCDF_PARALLEL="OFF"
cmake --build pdi/build --target install --parallel

# It fully integrates in its environment, detects HIP, builds without warnings,
# perfect compared to the softwares above..
git clone https://github.com/ginkgo-project/ginkgo.git -b v1.7.0 || true

rm -rf ginkgo/build
cmake -B ginkgo/build -S ginkgo -GNinja -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn \
    -DCMAKE_INSTALL_PREFIX="${_CMAKE_INSTALL_PREFIX}" \
    -DBUILD_SHARED_LIBS="${BUILD_SHARED_LIBS}" \
    -DGINKGO_BUILD_EXAMPLES="OFF" \
    -DGINKGO_BUILD_MPI="OFF" \
    -DGINKGO_BUILD_REFERENCE="OFF" \
    -DGINKGO_BUILD_BENCHMARKS="OFF" \
    -DGINKGO_BUILD_TESTS="OFF" \
    -DHIP_PATH="${ROCM_PATH}/hip" \
    -DGINKGO_BUILD_HIP="ON" \
    -DGINKGO_HIP_AMDGPU="gfx90a" \
    -DGINKGO_BUILD_OMP="ON" \
    -DCMAKE_EXE_LINKER_FLAGS="-fopenmp"
cmake --build ginkgo/build --target install --parallel

git clone https://gitlab.com/libeigen/eigen.git -b 3.4.0 || true

rm -rf eigen/build
cmake -B eigen/build -S eigen -GNinja -DCMAKE_C_COMPILER=cc -DCMAKE_CXX_COMPILER=CC -DCMAKE_Fortran_COMPILER=ftn \
    -DCMAKE_INSTALL_PREFIX="${_CMAKE_INSTALL_PREFIX}" \
    -DBUILD_SHARED_LIBS="${BUILD_SHARED_LIBS}"
cmake --build eigen/build --target install --parallel
