// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include <gtest/gtest.h>

#include "fft.hpp"

struct RDimX;
struct RDimY;
struct RDimZ;

TEST(FFTNorm, OFF)
{
    test_fft_norm<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX>(ddc::FFT_Normalization::OFF);
}

TEST(FFTNorm, BACKWARD)
{
    test_fft_norm<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX>(ddc::FFT_Normalization::BACKWARD);
}

TEST(FFTNorm, FORWARD)
{
    test_fft_norm<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX>(ddc::FFT_Normalization::FORWARD);
}

TEST(FFTNorm, ORTHO)
{
    test_fft_norm<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX>(ddc::FFT_Normalization::ORTHO);
}

TEST(FFTNorm, FULL)
{
    test_fft_norm<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX>(ddc::FFT_Normalization::FULL);
}

TEST(FFTSerialHost, R2C_1D)
{
    test_fft<Kokkos::Serial, Kokkos::Serial::memory_space, float, Kokkos::complex<float>, RDimX>();
}

TEST(FFTSerialHost, R2C_2D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX,
            RDimY>();
}

TEST(FFTSerialHost, R2C_3D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX,
            RDimY,
            RDimZ>();
}

#if fftw_omp_AVAIL
TEST(FFTParallelHost, R2C_1D)
{
    test_fft<Kokkos::OpenMP, Kokkos::OpenMP::memory_space, float, Kokkos::complex<float>, RDimX>();
}

TEST(FFTParallelHost, R2C_2D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX,
            RDimY>();
}

TEST(FFTParallelHost, R2C_3D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX,
            RDimY,
            RDimZ>();
}
#endif

TEST(FFTParallelDevice, R2C_1D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX>();
}

TEST(FFTParallelDevice, R2C_2D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX,
            RDimY>();
}

TEST(FFTParallelDevice, R2C_3D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            float,
            Kokkos::complex<float>,
            RDimX,
            RDimY,
            RDimZ>();
}

TEST(FFTSerialHost, C2C_1D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX>();
}

TEST(FFTSerialHost, C2C_2D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX,
            RDimY>();
}

TEST(FFTSerialHost, C2C_3D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX,
            RDimY,
            RDimZ>();
}

#if fftw_omp_AVAIL
TEST(FFTParallelHost, C2C_1D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX>();
}

TEST(FFTParallelHost, C2C_2D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX,
            RDimY>();
}

TEST(FFTParallelHost, C2C_3D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX,
            RDimY,
            RDimZ>();
}
#endif

TEST(FFTParallelDevice, C2C_1D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX>();
}

TEST(FFTParallelDevice, C2C_2D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX,
            RDimY>();
}

TEST(FFTParallelDevice, C2C_3D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::complex<float>,
            Kokkos::complex<float>,
            RDimX,
            RDimY,
            RDimZ>();
}

TEST(FFTSerialHost, D2Z_1D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX>();
}

TEST(FFTSerialHost, D2Z_2D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX,
            RDimY>();
}

TEST(FFTSerialHost, D2Z_3D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX,
            RDimY,
            RDimZ>();
}

#if fftw_omp_AVAIL
TEST(FFTParallelHost, D2Z_1D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX>();
}

TEST(FFTParallelHost, D2Z_2D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX,
            RDimY>();
}

TEST(FFTParallelHost, D2Z_3D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX,
            RDimY,
            RDimZ>();
}
#endif

TEST(FFTParallelDevice, D2Z_1D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX>();
}

TEST(FFTParallelDevice, D2Z_2D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX,
            RDimY>();
}

TEST(FFTParallelDevice, D2Z_3D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            double,
            Kokkos::complex<double>,
            RDimX,
            RDimY,
            RDimZ>();
}

TEST(FFTSerialHost, Z2Z_1D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX>();
}

TEST(FFTSerialHost, Z2Z_2D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX,
            RDimY>();
}

TEST(FFTSerialHost, Z2Z_3D)
{
    test_fft<
            Kokkos::Serial,
            Kokkos::Serial::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX,
            RDimY,
            RDimZ>();
}

#if fftw_omp_AVAIL
TEST(FFTParallelHost, Z2Z_1D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX>();
}

TEST(FFTParallelHost, Z2Z_2D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX,
            RDimY>();
}

TEST(FFTParallelHost, Z2Z_3D)
{
    test_fft<
            Kokkos::OpenMP,
            Kokkos::OpenMP::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX,
            RDimY,
            RDimZ>();
}
#endif

TEST(FFTParallelDevice, Z2Z_1D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX>();
}

TEST(FFTParallelDevice, Z2Z_2D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX,
            RDimY>();
}

TEST(FFTParallelDevice, Z2Z_3D)
{
    test_fft<
            Kokkos::DefaultExecutionSpace,
            Kokkos::DefaultExecutionSpace::memory_space,
            Kokkos::complex<double>,
            Kokkos::complex<double>,
            RDimX,
            RDimY,
            RDimZ>();
}
