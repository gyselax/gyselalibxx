#include <gtest/gtest.h>

#include "fftw.hpp"
#include "geometry.hpp"
#include "ifftw.hpp"

using namespace std;

double compute_f(const double xpts)
{
    return cos(4. * xpts);
}

double compute_deriv_f(const double xpts)
{
    return -4. * sin(4. * xpts);
}

// Check the computation of the wave vector when the number of points is even
TEST(FFT, DomainEven)
{
    FftwFourierTransform<Dim::X> fft;

    constexpr std::size_t N = 10;

    MeshX mesh_x(RCoordX(0.), RCoordX(20. / N));
    ProductMDomain<MeshX> domx(mesh_x, MCoordX(0), MLengthX(N));
    auto meshfx = fft.compute_fourier_domain(domx);

    //   f = [0, 1, ...,   n/2, -n/2+1, ..., -1] / (d*n)   if n is even
    std::array<double, N> expected_freqs {0., 1., 2., 3., 4., 5., -4., -3., -2., -1.};
    for (auto& f : expected_freqs) {
        f /= domx.size() * domx.mesh<MeshX>().step();
    }

    constexpr double tol = 2.e-16;
    for (std::size_t i = 0; i < meshfx.size(); ++i) {
        EXPECT_LT(std::fabs(meshfx.to_real(MCoordFx(i)) - expected_freqs[i]), tol);
    }
}


// Check the computation of the wave vector when the number of points is odd
TEST(FFT, DomainOdd)
{
    FftwFourierTransform<Dim::X> fft;

    constexpr std::size_t N = 9;

    MeshX mesh_x(RCoordX(0.), RCoordX(20. / N));
    ProductMDomain<MeshX> domx(mesh_x, MCoordX(0), MLengthX(N));
    auto meshfx = fft.compute_fourier_domain(domx);

    //   f = [0, 1, ..., (n-1)/2, -n/2, ..., -1] / (d*n)   if n is odd
    std::array<double, N> expected_freqs {0., 1., 2., 3., 4., -4., -3., -2., -1.};
    for (auto& f : expected_freqs) {
        f /= domx.size() * domx.mesh<MeshX>().step();
    }

    constexpr double tol = 2.e-16;
    for (std::size_t i = 0; i < meshfx.size(); ++i) {
        EXPECT_LT(std::fabs(meshfx.to_real(MCoordFx(i)) - expected_freqs[i]), tol);
    }
}


// Check if IFFT(FFT(f(x)))=f(x) by using real-to-complex FFT
//  and complex-to-real inverse FFT
TEST(FFT, Identity_UsingReal)
{
    // Construct a 1D mesh
    constexpr std::size_t N = 32;
    MeshX mesh_x(RCoordX(0.), RCoordX((2. * M_PI) / N));
    ProductMDomain<MeshX> domx(mesh_x, MCoordX(0), MLengthX(N));

    // Compute f(x) on the mesh mesh_x
    BlockX<double> values(domx);
    for (auto i : domx) {
        values(i) = compute_f(domx.to_real(i));
    }

    // Compute FFT(f(x))
    FftwFourierTransform<Dim::X> fft;
    auto meshfx = fft.compute_fourier_domain(domx);
    ProductMDomain<typeof(meshfx)> domfx(meshfx, MLengthFx(meshfx.size()));
    Block<std::complex<double>, MDomainFx> fft_values(domfx);
    fft(fft_values, values);

    // Compute IFFT(FFT(f(x))) and check if it is equivalent to f(x)
    FftwInverseFourierTransform<Dim::X> ifft;
    BlockX<double> ifft_fft_values(domx);
    ifft(ifft_fft_values, fft_values);

    double max_error = 0.;
    for (auto ii : domx) {
        max_error = std::fmax(max_error, std::abs(values(ii) - ifft_fft_values(ii)));
    }

    constexpr double tol = 1e-15;
    EXPECT_LE(max_error, tol) << "While evaluating iFFT(FFT(f))=f \n"
                              << " -> Tolerance= " << tol << "\n"
                              << " -> Max_error=" << max_error << endl;
}


// Check if IFFT(FFT(f(x)))=f(x) by using complex-to-complex FFT
//  and complex-to-complex inverse FFT
TEST(FFT, Identity_UsingComplex)
{
    // Construct a 1D mesh
    constexpr std::size_t N = 32;
    MeshX mesh_x(RCoordX(0.), RCoordX((2. * M_PI) / N));
    ProductMDomain<MeshX> domx(mesh_x, MCoordX(0), MLengthX(N));

    // Compute f(x) on the mesh mesh_x
    BlockX<std::complex<double>> values(domx);
    std::cout << domx.size() << std::endl;
    for (auto i : domx) {
        values(i) = compute_f(domx.to_real(i));
    }

    // Compute FFT(f(x))
    FftwFourierTransform<Dim::X> fft;
    auto meshfx = fft.compute_fourier_domain(domx);
    ProductMDomain<typeof(meshfx)> domfx(meshfx, MLengthFx(meshfx.size()));
    Block<std::complex<double>, MDomainFx> fft_values(domfx);
    fft(fft_values, values);

    // Compute IFFT(FFT(f(x))) and check if it is equivalent to f(x)
    BlockX<std::complex<double>> ifft_fft_values(domx);
    FftwInverseFourierTransform<Dim::X> ifft;
    ifft(ifft_fft_values, fft_values);

    double max_error = 0.;
    for (auto ii : domx) {
        max_error = std::fmax(max_error, std::abs(values(ii) - ifft_fft_values(ii)));
    }

    constexpr double tol = 1e-15;
    EXPECT_LE(max_error, tol) << "While evaluating iFFT(FFT(f))=f \n"
                              << " -> Tolerance= " << tol << "\n"
                              << " -> Max_error=" << max_error << endl;
}


// Compute first derivative of f(x) by using the equality
// df(x)/dx = IFFT(i*kx*FFT(f(x)))
TEST(FFT, FirstDerivative)
{
    // Construct a 1D mesh
    constexpr std::size_t N = 32;
    MeshX mesh_x(RCoordX(0.), RCoordX((2. * M_PI) / N));
    ProductMDomain<MeshX> domx(mesh_x, MCoordX(0), MLengthX(N));

    // Compute f(x) on the mesh mesh_x
    BlockX<double> values(domx);
    for (auto i : domx) {
        values(i) = compute_f(domx.to_real(i));
    }

    // Compute df(x)/dx on the mesh mesh_x
    BlockX<double> firstderiv_values(domx);
    for (auto i : domx) {
        firstderiv_values(i) = compute_deriv_f(domx.to_real(i));
    }

    // Compute FFT(f(x))
    FftwFourierTransform<Dim::X> fft;
    auto meshfx = fft.compute_fourier_domain(domx);
    ProductMDomain<typeof(meshfx)> domfx(meshfx, MLengthFx(meshfx.size()));
    Block<std::complex<double>, MDomainFx> fft_values(domfx);
    fft(fft_values, values);

    // Compute i*kx*FFT(f(x))
    Block<std::complex<double>, MDomainFx> ikx_fft_values(domfx);
    ikx_fft_values(domfx.front()) = 0.;
    for (auto it_freq = domfx.cbegin() + 1; it_freq != domfx.cend(); ++it_freq) {
        double const kx = 2. * M_PI * meshfx.to_real(*it_freq);
        ikx_fft_values(*it_freq) = 1.0i * kx * fft_values(*it_freq);
    }

    // Compute IFFT(i*kx*FFT(f(x))) and check if it is equivalent to df/dx
    FftwInverseFourierTransform<Dim::X> ifft;
    BlockX<double> firstderiv_values_computed(domx);
    ifft(firstderiv_values_computed, ikx_fft_values);

    double max_error = 0.;
    for (auto ii : domx) {
        max_error = std::
                fmax(max_error, std::abs(firstderiv_values(ii) - firstderiv_values_computed(ii)));
    }

    constexpr double tol = 1e-13;
    EXPECT_LE(max_error, tol) << "While evaluating df/dx: \n"
                              << " -> Tolerance= " << tol << "\n"
                              << " -> Max_error=" << max_error << endl;
}


// A pure sinusoidal signal of frequency f sampled on a domain of length 2 * T.
// The frequency spectrum should have only two peaks, the first one being located at the second
// position.
TEST(FFT, Simple)
{
    FftwFourierTransform<Dim::X> fft;

    constexpr double T = 5.3;
    constexpr double f = 1.0 / T;
    constexpr std::size_t M = 2;
    constexpr std::size_t N = 50;

    MeshX mesh_x(RCoordX(0.), RCoordX(M * T / N));
    ProductMDomain<MeshX> domx(mesh_x, MCoordX(0), MLengthX(N));
    MeshFx meshfx = fft.compute_fourier_domain(domx);
    ProductMDomain<MeshFx> domfx(meshfx, MLengthFx(meshfx.size()));

    BlockX<std::complex<double>> values(domx);
    for (auto i : domx) {
        values(i) = std::cos(2. * M_PI * f * domx.to_real(i));
    }

    Block<std::complex<double>, MDomainFx> fft_values(domfx);
    fft(fft_values, values);

    constexpr double tol = 1.0e-13;
    for (std::size_t i = 0; i < domfx.size(); ++i) {
        if ((i != M) && (i != domfx.size() - M)) {
            EXPECT_LE(std::fabs(fft_values(MCoordFx(i))), tol);
        }
    }
}
