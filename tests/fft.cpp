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
    FftwFourierTransform<RDimX> fft;

    constexpr std::size_t N = 10;

    IDimX mesh_x(CoordX(0.), CoordX(20. / N));
    DiscreteDomain<IDimX> domx(mesh_x, IndexX(0), IVectX(N));
    auto meshfx = fft.compute_fourier_domain(domx);

    //   f = [0, 1, ...,   n/2, -n/2+1, ..., -1] / (d*n)   if n is even
    std::array<double, N> expected_freqs {0., 1., 2., 3., 4., 5., -4., -3., -2., -1.};
    for (auto& f : expected_freqs) {
        f /= domx.size() * domx.mesh<IDimX>().step();
    }

    constexpr double tol = 2.e-16;
    for (std::size_t i = 0; i < meshfx.size(); ++i) {
        EXPECT_LT(std::fabs(meshfx.to_real(IndexFx(i)) - expected_freqs[i]), tol);
    }
}


// Check the computation of the wave vector when the number of points is odd
TEST(FFT, DomainOdd)
{
    FftwFourierTransform<RDimX> fft;

    constexpr std::size_t N = 9;

    IDimX mesh_x(CoordX(0.), CoordX(20. / N));
    DiscreteDomain<IDimX> domx(mesh_x, IndexX(0), IVectX(N));
    auto meshfx = fft.compute_fourier_domain(domx);

    //   f = [0, 1, ..., (n-1)/2, -n/2, ..., -1] / (d*n)   if n is odd
    std::array<double, N> expected_freqs {0., 1., 2., 3., 4., -4., -3., -2., -1.};
    for (auto& f : expected_freqs) {
        f /= domx.size() * domx.mesh<IDimX>().step();
    }

    constexpr double tol = 2.e-16;
    for (std::size_t i = 0; i < meshfx.size(); ++i) {
        EXPECT_LT(std::fabs(meshfx.to_real(IndexFx(i)) - expected_freqs[i]), tol);
    }
}


// Check if IFFT(FFT(f(x)))=f(x) by using real-to-complex FFT
//  and complex-to-real inverse FFT
TEST(FFT, IdentityUsingReal)
{
    // Construct a 1D mesh
    constexpr std::size_t N = 32;
    IDimX mesh_x(CoordX(0.), CoordX((2. * M_PI) / N));
    DiscreteDomain<IDimX> domx(mesh_x, IndexX(0), IVectX(N));

    // Compute f(x) on the mesh mesh_x
    FieldX<double> values(domx);
    for (auto i : domx) {
        values(i) = compute_f(domx.to_real(i));
    }

    // Compute FFT(f(x))
    FftwFourierTransform<RDimX> fft;
    auto meshfx = fft.compute_fourier_domain(domx);
    DiscreteDomain<typeof(meshfx)> domfx(meshfx, IVectFx(meshfx.size()));
    Chunk<std::complex<double>, IDomainFx> fft_values(domfx);
    fft(fft_values, values);

    // Compute IFFT(FFT(f(x))) and check if it is equivalent to f(x)
    FftwInverseFourierTransform<RDimX> ifft;
    FieldX<double> ifft_fft_values(domx);
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
TEST(FFT, IdentityUsingComplex)
{
    // Construct a 1D mesh
    constexpr std::size_t N = 32;
    IDimX mesh_x(CoordX(0.), CoordX((2. * M_PI) / N));
    DiscreteDomain<IDimX> domx(mesh_x, IndexX(0), IVectX(N));

    // Compute f(x) on the mesh mesh_x
    FieldX<std::complex<double>> values(domx);
    std::cout << domx.size() << std::endl;
    for (auto i : domx) {
        values(i) = compute_f(domx.to_real(i));
    }

    // Compute FFT(f(x))
    FftwFourierTransform<RDimX> fft;
    auto meshfx = fft.compute_fourier_domain(domx);
    DiscreteDomain<typeof(meshfx)> domfx(meshfx, IVectFx(meshfx.size()));
    Chunk<std::complex<double>, IDomainFx> fft_values(domfx);
    fft(fft_values, values);

    // Compute IFFT(FFT(f(x))) and check if it is equivalent to f(x)
    FieldX<std::complex<double>> ifft_fft_values(domx);
    FftwInverseFourierTransform<RDimX> ifft;
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
    IDimX mesh_x(CoordX(0.), CoordX((2. * M_PI) / N));
    DiscreteDomain<IDimX> domx(mesh_x, IndexX(0), IVectX(N));

    // Compute f(x) on the mesh mesh_x
    FieldX<double> values(domx);
    for (auto i : domx) {
        values(i) = compute_f(domx.to_real(i));
    }

    // Compute df(x)/dx on the mesh mesh_x
    FieldX<double> firstderiv_values(domx);
    for (auto i : domx) {
        firstderiv_values(i) = compute_deriv_f(domx.to_real(i));
    }

    // Compute FFT(f(x))
    FftwFourierTransform<RDimX> fft;
    auto meshfx = fft.compute_fourier_domain(domx);
    DiscreteDomain<typeof(meshfx)> domfx(meshfx, IVectFx(meshfx.size()));
    Chunk<std::complex<double>, IDomainFx> fft_values(domfx);
    fft(fft_values, values);

    // Compute i*kx*FFT(f(x))
    Chunk<std::complex<double>, IDomainFx> ikx_fft_values(domfx);
    ikx_fft_values(domfx.front()) = 0.;
    for (auto it_freq = domfx.cbegin() + 1; it_freq != domfx.cend(); ++it_freq) {
        double const kx = 2. * M_PI * meshfx.to_real(*it_freq);
        ikx_fft_values(*it_freq) = 1.0i * kx * fft_values(*it_freq);
    }

    // Compute IFFT(i*kx*FFT(f(x))) and check if it is equivalent to df/dx
    FftwInverseFourierTransform<RDimX> ifft;
    FieldX<double> firstderiv_values_computed(domx);
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
    FftwFourierTransform<RDimX> fft;

    constexpr double T = 5.3;
    constexpr double f = 1.0 / T;
    constexpr std::size_t M = 2;
    constexpr std::size_t N = 50;

    IDimX mesh_x(CoordX(0.), CoordX(M * T / N));
    DiscreteDomain<IDimX> domx(mesh_x, IndexX(0), IVectX(N));
    IDimFx meshfx = fft.compute_fourier_domain(domx);
    DiscreteDomain<IDimFx> domfx(meshfx, IVectFx(meshfx.size()));

    FieldX<std::complex<double>> values(domx);
    for (auto i : domx) {
        values(i) = std::cos(2. * M_PI * f * domx.to_real(i));
    }

    Chunk<std::complex<double>, IDomainFx> fft_values(domfx);
    fft(fft_values, values);

    constexpr double tol = 1.0e-13;
    for (std::size_t i = 0; i < domfx.size(); ++i) {
        if ((i != M) && (i != domfx.size() - M)) {
            EXPECT_LE(std::fabs(fft_values(IndexFx(i))), tol);
        }
    }
}
