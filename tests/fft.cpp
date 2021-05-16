#include <cmath>
#include <complex> //This library is declared before fftw3.h
#include <vector>

#include <fftw3.h>
#include <gtest/gtest.h>

#include "fft_1d.h"

using namespace std;


// Test the property iFFT(FFT(f))=f for f any function. Here f(x) = cos(8*pi*x)
TEST(FFT, fft_forward_backward)
{
    constexpr int Nx = 32;

    fftw_complex in[Nx];
    fftw_complex fft_in[Nx];
    fftw_complex ifft_fft_in[Nx];
    fftw_plan fft_plan;
    fftw_plan ifft_plan;

    fft_plan = fftw_plan_dft_1d(Nx, in, fft_in, FFTW_FORWARD, FFTW_ESTIMATE);
    ifft_plan = fftw_plan_dft_1d(Nx, fft_in, ifft_fft_in, FFTW_BACKWARD, FFTW_ESTIMATE);

    constexpr double x0 = 0.;
    constexpr double h = 1. / Nx;
    for (int ii = 0; ii < Nx; ++ii) {
        double xpts = x0 + ii * h;
        in[ii][0] = cos(8. * M_PI * xpts); // real part
        in[ii][1] = 0.; // imaginary part
    }

    // FFT: Forward Fourier transform
    fftw_execute(fft_plan);

    // iFFT: Backward Fourier transform
    fftw_execute(ifft_plan);
    for (int ii = 0; ii < Nx; ++ii) {
        ifft_fft_in[ii][0] = ifft_fft_in[ii][0] / Nx;
        ifft_fft_in[ii][1] = ifft_fft_in[ii][1] / Nx;
    }

    double max_error = 0.;
    for (int ii = 0; ii < Nx; ++ii) {
        max_error = max(max_error, abs(in[ii][0] - ifft_fft_in[ii][0]));
    }
    cout << "max_error=" << max_error << endl;

    constexpr double tol = 1e-15;
    EXPECT_LE(max_error, tol) << "While evaluating iFFT(FFT(f))=f";

    fftw_destroy_plan(fft_plan);
    fftw_destroy_plan(ifft_plan);
}


TEST(FFT, fft1d)
{
    const int Npoints = 32;
    DComplexVector values(Npoints);
    constexpr double x0 = 0.;
    constexpr double h = 1. / Npoints;
    for (int ii = 0; ii < values.size(); ++ii) {
        double xpts = x0 + ii * h;
        values[ii].real(cos(8. * M_PI * xpts));
        values[ii].imag(0.0);
    }

    // Computation of FFT(f)
    FFT_1D fft1d(Npoints, FFTW_FORWARD);
    fft1d.info();
    fft1d.set_values(values);

    double max_error = 0.;
    for (int ii = 0; ii < values.size(); ++ii) {
        max_error = max(max_error, abs(fft1d.get_value(ii)[0] - values[ii].real()));
    }
    constexpr double tol = 1e-15;
    EXPECT_LE(max_error, tol) << "Error in the copy of the values given by get_value";

    max_error = 0.;
    fftw_complex* in;
    in = &fft1d.get_values();

    for (int ii = 0; ii < values.size(); ++ii) {
        max_error = max(max_error, abs(in[ii][0] - values[ii].real()));
    }
    EXPECT_LE(max_error, tol) << "Error in the copy of the values obtained with get_values";

    // Computation of iFFT(FFT(f))
    fftw_complex* fft_in;
    fft_in = &fft1d.get_fourier_values();
    for (int ii = 0; ii < values.size(); ++ii) {
        values[ii].real(fft_in[ii][0]);
        values[ii].imag(fft_in[ii][1]);
    }

    FFT_1D inv_fft1d(Npoints, FFTW_BACKWARD);
    inv_fft1d.set_values(values);

    fftw_complex* ifft_fft_in;
    ifft_fft_in = &inv_fft1d.get_fourier_values();
    const int Npoints_tmp = inv_fft1d.get_Npoints();
    for (int ii = 0; ii < Npoints_tmp; ++ii) {
        ifft_fft_in[ii][0] = ifft_fft_in[ii][0] / Npoints_tmp;
        ifft_fft_in[ii][1] = ifft_fft_in[ii][1] / Npoints_tmp;
    }

    max_error = 0.;
    for (int ii = 0; ii < values.size(); ++ii) {
        max_error = max(max_error, abs(ifft_fft_in[ii][0] - fft1d.get_value(ii)[0]));
    }
    EXPECT_LE(max_error, tol) << "While evaluating iFFT(FFT(f))=f";
}
