#include <cassert>
#include <cmath>
#include <complex> //This library is declared before fftw3.h
#include <vector>

#include <fftw3.h>
#include <gtest/gtest.h>

#include "fft_1d.h"

using namespace std;

void compute_mesh(vector<double>& mesh1d, double const xmin, double const xmax, int const Nx)
{
    assert(mesh1d.size() == Nx);
    const double h = abs(xmax - xmin) / Nx;
    for (int ii = 0; ii < Nx; ++ii) {
        mesh1d[ii] = xmin + ii * h;
    }
}

double compute_f(const double xpts)
{
    //return cos(8. * M_PI * xpts);
    return cos(xpts);
}

void compute_f(vector<double>& fx, vector<double> mesh1d)
{
    assert(fx.size() == mesh1d.size());
    for (int ii = 0; ii < fx.size(); ++ii) {
        fx[ii] = compute_f(mesh1d[ii]);
    }
}

double compute_deriv_f(const double xpts)
{
    //return -8. * M_PI * sin(8. * M_PI * xpts);
    return -sin(xpts);
}

void compute_deriv_f(vector<double>& dfx_dx, vector<double> mesh1d)
{
    assert(dfx_dx.size() == mesh1d.size());
    for (int ii = 0; ii < dfx_dx.size(); ++ii) {
        dfx_dx[ii] = compute_deriv_f(mesh1d[ii]);
    }
}


/* -----------------------------------------------------------------------------
/> Test the property iFFT(FFT(f))=f for f any function. 
/>  Here f(x) = cos(8*pi*x)
!-----------------------------------------------------------------------------*/
TEST(FFT, fft_forward_backward)
{
    constexpr int Nx = 32;

    // Compute the mesh
    constexpr double xmin = 0.;
    constexpr double xmax = 1.;
    vector<double> mesh1d(Nx);
    compute_mesh(mesh1d, xmin, xmax, Nx);

    // Compute f(x) on the mesh
    vector<double> fx(Nx);
    compute_f(fx, mesh1d);

    // Prepare FFTW plan
    fftw_complex in[Nx];
    fftw_complex fft_in[Nx];
    fftw_complex ifft_fft_in[Nx];
    fftw_plan fft_plan;
    fftw_plan ifft_plan;

    fft_plan = fftw_plan_dft_1d(Nx, in, fft_in, FFTW_FORWARD, FFTW_ESTIMATE);
    ifft_plan = fftw_plan_dft_1d(Nx, fft_in, ifft_fft_in, FFTW_BACKWARD, FFTW_ESTIMATE);

    for (int ii = 0; ii < Nx; ++ii) {
        in[ii][0] = fx[ii]; // real part
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

    constexpr double tol = 1e-15;
    EXPECT_LE(max_error, tol) << "While evaluating iFFT(FFT(f))=f \n"
                              << " -> Tolerance= " << tol << "\n"
                              << " -> Max_error=" << max_error << endl;

    fftw_destroy_plan(fft_plan);
    fftw_destroy_plan(ifft_plan);
}

TEST(FFT, fft1d_get_values)
{
    const int Nx = 31;

    // Compute the mesh
    constexpr double xmin = 0.;
    constexpr double xmax = 1.;
    vector<double> mesh1d(Nx);
    compute_mesh(mesh1d, xmin, xmax, Nx);

    // Compute f(x) on the mesh
    vector<double> fx(Nx);
    compute_f(fx, mesh1d);

    // Prepare the FFTW plan + fill input values
    FFT_1D fft1d(Nx, FFTW_FORWARD);
    fft1d.info();
    DComplexVector values(Nx);
    for (int ii = 0; ii < values.size(); ++ii) {
        values[ii].real(fx[ii]);
        values[ii].imag(0.0);
    }
    fft1d.set_values(values);

    double max_error = 0.;
    for (int ii = 0; ii < values.size(); ++ii) {
        max_error = max(max_error, abs(fft1d.get_value(ii)[0] - values[ii].real()));
    }
    constexpr double tol = 1e-15;
    EXPECT_LE(max_error, tol) << "Error in the copy of the values given by get_value \n"
                              << " -> Tolerance= " << tol << "\n"
                              << " -> Max_error= " << max_error << endl;

    max_error = 0.;
    fftw_complex* in;
    in = &fft1d.get_values();

    for (int ii = 0; ii < values.size(); ++ii) {
        max_error = max(max_error, abs(in[ii][0] - values[ii].real()));
    }
    EXPECT_LE(max_error, tol) << "Error in the copy of the values obtained with get_values"
                              << " -> Tolerance= " << tol << "\n"
                              << " -> Max_error= " << max_error << endl;
}


TEST(FFT, fft1d_forward_backward)
{
    const int Nx = 32;

    // Compute the mesh
    constexpr double xmin = 0.;
    constexpr double xmax = 1.;
    vector<double> mesh1d(Nx);
    compute_mesh(mesh1d, xmin, xmax, Nx);

    // Compute f(x) on the mesh
    vector<double> fx(Nx);
    compute_f(fx, mesh1d);

    // Prepare the FFTW plan for Forward FFT + fill input values
    FFT_1D fft1d(Nx, FFTW_FORWARD);
    fft1d.info();
    DComplexVector complx_values(Nx);
    for (int ii = 0; ii < complx_values.size(); ++ii) {
        complx_values[ii].real(fx[ii]);
        complx_values[ii].imag(0.0);
    }
    fft1d.set_values(complx_values);

    // Computation of FFT(f)
    fftw_complex* fft_fx;
    fft_fx = &fft1d.get_fourier_values();

    // Prepare the FFTW plan for Backward FFT + fill input values
    FFT_1D inv_fft1d(Nx, FFTW_BACKWARD);
    for (int ii = 0; ii < complx_values.size(); ++ii) {
        complx_values[ii].real(fft_fx[ii][0]);
        complx_values[ii].imag(fft_fx[ii][1]);
    }
    inv_fft1d.set_values(complx_values);

    // Computation of iFFT(FFT(f))
    fftw_complex* ifft_fft_fx;
    constexpr int NORMALISED = 1;
    ifft_fft_fx = &inv_fft1d.get_fourier_values(NORMALISED);

    // Check if iFFT(FFT(fx)) = fx
    constexpr double tol = 1e-15;
    double max_error = 0.;
    for (int ii = 0; ii < complx_values.size(); ++ii) {
        max_error = max(max_error, abs(ifft_fft_fx[ii][0] - fft1d.get_value(ii)[0]));
    }
    EXPECT_LE(max_error, tol) << "While evaluating iFFT(FFT(f))=f"
                              << " -> Tolerance= " << tol << "\n"
                              << " -> Max_error= " << max_error << endl;
}


TEST(FFT, fft1d_deriv)
{
    const int Nx = 32;

    // Compute the mesh
    constexpr double xmin = 0.;
    const double xmax = 2. * M_PI;
    vector<double> mesh1d(Nx);
    compute_mesh(mesh1d, xmin, xmax, Nx);

    // Compute f(x) on the mesh
    vector<double> fx(Nx);
    compute_f(fx, mesh1d);

    // Compute df/dx(x) on the mesh
    vector<double> dfx_dx(Nx);
    compute_deriv_f(dfx_dx, mesh1d);

    // Prepare the FFTW plan for Forward FFT + fill input values
    FFT_1D fft1d(Nx, FFTW_FORWARD);
    fft1d.info();
    DComplexVector complx_values(Nx);
    for (int ii = 0; ii < complx_values.size(); ++ii) {
        complx_values[ii].real(fx[ii]);
        complx_values[ii].imag(0.0);
    }
    fft1d.set_values(complx_values);

    // Computation of FFT(f)
    fftw_complex* fft_fx;
    fft_fx = &fft1d.get_fourier_values();

    // Computation of dfx/dx = iFFT(i kx FFT(f))
    FFT_1D inv_fft1d(Nx, FFTW_BACKWARD);
    double kx_mode;
    constexpr complex I(0., 1.);
    double dx = mesh1d[1]-mesh1d[0];
    for (int ii = 0; ii < complx_values.size(); ++ii) {
        kx_mode = inv_fft1d.get_kx_mode(ii,dx);
        complx_values[ii].real(fft_fx[ii][0]);
        complx_values[ii].imag(fft_fx[ii][1]);
        complx_values[ii] = I * kx_mode * complx_values[ii];
    }
    inv_fft1d.set_values(complx_values);
    fftw_complex* ifft_ikx_fx;
    constexpr int NORMALISED = 1;
    ifft_ikx_fx = &inv_fft1d.get_fourier_values(NORMALISED);

    vector<double> dfx_dx_byfft(Nx);
    for (int ii = 0; ii < complx_values.size(); ++ii) {
        dfx_dx_byfft[ii] = ifft_ikx_fx[ii][0];
    }

    // Check if iFFT(FFT(fx)) = fx
    constexpr double tol = 1e-14;
    double max_error = 0.;
    for (int ii = 0; ii < dfx_dx.size(); ++ii) {
        max_error = max(max_error, abs(dfx_dx[ii] - dfx_dx_byfft[ii]));
    }
    EXPECT_LE(max_error, tol) << "While evaluating dfx_dx by FFT \n"
                              << " -> Tolerance= " << tol << "\n"
                              << " -> Max_error= " << max_error << endl;
}
