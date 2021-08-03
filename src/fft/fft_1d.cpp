#include <cassert>
#include <cmath>
#include <complex> //This library is declared before fftw3.h
#include <iostream>

#include <fftw3.h>

#include "fft_1d.h"

using namespace std;

void FFT_1D::info()
{
    cout << "Npoints=" << get_Npoints() << endl;
    cout << "Direction=" << this->direction << endl;
}

const int FFT_1D::get_Npoints()
{
    return this->Npoints;
}

const fftw_plan FFT_1D::get_plan()
{
    return this->plan;
}

void FFT_1D::set_values(DComplexVector& values)
{
    assert(values.size() == get_Npoints());
    for (int ii = 0; ii < values.size(); ++ii) {
        this->values[ii][0] = real(values[ii]);
        this->values[ii][1] = imag(values[ii]);
    }
}

const fftw_complex& FFT_1D::get_value(int const ival)
{
    return this->values[ival];
}

fftw_complex& FFT_1D::get_values()
{
    return this->values[0];
}

fftw_complex& FFT_1D::get_fourier_values(int const normalised)
{
    fftw_execute(this->plan);
    if (normalised == 1) {
        const int Nx = get_Npoints();
        const double inv_Nx = 1. / Nx;
        for (int ii = 0; ii < Nx; ++ii) {
            fourier_values[ii][0] = fourier_values[ii][0] * inv_Nx;
            fourier_values[ii][1] = fourier_values[ii][1] * inv_Nx;
        }
    }
    return this->fourier_values[0];
}

const double FFT_1D::get_freq(int const ix)
{
    return this->freqs[ix];
}

const double FFT_1D::get_kx_mode(int const ix, double const dx)
{
    // 2*pi*freq/dx
    return 2 * M_PI * get_freq(ix) / dx;
}
