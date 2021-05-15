#include <cassert>
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

const fftw_plan FFT_1D::get_plan(){
    return this->plan;
}

void FFT_1D::set_values(DComplexVector& values)
{
    assert(values.size() == get_Npoints());
    cout << "set_values Npoints=" << values.size() << endl;
    for (int ii=0; ii < values.size(); ++ii) {
        this->values[ii][0] = real(values[ii]);
        this->values[ii][1] = imag(values[ii]);
    }
}

const fftw_complex &FFT_1D::get_value(const int ival){
    return this->values[ival];
}

fftw_complex &FFT_1D::get_values(){
    return this->values[0];
}

fftw_complex &FFT_1D::get_fourier_values(){
    fftw_execute(this->plan);
    return this->fourier_values[0];
}

