#include <complex>
#include <iostream>
#include <vector>

#include <fftw3.h>

typedef std::vector<std::complex<double>> DComplexVector;

class FFT_1D
{
public:
    FFT_1D(const int& Nx, const int& fftw_direction)
        : Npoints {Nx}
        , direction {fftw_direction}
        , freqs(Nx)
    {
        values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx);
        fourier_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx);
        plan = fftw_plan_dft_1d(Npoints, values, fourier_values, fftw_direction, FFTW_ESTIMATE);
        set_freqs_();
    };

    ~FFT_1D()
    {
        fftw_destroy_plan(plan);
        fftw_free(values);
        fftw_free(fourier_values);
    };

    void info();
    const int get_Npoints();
    const fftw_plan get_plan();
    void set_values(DComplexVector& in_values);
    const fftw_complex& get_value(int const ival);
    fftw_complex& get_values();
    fftw_complex& get_fourier_values(int const normalised = 0);
    const double get_freq(int const ix);
    const double get_kx_mode(int const ix, const double dx);

private:
    void set_freqs_(const double dspace = 1.)
    /*-----------------------------------------------------------------------------
    !> Return the Discrete Fourier Transform sample frequencies.
    !> 
    !> The returned float array `f` contains the frequency bin centers in cycles
    !> per unit of the sample spacing (with zero at the start).  For instance, if
    !> the sample spacing is in seconds, then the frequency unit is cycles/second.
    !> 
    !> Given a window length `n` and a sample spacing `d`::
    !> 
    !>   f = [0, 1, ...,   n/2, -n/2+1, ..., -1] / (d*n)   if n is even
    !>   f = [0, 1, ..., (n-1)/2, -n/2, ..., -1] / (d*n)   if n is odd
    !> 
    !> Parameters
    !> ----------
    !> n : int
    !>     Window length.
    !> d : scalar, optional
    !>     Sample spacing (inverse of the sampling rate). Defaults to 1.
    !-----------------------------------------------------------------------------*/
    {
        const double inv_Nd = 1. / (Npoints * dspace);
        for (int ii = 0; ii <= Npoints / 2; ++ii) {
            freqs[ii] = (ii)*inv_Nd;
        }
        for (int ii = Npoints / 2 + 1; ii < Npoints; ++ii) {
            freqs[ii] = (ii - Npoints) * inv_Nd;
        }
    }

protected:
    int Npoints;
    int direction;
    fftw_complex* values;
    fftw_complex* fourier_values;
    fftw_plan plan;
    std::vector<double> freqs;
};
