#include <complex>
#include <vector>

#include <fftw3.h>

typedef std::vector<std::complex<double>> DComplexVector;

class FFT_1D
{
public:
    FFT_1D(const int& Nx, const int& fftw_direction) : Npoints {Nx}, direction {fftw_direction}
    {
        values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx);
        fourier_values = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * Nx);
        plan = fftw_plan_dft_1d(Npoints, values, fourier_values, fftw_direction, FFTW_ESTIMATE);
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
    void set_values(DComplexVector& values);
    const fftw_complex& get_value(const int ival);
    fftw_complex& get_values();
    fftw_complex& get_fourier_values();


protected:
    int Npoints;
    int direction;
    fftw_complex* values;
    fftw_complex* fourier_values;
    fftw_plan plan;
};
