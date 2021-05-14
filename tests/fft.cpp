#include <gtest/gtest.h>

#include <fftw3.h>

using namespace std;

TEST(FFT, fft_ifft)
{
   constexpr int Nx=30;

   fftw_complex *in, *out;
   fftw_plan my_plan;
   
   in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*Nx);
   my_plan = fftw_plan_dft_1d(Nx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
   fftw_execute(my_plan);

   fftw_destroy_plan(my_plan);
   fftw_free(in);
   fftw_free(out);
}
