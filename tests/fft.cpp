#include <gtest/gtest.h>

#include<cmath>
#include<complex.h>   //This library is declared before fftw3.h
#include <fftw3.h>
#include<vector>

using namespace std;


TEST(FFT,fft_forward_backward)
// Test the property iFFT(FFT(f))=f for f any function. Here f(x) = cos(8*pi*x)
{
   constexpr int Nx=32;


   fftw_complex in[Nx];
   fftw_complex fft_in[Nx];
   fftw_complex ifft_fft_in[Nx];
   fftw_plan fft_plan;
   fftw_plan ifft_plan;

   fft_plan = fftw_plan_dft_1d(Nx, in, fft_in, FFTW_FORWARD, FFTW_ESTIMATE);
   ifft_plan = fftw_plan_dft_1d(Nx, fft_in, ifft_fft_in, FFTW_BACKWARD, FFTW_ESTIMATE);

   constexpr double x0 = 0.;
   constexpr double h = 1./Nx;
   for (int ii=0;ii<Nx;++ii){
       double xpts = x0 + ii * h;
       in[ii][0] = cos(8.*M_PI*xpts);   // real part
       in[ii][1] = 0.;                  // imaginary part  
   }

   // FFT: Forward Fourier transform
   fftw_execute(fft_plan);

   // iFFT: Backward Fourier transform
   fftw_execute(ifft_plan);
   for (int ii=0;ii<Nx;++ii){
       ifft_fft_in[ii][0] = ifft_fft_in[ii][0]/Nx; 
       ifft_fft_in[ii][1] = ifft_fft_in[ii][1]/Nx; 
   }

   double max_error = 0.;
   for (int ii=0;ii<Nx;++ii){
       max_error = max(max_error,abs(in[ii][0]-ifft_fft_in[ii][0]));
   }
   cout << "max_error=" << max_error << endl;

   constexpr double tol = 1e-14;
   EXPECT_LE(max_error, tol) << "While evaluating iFFT(FFT(f))=f";

   fftw_destroy_plan(fft_plan);
   fftw_destroy_plan(ifft_plan);
}
