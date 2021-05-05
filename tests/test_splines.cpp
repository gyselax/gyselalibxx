#include <iostream>
using namespace std;

using real = double ;

real eval_cos_scalar(real &x, coeffs, derivative) 
{
  // real(F64)                , intent(in   ) :: x
  //vector<real> coeffs[2];  //real(F64), dimension(0:1), intent(in   ) :: coeffs
  //integer  , optional      , intent(in   ) :: derivative

  int d;                                       //integer :: d
  constexpr real pi = 3.141592653589793;       //real(F64), parameter :: pi = 3.141592653589793_f64
  constexpr real pi_2 = 6.283185307179586_f64; //real(F64), parameter :: pi_2 = 6.283185307179586_f64

  //if (present(derivative)) then; d = derivative; else; d = 0; end if

  y = (pi_2 * coeffs(0))**d * cos( 0.5_F64*pi*d + pi_2 * (coeffs(0)*x + coeffs(1) ) );
  return y;
}


int main()
{

  //real(F64)                 :: h, x0, xN, ncells
  //integer                   :: N, bc_xmin, bc_xmax
  //
  //x0 = 0.0_F64
  //xN = 1.0_F64
  //N  = 22
  //
  //!*** coordinate array allocation ***
  //allocate(eval_pts(0:N))
  //
  //!*** coordinate array initialisation ***
  //h = (xN-x0) / real(N, kind=F64)
  //do i = 0,N
  //  eval_pts(i) = x0+i*h
  //enddo
  //
  //max_norm_profile = 1.0_F64
  //
  //write(*,*)
  //write(*,'(a)') '---------------------------------------------------------------------------'
  //write(*,'(a)') ' TEST: convergence analysis on cos profile (with absolute error bound)     '
  //write(*,'(a)') '---------------------------------------------------------------------------'
  //write(*,*)
}
