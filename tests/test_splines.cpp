#include <iostream>
#include <math.h>
#include <vector>

#include <boundary_conditions.h>

using namespace std;

using real = double ;

//BoundaryCondition::sll_p_periodic;
//constexpr int BC_GREVILLE = BoundaryCondition::sll_p_greville;
//constexpr int BC_HERMITE  = BoundaryCondition::sll_p_hermite;


real eval_cos_scalar(real &x, vector<real> &coeffs, int *derivative=NULL) 
{
    real y;

    int d;
    constexpr real pi = 3.141592653589793;
    constexpr real pi_2 = 6.283185307179586;

    if ( derivative==NULL )
        d = 0 ;
    else
        d = *derivative; 

    y = pow(pi_2 * coeffs[0],d)*cos( 0.5*pi*d + pi_2 * (coeffs[0]*x + coeffs[1]) ) ;

    cout << "y=" << y << endl;
    return y;
}


void cos_splines_test(
      real &max_norm_error,
      real &max_norm_error_diff,
      real &max_norm_error_int,
      int  const &degree,
      real const &h,
      int  const &N,
      real const &x0,
      real const &xN,
      int  const &bc_xmin,
      int  const &bc_xmax,
      real const *coeffs,
      vector<real> *eval_pts_input)
{
    
}


int main()
{
  //const int ncells;
  //const int bc_xmin, bc_xmax;

  constexpr real x0 = 0.0;
  constexpr real xN = 1.0;
  constexpr int N = 22;

  constexpr real h = (xN-x0) / float(N);

  cout << "h=" << h << endl;

  // coordinate array initialization
  std::vector <real> eval_pts(N);
  //for (int i=0; i<eval_pts.size();i++) eval_pts.operator[](i) = x0+i*h;
  for (int i=0; i<eval_pts.size();i++) eval_pts.at(i) = x0+i*h;
  for (real pts: eval_pts) cout << pts << endl;

  //max_norm_profile = 1.0
  //
  //write(*,*)
  //write(*,'(a)') '---------------------------------------------------------------------------'
  //write(*,'(a)') ' TEST: convergence analysis on cos profile (with absolute error bound)     '
  //write(*,'(a)') '---------------------------------------------------------------------------'
  //write(*,*)
}
