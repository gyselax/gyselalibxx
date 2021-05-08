#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary_conditions.h"
#include "bsplines.h"
#include <spline_interpolator_1d.h>

using namespace std;
using namespace std::experimental;

//BoundaryCondition::sll_p_periodic;
constexpr BoundCond BC_GREVILLE = BoundCond::GREVILLE;
//constexpr int BC_HERMITE  = BoundaryCondition::sll_p_hermite;

constexpr double TWO_PI = 2. * M_PI;

static inline double eval_cos(
        double const x,
        array<double, 2> const& coeffs,
        int const derivative = 0);

vector<double>& eval_cos(
        vector<double>& y,
        vector<double> const& x,
        array<double, 2> const& coeffs,
        int const derivative = 0);


void cos_splines_test(
        double& max_norm_error,
        double& max_norm_error_diff,
        double& max_norm_error_int,
        int const degree,
        double const h,
        int const N,
        double const x0,
        double const xN,
        int const bc_xmin,
        int const bc_xmax,
        array<double, 2> const& coeffs,
        mdspan_1d const& eval_pts_input = {});


static inline double eval_cos(double const x, vector<double> const& coeffs, int const derivative)
{
    return pow(TWO_PI * coeffs[0], derivative)
         * cos(M_PI_2 * derivative + TWO_PI * (coeffs[0] * x + coeffs[1]));
}

vector<double>& eval_cos(
        vector<double>& y,
        vector<double> const& x,
        vector<double> const& coeffs,
        int const derivative)
{
    assert(y.size() == x.size());
    for (int ii = 0; ii < y.size(); ++ii) {
        y[ii] = eval_cos(x[ii], coeffs, derivative);
    }
    return y;
}


void cos_splines_test(
        double& max_norm_error,
        double& max_norm_error_diff,
        double& max_norm_error_int,
        int const degree,
        double const h,
        int const N,
        double const x0,
        double const xN,
        BoundCond const bc_xmin,
        BoundCond const bc_xmax,
        array<double, 2> const coeffs,
        mdspan_1d const& eval_pts_input)
{
    // Create B-splines (uniform or non-uniform depending on input)
    BSplines* bspline = BSplines::new_bsplines(
        degree,
        (bc_xmin == BoundCond::PERIODIC),
        x0,
        xN,
        Spline_interpolator_1D::compute_num_cells(degree, bc_xmin, bc_xmax , N));

      // Initialize 1D spline
      Spline_1D spline( *bspline );

      // Initialize 1D spline interpolator
      Spline_interpolator_1D spline_interpolator( *bspline, bc_xmin, bc_xmax );

       mdspan_1d const xgrid = spline_interpolator.get_interp_points();
       mdspan_1d eval_pts;
       if ( 0 != eval_pts_input.extent(0) ) {
           eval_pts = eval_pts_input;
       } else {
          eval_pts = xgrid;
       }

//TODO:       auto&& yvals = eval_cos( vector<double>(N), xgrid, coeffs )

      
      /* TODO
      !*** computation of rhs'(0) and rhs'(n)  ***
      !-> deriv_rhs(0) = rhs'(n) and deriv_rhs(1) = rhs'(0)
      s = 1-modulo(degree,2) ! shift = 1 for even order, 0 for odd order
      if ( bc_xmin == BC_HERMITE ) then
          do i = 1, degree/2
            Sderiv_lhs(i-1) = eval_cos(x0, coeffs, i-s)
          end do
          if ( bc_xmax == BC_HERMITE ) then
              do i = 1, degree/2
                Sderiv_rhs(i-1) = eval_cos(xN, coeffs, i-s)
              end do

              call spline1d_interpolator__interpolate( spline_interpolator, spline, yvals, bc_xmin, bc_xmax, Sderiv_lhs, Sderiv_rhs)
          else
              call spline1d_interpolator__interpolate( spline_interpolator, spline, yvals, bc_xmin, bc_xmax, Sderiv_lhs)
          endif
      else
          if ( bc_xmax == BC_HERMITE ) then
              do i = 1, degree/2
                Sderiv_rhs(i-1) = eval_cos(xN, coeffs, i-s)
              end do

              call spline1d_interpolator__interpolate( spline_interpolator, spline, yvals, bc_xmin, bc_xmax, derivs_xmax = Sderiv_rhs)
          else
              call spline1d_interpolator__interpolate( spline_interpolator, spline, yvals, bc_xmin, bc_xmax)
          endif
      endif

      max_norm_error = 0.0_F64
      max_norm_error_diff = 0.0_F64
      do i = lbound(eval_pts,1), ubound(eval_pts,1)
          !*** Check eval function ***
          spline_value = spline1d_eval( spline, eval_pts(i) )

          !*** Compute error ***
          error = spline_value - eval_cos(eval_pts(i), coeffs)
          max_norm_error = max( max_norm_error, abs( error ) )


          !*** Check eval_deriv function ***
          spline_value = spline1d_eval_deriv( spline, eval_pts(i) )

          !*** Compute error ***
          error = spline_value - eval_cos(eval_pts(i), coeffs, 1)
          max_norm_error_diff = max( max_norm_error_diff, abs( error ) )
      end do

      max_norm_error_int = abs( spline1d_integrate( spline ) )

      !*** array deallocation ***
      call bspline1d__free(bspline)
      call spline1d_del(spline)
      call spline1d_interpolator__del(spline_interpolator)
      */
}


int main()
{
    //const int ncells;
    //const int bc_xmin, bc_xmax;

    constexpr double x0 = 0.0;
    constexpr double xN = 1.0;
    constexpr int N = 22;

    constexpr double h = (xN - x0) / N;

    cout << "h=" << h << endl;

    // coordinate array initialization
    std::vector<double> eval_pts(N);
    for (int i = 0; i < eval_pts.size(); i++) {
        eval_pts[i] = x0 + i * h;
    }

    for (double pts : eval_pts) {
        cout << pts << endl;
    }

    //max_norm_profile = 1.0
    //
    //write(*,*)
    //write(*,'(a)') '---------------------------------------------------------------------------'
    //write(*,'(a)') ' TEST: convergence analysis on cos profile (with absolute error bound)     '
    //write(*,'(a)') '---------------------------------------------------------------------------'
    //write(*,*)
}

