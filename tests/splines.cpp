#include <array>
#include <cassert>
#include <cmath>
#include <vector>

#include <gtest/gtest.h>

#include "boundary_conditions.h"
#include "bsplines.h"
#include "spline_interpolator_1d.h"

using namespace std;
using namespace std::experimental;

constexpr double TWO_PI = 2. * M_PI;

static inline double eval_cos(
        double const x,
        View1D<double> const& coeffs,
        int const derivative = 0);

View1D<double>& eval_cos(
        View1D<double>& y,
        View1D<double> const& x,
        View1D<double> const& coeffs,
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
        View1D<double> const& coeffs,
        mdspan_1d const& eval_pts_input = {});


static inline double eval_cos(double const x, View1D<double> const& coeffs, int const derivative)
{
    return pow(TWO_PI * coeffs[0], derivative)
           * cos(M_PI_2 * derivative + TWO_PI * (coeffs[0] * x + coeffs[1]));
}

View1D<double>& eval_cos(
        View1D<double>& y,
        View1D<double> const& x,
        View1D<double> const& coeffs,
        int const derivative)
{
    assert(y.extent(0) == x.extent(0));
    for (int ii = 0; ii < y.extent(0); ++ii) {
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
        View1D<double> const& coeffs,
        mdspan_1d const& eval_pts_input)
{
    // Create B-splines (uniform or non-uniform depending on input)
    BSplines* bspline = BSplines::new_bsplines(
            degree,
            (bc_xmin == BoundCond::PERIODIC),
            x0,
            xN,
            Spline_interpolator_1D::compute_num_cells(degree, bc_xmin, bc_xmax, N));

    // Initialize 1D spline
    Spline_1D spline(*bspline);

    // Initialize 1D spline interpolator
    Spline_interpolator_1D spline_interpolator(*bspline, bc_xmin, bc_xmax);

    mdspan_1d const xgrid = spline_interpolator.get_interp_points();
    mdspan_1d eval_pts;
    if (0 != eval_pts_input.extent(0)) {
        eval_pts = eval_pts_input;
    } else {
        eval_pts = xgrid;
    }

    vector<double> yvals_data(N);
    View1D<double> yvals(yvals_data.data(), yvals_data.size());
    eval_cos(yvals, xgrid, coeffs);

    // computation of rhs'(0) and rhs'(n)
    // -> deriv_rhs(0) = rhs'(n) and deriv_rhs(1) = rhs'(0)
    int const shift = 1 - (degree % 2); // shift = 1 for even order, 0 for odd order

    if (bc_xmin == BoundCond::HERMITE) {
        vector<double> Sderiv_lhs_data(degree / 2);
        View1D<double> Sderiv_lhs(Sderiv_lhs_data.data(), Sderiv_lhs_data.size());
        for (int i = 1; i <= degree / 2; ++i) {
            Sderiv_lhs(i - 1) = eval_cos(x0, coeffs, i - shift);
        }
        if (bc_xmax == BoundCond::HERMITE) {
            vector<double> Sderiv_rhs_data(degree / 2);
            View1D<double> Sderiv_rhs(Sderiv_rhs_data.data(), Sderiv_rhs_data.size());

            for (int i = 1; i <= degree / 2; ++i) {
                Sderiv_rhs(i - 1) = eval_cos(xN, coeffs, i - shift);
            }

            spline_interpolator.compute_interpolant(spline, yvals, &Sderiv_lhs, &Sderiv_rhs);
        } else {
            spline_interpolator.compute_interpolant(spline, yvals, &Sderiv_lhs);
        }
    } else {
        if (bc_xmax == BoundCond::HERMITE) {
            vector<double> Sderiv_rhs_data(degree / 2);
            View1D<double> Sderiv_rhs(Sderiv_rhs_data.data(), Sderiv_rhs_data.size());
            for (int i = 1; i <= degree / 2; ++i) {
                Sderiv_rhs(i - 1) = eval_cos(xN, coeffs, i - shift);
            }
            spline_interpolator.compute_interpolant(spline, yvals, nullptr, &Sderiv_rhs);

        } else {
            spline_interpolator.compute_interpolant(spline, yvals);
        }
    }

    max_norm_error = 0.;
    max_norm_error_diff = 0.;

    for (int i = 0; i < eval_pts.extent(0); ++i) {
        // Check eval function
        double const spline_value = spline.eval(eval_pts[i]);

        // Compute error
        double const error = spline_value - eval_cos(eval_pts[i], coeffs);
        max_norm_error = max(max_norm_error, abs(error));

        // Check eval_deriv function
        double const spline_deriv_value = spline.eval_deriv(eval_pts[i]);

        // Compute error
        double const error_deriv = spline_deriv_value - eval_cos(eval_pts[i], coeffs, 1);
        max_norm_error_diff = max(max_norm_error_diff, abs(error_deriv));
    }
    /*        
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
      */
    /* TODO
      max_norm_error_int = abs( spline1d_integrate( spline ) )

      !*** array deallocation ***
      call bspline1d__free(bspline)
      call spline1d_del(spline)
      call spline1d_interpolator__del(spline_interpolator)
      */
}


TEST(Splines, test)
{
    //const int ncells;
    //const int bc_xmin, bc_xmax;

    constexpr double x0 = 0.0;
    constexpr double xN = 1.0;
    constexpr int N = 22;

    constexpr double h = (xN - x0) / N;

    // coordinate array initialization
    std::vector<double> eval_pts(N);
    for (int i = 0; i < eval_pts.size(); i++) {
        eval_pts[i] = x0 + i * h;
    }

    //View1D<double> myview(eval_pts.data(), eval_pts.size());

    //max_norm_profile = 1.0
}
