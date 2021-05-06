#include <array>
#include <cmath>
#include <iostream>
#include <vector>

#include "boundary_conditions.h"

using namespace std;

//BoundaryCondition::sll_p_periodic;
constexpr BoundaryCondition BC_GREVILLE = BoundaryCondition::sll_p_greville;
//constexpr int BC_HERMITE  = BoundaryCondition::sll_p_hermite;


double eval_cos_scalar(double const x, vector<double> const & coeffs, int const derivative=0)
{
    constexpr double TWOPI = 2.*M_PI;

    int const d = derivative;

    double const y = pow(TWOPI * coeffs[0], d) * cos(M_PI_2 * d + M_PI * (coeffs[0] * x + coeffs[1]));

    cout << "y=" << y << endl;
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
        int const bc_xmin,
        int const bc_xmax,
        array<double, 2> const coeffs,
        vector<double> const& eval_pts_input)
{
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
    for (int i = 0; i < eval_pts.size(); i++){
        eval_pts[i] = x0 + i * h;
    }

    for (double pts : eval_pts){
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
