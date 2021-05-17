#include <cassert>

#include "bsplines_non_uniform.h"

BSplines_non_uniform::BSplines_non_uniform(int degree, bool periodic, std::vector<double> breaks)
    : BSplines(
            degree,
            periodic,
            false, // uniform
            breaks.size() - 1, // ncells
            periodic ? (breaks.size() - 1) : (breaks.size() - 1 + degree), // nbasis
            breaks[0], // xmin
            breaks.back(), // xmax
            false) // radial
{
    assert(degree > 0);
    assert(ncells > 0);
    assert(xmin < xmax);
    assert(breaks.size() == ncells + 1);

    npoints = breaks.size();
    knots = new double[npoints + 2 * degree];

    for (int i(0); i < npoints; ++i) {
        get_knot(i) = breaks[i];
    }

    // Fill out the extra nodes
    if (periodic) {
        double period = breaks[npoints - 1] - breaks[0];
        for (int i(1); i < degree + 1; ++i) {
            get_knot(-i) = breaks[npoints - 1 - i] - period;
            get_knot(npoints - 1 + i) = breaks[i] + period;
        }
    } else // open
    {
        for (int i(1); i < degree + 1; ++i) {
            get_knot(-i) = breaks[0];
            get_knot(npoints - 1 + i) = breaks[npoints - 1];
        }
    }
}

BSplines_non_uniform::~BSplines_non_uniform()
{
    delete[] knots;
}

BSplines* BSplines::new_bsplines(int degree, bool periodic, std::vector<double> const& breaks)
{
    return new BSplines_non_uniform(degree, periodic, breaks);
}

void BSplines_non_uniform::eval_basis(double x, mdspan_1d& values, int& jmin) const
{
    double left[degree];
    double right[degree];

    assert(x >= xmin);
    assert(x <= xmax);
    assert(values.extent(0) == degree + 1);

    // 1. Compute cell index 'icell'
    int icell = find_cell(x);

    assert(icell >= 0);
    assert(icell <= ncells - 1);
    assert(get_knot(icell) <= x);
    assert(get_knot(icell + 1) >= x);

    // 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell;

    // 3. Compute values of aforementioned B-splines
    double temp;
    values(0) = 1.0;
    for (int j(0); j < degree; ++j) {
        left[j] = x - get_knot(icell - j);
        right[j] = get_knot(icell + j + 1) - x;
        double saved(0.0);
        for (int r(0); r < j + 1; ++r) {
            temp = values(r) / (right[r] + left[j - r]);
            values(r) = saved + right[r] * temp;
            saved = left[j - r] * temp;
        }
        values(j + 1) = saved;
    }
}

void BSplines_non_uniform::eval_deriv(double x, mdspan_1d& derivs, int& jmin) const
{
    double left[degree];
    double right[degree];

    assert(x >= xmin);
    assert(x <= xmax);
    assert(derivs.extent(0) == degree + 1);

    // 1. Compute cell index 'icell'
    int icell = find_cell(x);

    assert(icell >= 0);
    assert(icell <= ncells - 1);
    assert(get_knot(icell) <= x);
    assert(get_knot(icell + 1) >= x);

    // 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell;

    // 3. Compute values of aforementioned B-splines

    /*
     * Compute nonzero basis functions and knot differences
     * for splines up to degree degree-1 which are needed to compute derivative
     * First part of Algorithm  A3.2 of NURBS book
     */
    double saved, temp;
    derivs(0) = 1.0;
    for (int j(0); j < degree - 1; ++j) {
        left[j] = x - get_knot(icell - j);
        right[j] = get_knot(icell + j + 1) - x;
        saved = 0.0;
        for (int r(0); r < j + 1; ++r) {
            temp = derivs(r) / (right[r] + left[j - r]);
            derivs(r) = saved + right[r] * temp;
            saved = left[j - r] * temp;
        }
        derivs(j + 1) = saved;
    }

    /*
     * Compute derivatives at x using values stored in bsdx and formula
     * for spline derivative based on difference of splines of degree degree-1
     */
    saved = degree * derivs(0) / (get_knot(icell + 1) - get_knot(icell + 1 - degree));
    derivs(0) = -saved;
    for (int j(1); j < degree; ++j) {
        temp = saved;
        saved = degree * derivs(j) / (get_knot(icell + j + 1) - get_knot(icell + j + 1 - degree));
        derivs(j) = temp - saved;
    }
    derivs(degree) = saved;
}
void BSplines_non_uniform::eval_basis_and_n_derivs(double x, int n, mdspan_2d& derivs, int& jmin)
        const
{
    double left[degree];
    double right[degree];

    double a_ptr[2 * (degree + 1)];
    std::experimental::mdspan<double, std::experimental::dynamic_extent, 2> a(a_ptr, degree + 1);

    double ndu_ptr[(degree + 1) * (degree + 1)];
    mdspan_2d ndu(ndu_ptr, degree + 1, degree + 1);

    assert(x >= xmin);
    assert(x <= xmax);
    assert(n >= 0);
    assert(n <= degree);
    assert(derivs.extent(0) == 1 + degree);
    assert(derivs.extent(1) == 1 + n);

    // 1. Compute cell index 'icell' and x_offset
    int icell(find_cell(x));

    assert(icell >= 0);
    assert(icell <= ncells - 1);
    assert(get_knot(icell) <= x);
    assert(get_knot(icell + 1) >= x);

    // 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell;

    // 3. Compute nonzero basis functions and knot differences for splines
    //    up to degree (degree-1) which are needed to compute derivative
    //    Algorithm  A2.3 of NURBS book
    //
    //    21.08.2017: save inverse of knot differences to avoid unnecessary
    //    divisions
    //                [Yaman Güçlü, Edoardo Zoni]

    double saved, temp;
    ndu(0, 0) = 1.0;
    for (int j(0); j < degree; ++j) {
        left[j] = x - get_knot(icell - j);
        right[j] = get_knot(icell + j + 1) - x;
        saved = 0.0;
        for (int r(0); r < j + 1; ++r) {
            // compute inverse of knot differences and save them into lower
            // triangular part of ndu
            ndu(r, j + 1) = 1.0 / (right[r] + left[j - r]);
            // compute basis functions and save them into upper triangular part
            // of ndu
            temp = ndu(j, r) * ndu(r, j + 1);
            ndu(j + 1, r) = saved + right[r] * temp;
            saved = left[j - r] * temp;
        }
        ndu(j + 1, j + 1) = saved;
    }
    // Save 0-th derivative
    for (int j(0); j < degree + 1; ++j) {
        derivs(j, 0) = ndu(degree, j);
    }

    for (int r(0); r < degree + 1; ++r) {
        int s1 = 0;
        int s2 = 1;
        a(0, 0) = 1.0;
        for (int k(1); k < n + 1; ++k) {
            double d(0.0);
            int rk = r - k;
            int pk = degree - k;
            if (r >= k) {
                a(0, s2) = a(0, s1) * ndu(rk, pk + 1);
                d = a(0, s2) * ndu(pk, rk);
            }
            int j1 = rk > -1 ? 1 : (-rk);
            int j2 = (r - 1) <= pk ? k : (degree - r + 1);
            for (int j(j1); j < j2; ++j) {
                a(j, s2) = (a(j, s1) - a(j - 1, s1)) * ndu(rk + j, pk + 1);
                d += a(j, s2) * ndu(pk, rk + j);
            }
            if (r <= pk) {
                a(k, s2) = -a(k - 1, s1) * ndu(r, pk + 1);
                d += a(k, s2) * ndu(pk, r);
            }
            derivs(r, k) = d;
            int tmp(s1);
            s1 = s2;
            s2 = tmp;
        }
    }

    int r(degree);
    for (int k(1); k < n + 1; ++k) {
        for (int i(0); i < derivs.extent(0); i++) {
            derivs(i, k) *= r;
        }
        r *= (degree - k);
    }
}

int BSplines_non_uniform::find_cell(double x) const
{
    if (x > xmax)
        return -1;
    if (x < xmin)
        return -1;

    if (x == xmin)
        return 0;
    if (x == xmax)
        return ncells - 1;

    // Binary search
    int low(0), high(ncells);
    int icell((low + high) / 2);
    while (x < get_knot(icell) or x >= get_knot(icell + 1)) {
        if (x < get_knot(icell)) {
            high = icell;
        } else {
            low = icell;
        }
        icell = (low + high) / 2;
    }
    return icell;
}

void BSplines_non_uniform::integrals(mdspan_1d& int_vals) const
{
    assert(int_vals.extent(0) == nbasis + degree * periodic);

    double inv_deg(1.0 / (degree + 1));

    for (int i(0); i < nbasis; ++i) {
        int_vals(i) = (get_knot(i + 1) - get_knot(i - degree)) * inv_deg;
    }

    if (periodic) {
        for (int i(0); i < degree; ++i) {
            int_vals(nbasis + i) = 0;
        }
    }
}
