#include <cassert>
#include <memory>
#include <vector>

#include "deprecated/bsplines_non_uniform.h"

namespace stdex = std::experimental;

namespace deprecated {

NonUniformBSplines::NonUniformBSplines(int degree, bool periodic, const std::vector<double>& breaks)
    : BSplines(
            degree,
            periodic,
            breaks.size() - 1, // ncells
            periodic ? (breaks.size() - 1) : (breaks.size() - 1 + degree), // nbasis
            breaks.front(), // xmin
            breaks.back(), // xmax
            false) // radial
{
    assert(degree > 0);
    assert(m_ncells > 0);
    assert(m_xmin < m_xmax);
    assert(breaks.size() == m_ncells + 1);

    m_npoints = breaks.size();
    m_knots = std::make_unique<double[]>(m_npoints + 2 * degree);

    for (int i(0); i < m_npoints; ++i) {
        get_knot(i) = breaks[i];
    }

    // Fill out the extra nodes
    if (periodic) {
        double period = breaks[m_npoints - 1] - breaks[0];
        for (int i(1); i < degree + 1; ++i) {
            get_knot(-i) = breaks[m_npoints - 1 - i] - period;
            get_knot(m_npoints - 1 + i) = breaks[i] + period;
        }
    } else // open
    {
        for (int i(1); i < degree + 1; ++i) {
            get_knot(-i) = breaks[0];
            get_knot(m_npoints - 1 + i) = breaks[m_npoints - 1];
        }
    }
}

void NonUniformBSplines::eval_basis(double x, DSpan1D& values, int& jmin) const
{
    std::vector<double> left(m_degree);
    std::vector<double> right(m_degree);

    assert(x >= m_xmin);
    assert(x <= m_xmax);
    assert(values.extent(0) == m_degree + 1);

    // 1. Compute cell index 'icell'
    int icell = find_cell(x);

    assert(icell >= 0);
    assert(icell <= m_ncells - 1);
    assert(get_knot(icell) <= x);
    assert(get_knot(icell + 1) >= x);

    // 2. Compute index range of B-splines with support over cell 'icell'
    jmin = icell;

    // 3. Compute values of aforementioned B-splines
    double temp;
    values(0) = 1.0;
    for (int j(0); j < m_degree; ++j) {
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

void NonUniformBSplines::eval_deriv(double x, DSpan1D& derivs, int& jmin) const
{
    std::vector<double> left(m_degree);
    std::vector<double> right(m_degree);

    assert(x >= m_xmin);
    assert(x <= m_xmax);
    assert(derivs.extent(0) == m_degree + 1);

    // 1. Compute cell index 'icell'
    int icell = find_cell(x);

    assert(icell >= 0);
    assert(icell <= m_ncells - 1);
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
    for (int j(0); j < m_degree - 1; ++j) {
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
    saved = m_degree * derivs(0) / (get_knot(icell + 1) - get_knot(icell + 1 - m_degree));
    derivs(0) = -saved;
    for (int j(1); j < m_degree; ++j) {
        temp = saved;
        saved = m_degree * derivs(j)
                / (get_knot(icell + j + 1) - get_knot(icell + j + 1 - m_degree));
        derivs(j) = temp - saved;
    }
    derivs(m_degree) = saved;
}
void NonUniformBSplines::eval_basis_and_n_derivs(double x, int n, DSpan2D& derivs, int& jmin) const
{
    std::vector<double> left(m_degree);
    std::vector<double> right(m_degree);

    std::vector<double> a_ptr(2 * (m_degree + 1));
    stdex::mdspan<double, stdex::dynamic_extent, 2> a(a_ptr.data(), m_degree + 1);

    std::vector<double> ndu_ptr((m_degree + 1) * (m_degree + 1));
    DSpan2D ndu(ndu_ptr.data(), m_degree + 1, m_degree + 1);

    assert(x >= m_xmin);
    assert(x <= m_xmax);
    assert(n >= 0);
    assert(n <= m_degree);
    assert(derivs.extent(0) == 1 + m_degree);
    assert(derivs.extent(1) == 1 + n);

    // 1. Compute cell index 'icell' and x_offset
    int icell(find_cell(x));

    assert(icell >= 0);
    assert(icell <= m_ncells - 1);
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
    for (int j(0); j < m_degree; ++j) {
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
    for (int j(0); j < m_degree + 1; ++j) {
        derivs(j, 0) = ndu(m_degree, j);
    }

    for (int r(0); r < m_degree + 1; ++r) {
        int s1 = 0;
        int s2 = 1;
        a(0, 0) = 1.0;
        for (int k(1); k < n + 1; ++k) {
            double d(0.0);
            int rk = r - k;
            int pk = m_degree - k;
            if (r >= k) {
                a(0, s2) = a(0, s1) * ndu(rk, pk + 1);
                d = a(0, s2) * ndu(pk, rk);
            }
            int j1 = rk > -1 ? 1 : (-rk);
            int j2 = (r - 1) <= pk ? k : (m_degree - r + 1);
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

    int r(m_degree);
    for (int k(1); k < n + 1; ++k) {
        for (int i(0); i < derivs.extent(0); i++) {
            derivs(i, k) *= r;
        }
        r *= (m_degree - k);
    }
}

int NonUniformBSplines::find_cell(double x) const
{
    if (x > m_xmax)
        return -1;
    if (x < m_xmin)
        return -1;

    if (x == m_xmin)
        return 0;
    if (x == m_xmax)
        return m_ncells - 1;

    // Binary search
    int low(0), high(m_ncells);
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

void NonUniformBSplines::integrals(DSpan1D& int_vals) const
{
    assert(int_vals.extent(0) == m_nbasis + m_degree * m_periodic);

    double inv_deg(1.0 / (m_degree + 1));

    for (int i(0); i < m_nbasis; ++i) {
        int_vals(i) = (get_knot(i + 1) - get_knot(i - m_degree)) * inv_deg;
    }

    if (m_periodic) {
        for (int i(0); i < m_degree; ++i) {
            int_vals(m_nbasis + i) = 0;
        }
    }
}

bool NonUniformBSplines::is_uniform() const
{
    return false;
}


} // namespace deprecated
