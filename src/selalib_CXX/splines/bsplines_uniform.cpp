#include <cassert>
#include "bsplines_uniform.h"
#include "math_tools.h"

BSplines_uniform::BSplines_uniform(int degree, bool periodic, double xmin, double xmax, int ncells)
    : BSplines(degree,
               periodic,
               true,
               ncells,
               periodic ? ncells : ncells+degree,
               xmin,
               xmax,
               false)
{
    assert( degree > 0 );
    assert( ncells > 0 );
    assert( xmin   < xmax );
    inv_dx = ncells / (xmax - xmin);
    dx     = (xmax - xmin) / ncells;
}

BSplines_uniform::~BSplines_uniform() {}

void BSplines_uniform::eval_basis(double x, mdspan_1d& values, int& jmin, int deg) const
{
    assert( values.extent(0) == deg + 1 );

    double offset;
    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(x, jmin, offset);

    // 3. Compute values of aforementioned B-splines
    double xx, temp, saved;
    values(0) = 1.0;
    for (int j(1); j<deg+1; ++j)
    {
        xx = -offset;
        saved = 0.0;
        for (int r(0); r<j; ++r)
        {
            xx += 1;
            temp = values(r) / j;
            values(r) = saved + xx* temp;
            saved     = (j - xx) * temp;
        }
        values(j) = saved;
    }
}

void BSplines_uniform::eval_deriv(double x, mdspan_1d& derivs, int& jmin) const
{
    assert( derivs.extent(0) == degree + 1 );

    double offset;
    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(x, jmin, offset);

    // 3. Compute derivatives of aforementioned B-splines
    //    Derivatives are normalized, hence they should be divided by dx
    double xx, temp, saved;
    derivs(0) = inv_dx;
    for (int j(1); j<degree; ++j)
    {
        xx = -offset;
        saved = 0.0;
        for (int r(0); r < j; ++r)
        {
            xx += 1.0;
            temp = derivs(r) / j;
            derivs(r) = saved + xx * temp;
            saved = (j - xx) * temp;
        }
        derivs(j) = saved;
    }

    // Compute derivatives
    double bjm1 = derivs(0);
    double bj   = bjm1;
    derivs(0) = -bjm1;
    for (int j(1); j < degree; ++j)
    {
        bj = derivs(j);
        derivs(j) = bjm1 - bj;
        bjm1 = bj;
    }
    derivs(degree) = bj;
}

void BSplines_uniform::eval_basis_and_n_derivs(double x, int n, mdspan_2d& derivs, int& jmin) const
{
    double ndu_ptr[(degree + 1) * (degree + 1)];
    mdspan_2d ndu(ndu_ptr, degree + 1, degree + 1);
    double a_ptr[2*(degree + 1)];
    std::experimental::mdspan<double,std::experimental::dynamic_extent,2> a(a_ptr,degree + 1);
    double offset;

    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(x, jmin, offset);

    // 3. Recursively evaluate B-splines (see "sll_s_uniform_BSplines_eval_basis")
    //    up to self%degree, and store them all in the upper-right triangle of ndu
    double xx, temp, saved;
    ndu(0,0) = 1.0;
    for (int j(1); j<degree + 1; ++j)
    {
        xx = -offset;
        saved = 0.0;
        for (int r(0); r<j; ++r)
        {
            xx += 1.0;
            temp = ndu(j-1,r) / j;
            ndu(j,r) = saved + xx * temp;
            saved = (j - xx) * temp;
        }
        ndu(j,j) = saved;
    }
    for (int i(0); i<ndu.extent(1); ++i)
    {
        derivs(i,0) = ndu(degree,i);
    }

    for (int r(0); r < degree + 1; ++r)
    {
        int s1 = 0;
        int s2 = 1;
        a(0,0) = 1.0;
        for (int k(1); k<n+1; ++k)
        {
            double d(0.0);
            int rk = r - k;
            int pk = degree - k;
            if (r >= k )
            {
                a(0, s2) = a(0, s1) / (pk + 1);
                d = a(0, s2) * ndu(pk, rk);
            }
            int j1 = rk > -1 ? 1 : ( -rk);
            int j2 = (r-1) <= pk ? k : (degree - r + 1);
            for (int j(j1); j<j2; ++j)
            {
                a(j, s2) = (a(j, s1) - a(j-1, s1)) / (pk + 1);
                d += a(j, s2) * ndu(pk, rk + j);
            }
            if (r <= pk)
            {
                a(k, s2) = -a(k-1, s1) / (pk + 1);
                d += a(k, s2) * ndu(pk, r);
            }
            derivs(r, k) = d;
            int tmp (s1);
            s1 = s2;
            s2 = tmp;
        }
    }

    // Multiply result by correct factors:
    // degree!/(degree-n)! = degree*(degree-1)*...*(degree-n+1)
    // k-th derivatives are normalized, hence they should be divided by dx^k
    double d = degree * inv_dx;
    for (int k(1); k<n+1; ++k)
    {
        for (int i(0); i<derivs.extent(0); ++i)
        {
            derivs(i,k) *= d;
        }
        d *= (degree - k) * inv_dx;
    }
}

void BSplines_uniform::get_icell_and_offset(double x, int& icell, double& offset) const
{
    assert( x >= xmin );
    assert( x <= xmax );

    if (x == xmin)
    {
        icell = 0;
        offset = 0.0;
    }
    else if (x == xmax)
    {
        icell = ncells-1;
        offset = 1.0;
    }
    else
    {
        offset = (x - xmin) * inv_dx;
        icell  = int(offset);
        offset = offset - icell;

        // When x is very close to xmax, round-off may cause the wrong answer
        // icell=ncells and x_offset=0, which we convert to the case x=xmax:
        if (icell == ncells and offset == 0.0)
        {
            icell = ncells-1;
            offset = 1.0;
        }
    }
}

void BSplines_uniform::integrate(mdspan_1d& int_vals) const
{
    assert(int_vals.extent(0) == nbasis + degree * periodic);
    for (int i(degree); i < nbasis - degree; ++i)
    {
        int_vals(i) = dx;
    }

    if (periodic)
    {
        // Periodic conditions lead to repeat spline coefficients
        for (int i(0); i < degree; ++i)
        {
            int_vals(i) = dx;
            int_vals(nbasis - i - 1) = dx;
            int_vals(nbasis + i) = 0;
        }
    }
    else
    {
        int jmin(0);
        double edge_vals_ptr[degree + 2];
        mdspan_1d edge_vals(edge_vals_ptr, degree + 2);

        eval_basis(xmin, edge_vals, jmin, degree + 1);

        double d_eval = sum(edge_vals);

        for (int i(0); i < degree; ++i)
        {
            double c_eval = sum(edge_vals, 0, degree - i);

            int_vals(i) = dx * (d_eval - c_eval);
            int_vals(nbasis - 1 - i) = int_vals(i);
        }
    }
}
