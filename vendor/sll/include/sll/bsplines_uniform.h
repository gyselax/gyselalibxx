#pragma once

#include <array>
#include <cassert>
#include <memory>

#include <ddc/MDomain>
#include <ddc/ProductMDomain>
#include <ddc/UniformMesh>

#include "sll/bsplines.h"
#include "sll/math_tools.h"

template <class Tag, std::size_t D>
using UniformBSplines = BSplines<UniformMesh<Tag>, D>;

template <class Tag, std::size_t D>
class BSplines<UniformMesh<Tag>, D>
{
    static_assert(D > 0, "Parameter `D` must be positive");

private:
    using mesh_type = UniformMesh<Tag>;

    using domain_type = ProductMDomain<mesh_type>;

public:
    using tag_type = Tag;

    using rcoord_type = RCoord<Tag>;

    using mcoord_type = MCoord<BSplines>;

public:
    static constexpr std::size_t degree() noexcept
    {
        return D;
    }

    static constexpr bool is_periodic() noexcept
    {
        return Tag::PERIODIC;
    }

    static constexpr bool is_radial() noexcept
    {
        return false;
    }

    static constexpr bool is_uniform() noexcept
    {
        return true;
    }

private:
    mesh_type m_mesh;
    // In the periodic case, it contains twice the periodic point!!!
    MDomain<mesh_type> m_domain;

public:
    BSplines() = default;

    /** Constructs a BSpline basis with n equidistant knots over \f$[a, b]\f$
     * 
     * @param rmin    the real coordinate of the first knot
     * @param rmax    the real coordinate of the last knot
     * @param n_knots the number of knots
     */
    explicit BSplines(rcoord_type rmin, rcoord_type rmax, std::size_t ncells)
        : m_mesh(rmin, rmax, ncells + 1) // Compute correct step, independent of periodicity
        , m_domain(m_mesh, 0, ncells + 1) // Create a mesh including the eventual periodic point
    {
    }

    BSplines(BSplines const& x) = default;

    BSplines(BSplines&& x) = default;

    ~BSplines() = default;

    BSplines& operator=(BSplines const& x) = default;

    BSplines& operator=(BSplines&& x) = default;

    void eval_basis(
            double x,
            std::experimental::mdspan<double, std::experimental::dextents<1>>& values,
            int& jmin) const
    {
        return eval_basis(x, values, jmin, degree());
    }

    void eval_deriv(
            double x,
            std::experimental::mdspan<double, std::experimental::dextents<1>>& derivs,
            int& jmin) const;

    void eval_basis_and_n_derivs(
            double x,
            int n,
            std::experimental::mdspan<double, std::experimental::dextents<2>>& derivs,
            int& jmin) const;

    void integrals(
            std::experimental::mdspan<double, std::experimental::dextents<1>>& int_vals) const;

    double get_knot(int idx) const noexcept
    {
        return m_domain.rmin() + idx * m_domain.mesh().step();
    }

    double rmin() const noexcept
    {
        return m_domain.rmin();
    }

    double rmax() const noexcept
    {
        return m_domain.rmax();
    }

    double length() const noexcept
    {
        return rmax() - rmin();
    }

    std::size_t size() const noexcept
    {
        return degree() + ncells();
    }

    std::size_t nbasis() const noexcept
    {
        return ncells() + !is_periodic() * degree();
    }

    std::size_t ncells() const noexcept
    {
        return m_domain.size() - 1;
    }

private:
    double inv_step() const noexcept
    {
        return 1.0 / m_domain.mesh().step();
    }

    void eval_basis(
            double x,
            std::experimental::mdspan<double, std::experimental::dextents<1>>& values,
            int& jmin,
            int degree) const;
    void get_icell_and_offset(double x, int& icell, double& offset) const;
};

template <class Tag, std::size_t D>
void BSplines<UniformMesh<Tag>, D>::eval_basis(
        double x,
        std::experimental::mdspan<double, std::experimental::dextents<1>>& values,
        int& jmin,
        int deg) const
{
    assert(values.extent(0) == deg + 1);

    double offset;
    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(x, jmin, offset);

    // 3. Compute values of aforementioned B-splines
    double xx, temp, saved;
    values(0) = 1.0;
    for (int j(1); j < deg + 1; ++j) {
        xx = -offset;
        saved = 0.0;
        for (int r(0); r < j; ++r) {
            xx += 1;
            temp = values(r) / j;
            values(r) = saved + xx * temp;
            saved = (j - xx) * temp;
        }
        values(j) = saved;
    }
}

template <class Tag, std::size_t D>
void BSplines<UniformMesh<Tag>, D>::eval_deriv(
        double x,
        std::experimental::mdspan<double, std::experimental::dextents<1>>& derivs,
        int& jmin) const
{
    assert(derivs.extent(0) == degree() + 1);

    double offset;
    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(x, jmin, offset);

    // 3. Compute derivatives of aforementioned B-splines
    //    Derivatives are normalized, hence they should be divided by dx
    double xx, temp, saved;
    derivs(0) = 1.0 / m_domain.mesh().step();
    for (int j(1); j < degree(); ++j) {
        xx = -offset;
        saved = 0.0;
        for (int r(0); r < j; ++r) {
            xx += 1.0;
            temp = derivs(r) / j;
            derivs(r) = saved + xx * temp;
            saved = (j - xx) * temp;
        }
        derivs(j) = saved;
    }

    // Compute derivatives
    double bjm1 = derivs(0);
    double bj = bjm1;
    derivs(0) = -bjm1;
    for (int j(1); j < degree(); ++j) {
        bj = derivs(j);
        derivs(j) = bjm1 - bj;
        bjm1 = bj;
    }
    derivs(degree()) = bj;
}

template <class Tag, std::size_t D>
void BSplines<UniformMesh<Tag>, D>::eval_basis_and_n_derivs(
        double x,
        int n,
        std::experimental::mdspan<double, std::experimental::dextents<2>>& derivs,
        int& jmin) const
{
    std::array<double, (degree() + 1) * (degree() + 1)> ndu_ptr;
    std::experimental::mdspan<double, std::experimental::dextents<2>>
            ndu(ndu_ptr.data(), degree() + 1, degree() + 1);
    std::array<double, 2 * (degree() + 1)> a_ptr;
    std::experimental::mdspan<double, std::experimental::extents<degree() + 1, 2>> a(a_ptr.data());
    double offset;

    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(x, jmin, offset);

    // 3. Recursively evaluate B-splines (see
    // "sll_s_uniform_BSplines_eval_basis")
    //    up to self%degree, and store them all in the upper-right triangle of
    //    ndu
    double xx, temp, saved;
    ndu(0, 0) = 1.0;
    for (int j(1); j < degree() + 1; ++j) {
        xx = -offset;
        saved = 0.0;
        for (int r(0); r < j; ++r) {
            xx += 1.0;
            temp = ndu(j - 1, r) / j;
            ndu(j, r) = saved + xx * temp;
            saved = (j - xx) * temp;
        }
        ndu(j, j) = saved;
    }
    for (int i(0); i < ndu.extent(1); ++i) {
        derivs(i, 0) = ndu(degree(), i);
    }

    for (int r(0); r < degree() + 1; ++r) {
        int s1 = 0;
        int s2 = 1;
        a(0, 0) = 1.0;
        for (int k(1); k < n + 1; ++k) {
            double d(0.0);
            int rk = r - k;
            int pk = degree() - k;
            if (r >= k) {
                a(0, s2) = a(0, s1) / (pk + 1);
                d = a(0, s2) * ndu(pk, rk);
            }
            int j1 = rk > -1 ? 1 : (-rk);
            int j2 = (r - 1) <= pk ? k : (degree() - r + 1);
            for (int j(j1); j < j2; ++j) {
                a(j, s2) = (a(j, s1) - a(j - 1, s1)) / (pk + 1);
                d += a(j, s2) * ndu(pk, rk + j);
            }
            if (r <= pk) {
                a(k, s2) = -a(k - 1, s1) / (pk + 1);
                d += a(k, s2) * ndu(pk, r);
            }
            derivs(r, k) = d;
            int tmp(s1);
            s1 = s2;
            s2 = tmp;
        }
    }

    // Multiply result by correct factors:
    // degree!/(degree-n)! = degree*(degree-1)*...*(degree-n+1)
    // k-th derivatives are normalized, hence they should be divided by dx^k
    double const inv_dx = inv_step();
    double d = degree() * inv_dx;
    for (int k(1); k < n + 1; ++k) {
        for (int i(0); i < derivs.extent(0); ++i) {
            derivs(i, k) *= d;
        }
        d *= (degree() - k) * inv_dx;
    }
}

template <class Tag, std::size_t D>
void BSplines<UniformMesh<Tag>, D>::get_icell_and_offset(double x, int& icell, double& offset) const
{
    assert(x >= rmin());
    assert(x <= rmax());

    double const inv_dx = inv_step();
    if (x == rmin()) {
        icell = 0;
        offset = 0.0;
    } else if (x == rmax()) {
        icell = ncells() - 1;
        offset = 1.0;
    } else {
        offset = (x - rmin()) * inv_dx;
        icell = int(offset);
        offset = offset - icell;

        // When x is very close to xmax, round-off may cause the wrong answer
        // icell=ncells and x_offset=0, which we convert to the case x=xmax:
        if (icell == ncells() and offset == 0.0) {
            icell = ncells() - 1;
            offset = 1.0;
        }
    }
}

template <class Tag, std::size_t D>
void BSplines<UniformMesh<Tag>, D>::integrals(
        std::experimental::mdspan<double, std::experimental::dextents<1>>& int_vals) const
{
    assert(int_vals.extent(0) == nbasis() + degree() * is_periodic());
    for (int i(degree()); i < nbasis() - degree(); ++i) {
        int_vals(i) = m_domain.mesh().step();
    }

    if constexpr (is_periodic()) {
        // Periodic conditions lead to repeat spline coefficients
        for (int i(0); i < degree(); ++i) {
            int_vals(i) = m_domain.mesh().step();
            int_vals(nbasis() - i - 1) = m_domain.mesh().step();
            int_vals(nbasis() + i) = 0;
        }
    } else {
        int jmin(0);
        std::array<double, degree() + 2> edge_vals_ptr;
        std::experimental::mdspan<double, std::experimental::dextents<1>>
                edge_vals(edge_vals_ptr.data(), degree() + 2);

        eval_basis(rmin(), edge_vals, jmin, degree() + 1);

        double d_eval = sum(edge_vals);

        for (int i(0); i < degree(); ++i) {
            double c_eval = sum(edge_vals, 0, degree() - i);

            int_vals(i) = m_domain.mesh().step() * (d_eval - c_eval);
            int_vals(nbasis() - 1 - i) = int_vals(i);
        }
    }
}
