#pragma once

#include <array>
#include <cassert>
#include <memory>

#include <ddc/DiscreteDomain>
#include <ddc/UniformDiscretization>

#include "sll/bspline.hpp"
#include "sll/math_tools.hpp"
#include "sll/view.hpp"

template <class Tag, std::size_t D>
class UniformBSplines
{
    static_assert(D > 0, "Parameter `D` must be positive");

private:
    template <class T>
    struct InternalTagGenerator;

    /// An internal tag necessary to allocate an internal discretization function.
    /// It must remain internal, for example it shall not be exposed when returning coordinates. Instead use `Tag`
    using internal_tag = InternalTagGenerator<Tag>;

    using mesh_type = UniformDiscretization<internal_tag>;

    using domain_type = DiscreteDomain<mesh_type>;

public:
    using rdim_type = BSpline<Tag>;

    using tag_type = Tag;

    using rcoord_type = Coordinate<UniformBSplines>;

    using mcoord_type = DiscreteCoordinate<UniformBSplines>;

public:
    static constexpr std::size_t rank()
    {
        return 1;
    }

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
    // In the periodic case, it contains twice the periodic point!!!
    DiscreteDomain<mesh_type> m_domain;

public:
    UniformBSplines() = default;

    /** Constructs a BSpline basis with n equidistant knots over \f$[a, b]\f$
     * 
     * @param rmin    the real coordinate of the first knot
     * @param rmax    the real coordinate of the last knot
     * @param n_knots the number of knots
     */
    explicit UniformBSplines(Coordinate<Tag> rmin, Coordinate<Tag> rmax, std::size_t ncells)
        : m_domain(
                DiscreteCoordinate<mesh_type>(0),
                DiscreteVector<mesh_type>(
                        ncells + 1)) // Create a mesh including the eventual periodic point
    {
        init_discretization<mesh_type>(
                Coordinate<internal_tag>(rmin.value()),
                Coordinate<internal_tag>(rmax.value()),
                DiscreteVector<mesh_type>(ncells + 1));
    }

    UniformBSplines(UniformBSplines const& x) = default;

    UniformBSplines(UniformBSplines&& x) = default;

    ~UniformBSplines() = default;

    UniformBSplines& operator=(UniformBSplines const& x) = default;

    UniformBSplines& operator=(UniformBSplines&& x) = default;

    void eval_basis(DSpan1D values, int& jmin, double x) const
    {
        return eval_basis(values, jmin, x, degree());
    }

    void eval_deriv(DSpan1D derivs, int& jmin, double x) const;

    void eval_basis_and_n_derivs(DSpan2D derivs, int& jmin, double x, int n) const;

    DSpan1D integrals(DSpan1D int_vals) const;

    double get_knot(int idx) const noexcept
    {
        return m_domain.rmin() + idx * step<mesh_type>();
    }

    double rmin() const noexcept
    {
        return to_real(m_domain.front());
    }

    double rmax() const noexcept
    {
        return to_real(m_domain.back());
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
        return 1.0 / step<mesh_type>();
    }

    void eval_basis(DSpan1D values, int& jmin, double x, int degree) const;
    void get_icell_and_offset(int& icell, double& offset, double x) const;
};

template <class Tag, std::size_t D>
void UniformBSplines<Tag, D>::eval_basis(
        DSpan1D const values,
        int& jmin,
        double const x,
        int const deg) const
{
    assert(values.extent(0) == deg + 1);

    double offset;
    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(jmin, offset, x);

    // 3. Compute values of aforementioned B-splines
    double xx, temp, saved;
    values(0) = 1.0;
    for (int j = 1; j < deg + 1; ++j) {
        xx = -offset;
        saved = 0.0;
        for (int r = 0; r < j; ++r) {
            xx += 1;
            temp = values(r) / j;
            values(r) = saved + xx * temp;
            saved = (j - xx) * temp;
        }
        values(j) = saved;
    }
}

template <class Tag, std::size_t D>
void UniformBSplines<Tag, D>::eval_deriv(DSpan1D const derivs, int& jmin, double const x) const
{
    assert(derivs.extent(0) == degree() + 1);

    double offset;
    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(jmin, offset, x);

    // 3. Compute derivatives of aforementioned B-splines
    //    Derivatives are normalized, hence they should be divided by dx
    double xx, temp, saved;
    derivs(0) = 1.0 / step<mesh_type>();
    for (int j = 1; j < degree(); ++j) {
        xx = -offset;
        saved = 0.0;
        for (int r = 0; r < j; ++r) {
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
    for (int j = 1; j < degree(); ++j) {
        bj = derivs(j);
        derivs(j) = bjm1 - bj;
        bjm1 = bj;
    }
    derivs(degree()) = bj;
}

template <class Tag, std::size_t D>
void UniformBSplines<Tag, D>::eval_basis_and_n_derivs(
        DSpan2D const derivs,
        int& jmin,
        double const x,
        int const n) const
{
    std::array<double, (degree() + 1) * (degree() + 1)> ndu_ptr;
    std::experimental::mdspan<double, std::experimental::extents<degree() + 1, degree() + 1>> const
            ndu(ndu_ptr.data());
    std::array<double, 2 * (degree() + 1)> a_ptr;
    std::experimental::mdspan<double, std::experimental::extents<degree() + 1, 2>> const a(
            a_ptr.data());
    double offset;

    // 1. Compute cell index 'icell' and x_offset
    // 2. Compute index range of B-splines with support over cell 'icell'
    get_icell_and_offset(jmin, offset, x);

    // 3. Recursively evaluate B-splines (see
    // "sll_s_uniform_BSplines_eval_basis")
    //    up to self%degree, and store them all in the upper-right triangle of
    //    ndu
    double xx, temp, saved;
    ndu(0, 0) = 1.0;
    for (int j = 1; j < degree() + 1; ++j) {
        xx = -offset;
        saved = 0.0;
        for (int r = 0; r < j; ++r) {
            xx += 1.0;
            temp = ndu(j - 1, r) / j;
            ndu(j, r) = saved + xx * temp;
            saved = (j - xx) * temp;
        }
        ndu(j, j) = saved;
    }
    for (int i = 0; i < ndu.extent(1); ++i) {
        derivs(i, 0) = ndu(degree(), i);
    }

    for (int r = 0; r < degree() + 1; ++r) {
        int s1 = 0;
        int s2 = 1;
        a(0, 0) = 1.0;
        for (int k = 1; k < n + 1; ++k) {
            double d = 0.0;
            int const rk = r - k;
            int const pk = degree() - k;
            if (r >= k) {
                a(0, s2) = a(0, s1) / (pk + 1);
                d = a(0, s2) * ndu(pk, rk);
            }
            int const j1 = rk > -1 ? 1 : (-rk);
            int const j2 = (r - 1) <= pk ? k : (degree() - r + 1);
            for (int j = j1; j < j2; ++j) {
                a(j, s2) = (a(j, s1) - a(j - 1, s1)) / (pk + 1);
                d += a(j, s2) * ndu(pk, rk + j);
            }
            if (r <= pk) {
                a(k, s2) = -a(k - 1, s1) / (pk + 1);
                d += a(k, s2) * ndu(pk, r);
            }
            derivs(r, k) = d;
            std::swap(s1, s2);
        }
    }

    // Multiply result by correct factors:
    // degree!/(degree-n)! = degree*(degree-1)*...*(degree-n+1)
    // k-th derivatives are normalized, hence they should be divided by dx^k
    double const inv_dx = inv_step();
    double d = degree() * inv_dx;
    for (int k = 1; k < n + 1; ++k) {
        for (int i = 0; i < derivs.extent(0); ++i) {
            derivs(i, k) *= d;
        }
        d *= (degree() - k) * inv_dx;
    }
}

template <class Tag, std::size_t D>
void UniformBSplines<Tag, D>::get_icell_and_offset(int& icell, double& offset, double const x) const
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
        icell = static_cast<int>(offset);
        offset = offset - icell;

        // When x is very close to xmax, round-off may cause the wrong answer
        // icell=ncells and x_offset=0, which we convert to the case x=xmax:
        if (icell == ncells() && offset == 0.0) {
            icell = ncells() - 1;
            offset = 1.0;
        }
    }
}

template <class Tag, std::size_t D>
DSpan1D UniformBSplines<Tag, D>::integrals(DSpan1D const int_vals) const
{
    assert(int_vals.extent(0) == nbasis() + degree() * is_periodic());
    for (int i = degree(); i < nbasis() - degree(); ++i) {
        int_vals(i) = step<mesh_type>();
    }

    if constexpr (is_periodic()) {
        // Periodic conditions lead to repeat spline coefficients
        for (int i = 0; i < degree(); ++i) {
            int_vals(i) = step<mesh_type>();
            int_vals(nbasis() - i - 1) = step<mesh_type>();
            int_vals(nbasis() + i) = 0;
        }
    } else {
        int jmin = 0;
        std::array<double, degree() + 2> edge_vals_ptr;
        std::experimental::mdspan<double, std::experimental::extents<degree() + 2>> const edge_vals(
                edge_vals_ptr.data());

        eval_basis(edge_vals, jmin, rmin(), degree() + 1);

        double const d_eval = sum(edge_vals);

        for (int i = 0; i < degree(); ++i) {
            double const c_eval = sum(edge_vals, 0, degree() - i);

            int_vals(i) = step<mesh_type>() * (d_eval - c_eval);
            int_vals(nbasis() - 1 - i) = int_vals(i);
        }
    }
    return int_vals;
}
