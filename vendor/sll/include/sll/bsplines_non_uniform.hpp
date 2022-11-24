#pragma once

#include <array>
#include <cassert>
#include <memory>
#include <vector>

#include <ddc/ddc.hpp>

#include "sll/bspline.hpp"
#include "sll/view.hpp"

/// NonUniformPointSampling specialization of BSplines
template <class Tag, std::size_t D>
class NonUniformBSplines
{
    static_assert(D > 0, "Parameter `D` must be positive");

public:
    using tag_type = Tag;

    using continuous_dimension_type = BSpline<Tag>;


    using discrete_dimension_type = NonUniformBSplines;

    using discrete_element_type = DiscreteElement<NonUniformBSplines>;

    using discrete_domain_type = DiscreteDomain<NonUniformBSplines>;

    using discrete_vector_type = DiscreteVector<NonUniformBSplines>;

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
        return false;
    }

    template <class MemorySpace>
    class Impl
    {
        template <class OMemorySpace>
        friend class Impl;

    private:
        Kokkos::View<Coordinate<Tag>*, MemorySpace> m_knots;

    public:
        using discrete_dimension_type = NonUniformBSplines;

        Impl() = default;

        template <class OriginMemorySpace>
        explicit Impl(Impl<OriginMemorySpace> const& impl) : m_knots("", impl.m_knots.extent(0))
        {
            Kokkos::deep_copy(m_knots, impl.m_knots);
        }

        /// @brief Construct a `Impl` using a brace-list, i.e. `Impl bsplines({0., 1.})`
        explicit Impl(std::initializer_list<Coordinate<Tag>> breaks)
            : Impl(breaks.begin(), breaks.end())
        {
        }

        /// @brief Construct a `Impl` using a C++20 "common range".
        explicit Impl(std::vector<Coordinate<Tag>> const& breaks)
            : Impl(breaks.begin(), breaks.end())
        {
        }

        /// @brief Construct a `Impl` using a pair of iterators.
        template <class RandomIt>
        Impl(RandomIt breaks_begin, RandomIt breaks_end);

        Impl(Impl const& x) = default;

        Impl(Impl&& x) = default;

        ~Impl() = default;

        Impl& operator=(Impl const& x) = default;

        Impl& operator=(Impl&& x) = default;

        discrete_element_type eval_basis(DSpan1D values, double x) const;

        discrete_element_type eval_deriv(DSpan1D derivs, double x) const;

        discrete_element_type eval_basis_and_n_derivs(DSpan2D derivs, double x, std::size_t n)
                const;

        ChunkSpan<double, discrete_domain_type> integrals(
                ChunkSpan<double, discrete_domain_type> int_vals) const;

        double get_knot(int break_idx) const noexcept
        {
            // TODO: assert break_idx >= 1 - degree
            // TODO: assert break_idx <= npoints + degree
            return m_knots(break_idx + degree());
        }

        double get_first_support_knot(discrete_element_type const& ix) const
        {
            return m_knots(ix.uid());
        }

        double get_last_support_knot(discrete_element_type const& ix) const
        {
            return m_knots(ix.uid() + degree() + 1);
        }

        double get_support_knot_n(discrete_element_type const& ix, int n) const
        {
            return m_knots(ix.uid() + n);
        }

        Coordinate<Tag> rmin() const noexcept
        {
            return Coordinate<Tag>(get_knot(0));
        }

        Coordinate<Tag> rmax() const noexcept
        {
            return Coordinate<Tag>(get_knot(ncells()));
        }

        double length() const noexcept
        {
            return rmax() - rmin();
        }

        std::size_t size() const noexcept
        {
            return degree() + ncells();
        }

        /// Returns the discrete domain including ghost bsplines
        discrete_domain_type full_domain() const
        {
            return discrete_domain_type(discrete_element_type(0), discrete_vector_type(size()));
        }

        std::size_t npoints() const noexcept
        {
            return m_knots.size() - 2 * degree();
        }

        std::size_t nbasis() const noexcept
        {
            return ncells() + !is_periodic() * degree();
        }

        std::size_t ncells() const noexcept
        {
            return npoints() - 1;
        }

    private:
        int find_cell(double x) const;

        double& get_knot(int break_idx)
        {
            // TODO: assert break_idx >= 1 - degree
            // TODO: assert break_idx <= npoints + degree
            return m_knots(break_idx + degree());
        }
    };
};

template <class Tag, std::size_t D>
template <class MemorySpace>
template <class RandomIt>
NonUniformBSplines<Tag, D>::Impl<MemorySpace>::Impl(
        RandomIt const break_begin,
        RandomIt const break_end)
    : m_knots("", (break_end - break_begin) + 2 * degree())
{
    assert(m_knots.size() > 2 * degree() + 1);

    // Fill the provided knots
    int ii = 0;
    for (RandomIt it = break_begin; it < break_end; ++it) {
        get_knot(ii) = *it;
        ++ii;
    }
    assert(rmin() < rmax());

    // Fill out the extra knots
    if constexpr (is_periodic()) {
        double const period = rmax() - rmin();
        for (std::size_t i = 1; i < degree() + 1; ++i) {
            get_knot(-i) = get_knot(ncells() - i) - period;
            get_knot(ncells() + i) = get_knot(i) + period;
        }
    } else // open
    {
        for (std::size_t i = 1; i < degree() + 1; ++i) {
            get_knot(-i) = rmin();
            get_knot(npoints() - 1 + i) = rmax();
        }
    }
}

template <class Tag, std::size_t D>
template <class MemorySpace>
DiscreteElement<NonUniformBSplines<Tag, D>> NonUniformBSplines<Tag, D>::Impl<
        MemorySpace>::eval_basis(DSpan1D const values, double const x) const
{
    std::array<double, degree()> left;
    std::array<double, degree()> right;

    assert(x >= rmin());
    assert(x <= rmax());
    assert(values.extent(0) == degree() + 1);

    // 1. Compute cell index 'icell'
    int const icell = find_cell(x);

    assert(icell >= 0);
    assert(icell <= int(ncells() - 1));
    assert(get_knot(icell) <= x);
    assert(get_knot(icell + 1) >= x);

    // 2. Compute values of B-splines with support over cell 'icell'
    double temp;
    values(0) = 1.0;
    for (std::size_t j = 0; j < degree(); ++j) {
        left[j] = x - get_knot(icell - j);
        right[j] = get_knot(icell + j + 1) - x;
        double saved = 0.0;
        for (std::size_t r = 0; r < j + 1; ++r) {
            temp = values(r) / (right[r] + left[j - r]);
            values(r) = saved + right[r] * temp;
            saved = left[j - r] * temp;
        }
        values(j + 1) = saved;
    }

    return discrete_element_type(icell);
}

template <class Tag, std::size_t D>
template <class MemorySpace>
DiscreteElement<NonUniformBSplines<Tag, D>> NonUniformBSplines<Tag, D>::Impl<
        MemorySpace>::eval_deriv(DSpan1D const derivs, double const x) const
{
    std::array<double, degree()> left;
    std::array<double, degree()> right;

    assert(x >= rmin());
    assert(x <= rmax());
    assert(derivs.extent(0) == degree() + 1);

    // 1. Compute cell index 'icell'
    int const icell = find_cell(x);

    assert(icell >= 0);
    assert(icell <= int(ncells() - 1));
    assert(get_knot(icell) <= x);
    assert(get_knot(icell + 1) >= x);

    // 2. Compute values of derivatives of B-splines with support over cell 'icell'

    /*
     * Compute nonzero basis functions and knot differences
     * for splines up to degree degree-1 which are needed to compute derivative
     * First part of Algorithm  A3.2 of NURBS book
     */
    double saved, temp;
    derivs(0) = 1.0;
    for (std::size_t j = 0; j < degree() - 1; ++j) {
        left[j] = x - get_knot(icell - j);
        right[j] = get_knot(icell + j + 1) - x;
        saved = 0.0;
        for (std::size_t r = 0; r < j + 1; ++r) {
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
    saved = degree() * derivs(0) / (get_knot(icell + 1) - get_knot(icell + 1 - degree()));
    derivs(0) = -saved;
    for (std::size_t j = 1; j < degree(); ++j) {
        temp = saved;
        saved = degree() * derivs(j)
                / (get_knot(icell + j + 1) - get_knot(icell + j + 1 - degree()));
        derivs(j) = temp - saved;
    }
    derivs(degree()) = saved;

    return discrete_element_type(icell);
}

template <class Tag, std::size_t D>
template <class MemorySpace>
DiscreteElement<NonUniformBSplines<Tag, D>> NonUniformBSplines<Tag, D>::Impl<MemorySpace>::
        eval_basis_and_n_derivs(DSpan2D const derivs, double const x, std::size_t const n) const
{
    std::array<double, degree()> left;
    std::array<double, degree()> right;

    std::array<double, 2 * (degree() + 1)> a_ptr;
    std::experimental::mdspan<double, std::experimental::extents<degree() + 1, 2>> const a(
            a_ptr.data());

    std::array<double, (degree() + 1) * (degree() + 1)> ndu_ptr;
    std::experimental::mdspan<double, std::experimental::extents<degree() + 1, degree() + 1>> const
            ndu(ndu_ptr.data());

    assert(x >= rmin());
    assert(x <= rmax());
    assert(n >= 0);
    assert(n <= degree());
    assert(derivs.extent(0) == 1 + degree());
    assert(derivs.extent(1) == 1 + n);

    // 1. Compute cell index 'icell' and x_offset
    int const icell = find_cell(x);

    assert(icell >= 0);
    assert(icell <= int(ncells() - 1));
    assert(get_knot(icell) <= x);
    assert(get_knot(icell + 1) >= x);

    // 2. Compute nonzero basis functions and knot differences for splines
    //    up to degree (degree-1) which are needed to compute derivative
    //    Algorithm  A2.3 of NURBS book
    //
    //    21.08.2017: save inverse of knot differences to avoid unnecessary
    //    divisions
    //                [Yaman Güçlü, Edoardo Zoni]

    double saved, temp;
    ndu(0, 0) = 1.0;
    for (std::size_t j = 0; j < degree(); ++j) {
        left[j] = x - get_knot(icell - j);
        right[j] = get_knot(icell + j + 1) - x;
        saved = 0.0;
        for (std::size_t r = 0; r < j + 1; ++r) {
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
    for (std::size_t j = 0; j < degree() + 1; ++j) {
        derivs(j, 0) = ndu(degree(), j);
    }

    for (int r = 0; r < int(degree() + 1); ++r) {
        int s1 = 0;
        int s2 = 1;
        a(0, 0) = 1.0;
        for (int k = 1; k < int(n + 1); ++k) {
            double d = 0.0;
            int const rk = r - k;
            int const pk = degree() - k;
            if (r >= k) {
                a(0, s2) = a(0, s1) * ndu(rk, pk + 1);
                d = a(0, s2) * ndu(pk, rk);
            }
            int const j1 = rk > -1 ? 1 : (-rk);
            int const j2 = (r - 1) <= pk ? k : (degree() - r + 1);
            for (int j = j1; j < j2; ++j) {
                a(j, s2) = (a(j, s1) - a(j - 1, s1)) * ndu(rk + j, pk + 1);
                d += a(j, s2) * ndu(pk, rk + j);
            }
            if (r <= pk) {
                a(k, s2) = -a(k - 1, s1) * ndu(r, pk + 1);
                d += a(k, s2) * ndu(pk, r);
            }
            derivs(r, k) = d;
            std::swap(s1, s2);
        }
    }

    int r = degree();
    for (int k = 1; k < int(n + 1); ++k) {
        for (std::size_t i = 0; i < derivs.extent(0); i++) {
            derivs(i, k) *= r;
        }
        r *= degree() - k;
    }

    return discrete_element_type(icell);
}

template <class Tag, std::size_t D>
template <class MemorySpace>
int NonUniformBSplines<Tag, D>::Impl<MemorySpace>::find_cell(double const x) const
{
    if (x > rmax())
        return -1;
    if (x < rmin())
        return -1;

    if (x == rmin())
        return 0;
    if (x == rmax())
        return ncells() - 1;

    // Binary search
    int low = 0, high = ncells();
    int icell = (low + high) / 2;
    while (x < get_knot(icell) || x >= get_knot(icell + 1)) {
        if (x < get_knot(icell)) {
            high = icell;
        } else {
            low = icell;
        }
        icell = (low + high) / 2;
    }
    return icell;
}

template <class Tag, std::size_t D>
template <class MemorySpace>
ChunkSpan<double, DiscreteDomain<NonUniformBSplines<Tag, D>>> NonUniformBSplines<Tag, D>::Impl<
        MemorySpace>::integrals(ChunkSpan<double, DiscreteDomain<NonUniformBSplines<Tag, D>>>
                                        int_vals) const
{
    if constexpr (is_periodic()) {
        assert(int_vals.size() == nbasis() || int_vals.size() == size());
    } else {
        assert(int_vals.size() == nbasis());
    }

    double const inv_deg = 1.0 / (degree() + 1);

    discrete_domain_type const dom_bsplines(
            full_domain().take_first(discrete_vector_type {nbasis()}));
    for (auto ix : dom_bsplines) {
        int_vals(ix) = (get_last_support_knot(ix) - get_first_support_knot(ix)) * inv_deg;
    }

    if constexpr (is_periodic()) {
        if (int_vals.size() == size()) {
            discrete_domain_type const dom_bsplines_wrap(
                    full_domain().take_last(discrete_vector_type {degree()}));
            for (auto ix : dom_bsplines_wrap) {
                int_vals(ix) = 0;
            }
        }
    }
    return int_vals;
}
