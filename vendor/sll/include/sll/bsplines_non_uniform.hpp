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
    // From nvcc: 'A type that is defined inside a class and has private or protected access cannot be used
    // in the template argument type of a variable template instantiation'
    template <class T>
    struct InternalTagGenerator;

    /// An internal tag necessary to allocate an internal ddc::discrete_space function.
    /// It must remain internal, for example it shall not be exposed when returning ddc::coordinates. Instead use `Tag`
    using KnotDim = InternalTagGenerator<Tag>;

    using mesh_type = ddc::NonUniformPointSampling<KnotDim>;

    static inline ddc::Coordinate<KnotDim> knot_from_coord(ddc::Coordinate<Tag> const& coord)
    {
        return ddc::Coordinate<KnotDim>(ddc::get<Tag>(coord));
    }
    static inline ddc::Coordinate<Tag> coord_from_knot(ddc::Coordinate<KnotDim> const& coord)
    {
        return ddc::Coordinate<Tag>(ddc::get<KnotDim>(coord));
    }

public:
    using tag_type = Tag;

    using continuous_dimension_type = BSpline<Tag>;


    using discrete_dimension_type = NonUniformBSplines;

    using discrete_element_type = ddc::DiscreteElement<NonUniformBSplines>;

    using discrete_domain_type = ddc::DiscreteDomain<NonUniformBSplines>;

    using discrete_vector_type = ddc::DiscreteVector<NonUniformBSplines>;

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
        ddc::DiscreteDomain<mesh_type> m_domain;
        const int m_nknots;

    public:
        using discrete_dimension_type = NonUniformBSplines;

        Impl() = default;

        template <class OriginMemorySpace>
        explicit Impl(Impl<OriginMemorySpace> const& impl)
            : m_domain(impl.m_domain)
            , m_nknots(impl.m_nknots)
        {
        }

        /// @brief Construct a `Impl` using a brace-list, i.e. `Impl bsplines({0., 1.})`
        explicit Impl(std::initializer_list<ddc::Coordinate<Tag>> breaks)
            : Impl(breaks.begin(), breaks.end())
        {
        }

        /// @brief Construct a `Impl` using a C++20 "common range".
        explicit Impl(std::vector<ddc::Coordinate<Tag>> const& breaks)
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

        discrete_element_type eval_basis(DSpan1D values, ddc::Coordinate<Tag> const& x) const;

        discrete_element_type eval_deriv(DSpan1D derivs, ddc::Coordinate<Tag> const& x) const;

        discrete_element_type eval_basis_and_n_derivs(
                DSpan2D derivs,
                ddc::Coordinate<Tag> const& x,
                std::size_t n) const;

        ddc::ChunkSpan<double, discrete_domain_type> integrals(
                ddc::ChunkSpan<double, discrete_domain_type> int_vals) const;

        ddc::Coordinate<Tag> get_knot(int knot_idx) const noexcept
        {
            return coord_from_knot(
                    ddc::coordinate(ddc::DiscreteElement<mesh_type>(knot_idx + degree())));
        }

        ddc::Coordinate<Tag> get_first_support_knot(discrete_element_type const& ix) const
        {
            return coord_from_knot(ddc::coordinate(ddc::DiscreteElement<mesh_type>(ix.uid())));
        }

        ddc::Coordinate<Tag> get_last_support_knot(discrete_element_type const& ix) const
        {
            return coord_from_knot(
                    ddc::coordinate(ddc::DiscreteElement<mesh_type>(ix.uid() + degree() + 1)));
        }

        ddc::Coordinate<Tag> get_support_knot_n(discrete_element_type const& ix, int n) const
        {
            return coord_from_knot(ddc::coordinate(ddc::DiscreteElement<mesh_type>(ix.uid() + n)));
        }

        ddc::Coordinate<Tag> rmin() const noexcept
        {
            return get_knot(0);
        }

        ddc::Coordinate<Tag> rmax() const noexcept
        {
            return get_knot(ncells());
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
            return m_nknots - 2 * degree();
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
        int find_cell(ddc::Coordinate<Tag> const& x) const;
    };
};

template <class Tag, std::size_t D>
template <class MemorySpace>
template <class RandomIt>
NonUniformBSplines<Tag, D>::Impl<MemorySpace>::Impl(
        RandomIt const break_begin,
        RandomIt const break_end)
    : m_domain(
            ddc::DiscreteElement<mesh_type>(0),
            ddc::DiscreteVector<mesh_type>(
                    (break_end - break_begin)
                    + 2 * degree())) // Create a mesh including the eventual periodic point
    , m_nknots((break_end - break_begin) + 2 * degree())
{
    std::vector<ddc::Coordinate<KnotDim>> knots((break_end - break_begin) + 2 * degree());
    // Fill the provided knots
    int ii = 0;
    for (RandomIt it = break_begin; it < break_end; ++it) {
        knots[degree() + ii] = knot_from_coord(*it);
        ++ii;
    }
    ddc::Coordinate<KnotDim> const rmin = knots[degree()];
    ddc::Coordinate<KnotDim> const rmax = knots[(break_end - break_begin) + degree() - 1];
    assert(rmin < rmax);

    // Fill out the extra knots
    if constexpr (is_periodic()) {
        double const period = rmax - rmin;
        for (std::size_t i = 1; i < degree() + 1; ++i) {
            knots[degree() + -i]
                    = ddc::Coordinate<KnotDim>(knots[degree() + ncells() - i] - period);
            knots[degree() + ncells() + i] = ddc::Coordinate<KnotDim>(knots[degree() + i] + period);
        }
    } else // open
    {
        for (std::size_t i = 1; i < degree() + 1; ++i) {
            knots[degree() + -i] = rmin;
            knots[degree() + npoints() - 1 + i] = rmax;
        }
    }
    ddc::init_discrete_space<mesh_type>(knots);
}

template <class Tag, std::size_t D>
template <class MemorySpace>
ddc::DiscreteElement<NonUniformBSplines<Tag, D>> NonUniformBSplines<Tag, D>::Impl<
        MemorySpace>::eval_basis(DSpan1D const values, ddc::Coordinate<Tag> const& x) const
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
ddc::DiscreteElement<NonUniformBSplines<Tag, D>> NonUniformBSplines<Tag, D>::Impl<
        MemorySpace>::eval_deriv(DSpan1D const derivs, ddc::Coordinate<Tag> const& x) const
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
ddc::DiscreteElement<NonUniformBSplines<Tag, D>> NonUniformBSplines<Tag, D>::Impl<MemorySpace>::
        eval_basis_and_n_derivs(
                DSpan2D const derivs,
                ddc::Coordinate<Tag> const& x,
                std::size_t const n) const
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
    // assert(n >= 0); as long as n is unsigned
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
int NonUniformBSplines<Tag, D>::Impl<MemorySpace>::find_cell(ddc::Coordinate<Tag> const& x) const
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
ddc::ChunkSpan<double, ddc::DiscreteDomain<NonUniformBSplines<Tag, D>>> NonUniformBSplines<Tag, D>::
        Impl<MemorySpace>::integrals(
                ddc::ChunkSpan<double, ddc::DiscreteDomain<NonUniformBSplines<Tag, D>>> int_vals)
                const
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
