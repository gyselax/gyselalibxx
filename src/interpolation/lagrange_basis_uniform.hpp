// SPDX-License-Identifier: MIT
#pragma once
#include <array>
#include <cstddef>
#include <type_traits>

#include "ddc_aliases.hpp"
#include "view.hpp"

namespace detail {

class UniformLagrangeBasisBase
{
};

} // namespace detail

template <class T>
struct UniformLagrangeKnots : UniformGridBase<typename T::continuous_dimension_type>
{
};

/**
 * @brief Class describing Lagrange polynomials on a uniform grid.
 *
 * This class uses the second barycentric formulation to evaluate the
 * polynomials. This formula is used for stability.
 * It is described in
 * Barycentric Lagrange Interpolation
 * Jean-Paul Berrut and Lloyd N. Trefethen
 * SIAM Review 2004 46:3, 501-517
 *
 * @tparam Dim The dimension on which the Lagrange polynomials
 *                  are defined.
 * @tparam D The degree of the polynomials, equal to the number of
 *                  cells over which the Lagrange polynomials are
 *                  defined.
 * @tparam DataType The data type used for the calculations.
 *                  Double by default.
 */
template <class Dim, std::size_t D, class DataType = double>
class UniformLagrangeBasis : detail::UniformLagrangeBasisBase
{
    static_assert(D > 0, "Parameter `D` must be positive");
    static_assert(std::is_floating_point_v<DataType>);

public:
    /// @brief The tag identifying the continuous dimension on which the Lagrange polynomials are defined.
    using continuous_dimension_type = Dim;

    /// @brief The type of the coordinates on which the Lagrange polynomials can be evaluated.
    using coord_type = Coord<continuous_dimension_type>;

    /// @brief The discrete dimension representing B-splines.
    using discrete_dimension_type = UniformLagrangeBasis;

    /** @brief The degree of the Lagrange polynomials.
     *
     * @return The degree.
     */
    static constexpr std::size_t degree() noexcept
    {
        return D;
    }

    /** @brief Indicates if the Lagrange polynomials are periodic or not.
     *
     * @return A boolean indicating if the Lagrange polynomials are periodic or not.
     */
    static constexpr bool is_periodic() noexcept
    {
        return continuous_dimension_type::PERIODIC;
    }

    /** @brief Indicates if the Lagrange polynomials are uniform or not (this is the case here).
     *
     * @return A boolean indicating if the Lagrange polynomials are uniform or not.
     */
    static constexpr bool is_uniform() noexcept
    {
        return true;
    }


    /** @brief Storage class of the static attributes of the discrete dimension.
     *
     * @tparam DDim The name of the discrete dimension.
     * @tparam MemorySpace The Kokkos memory space where the attributes are being stored.
     */
    template <class DDim, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

    public:
        /// @brief The type of the knots defining the B-splines.
        using knot_grid = UniformLagrangeKnots<DDim>;

    private:
        std::array<DataType, D + 1> m_weights;

        IdxRange<knot_grid> m_knot_domain;
        IdxRange<knot_grid> m_break_point_domain;

        Idx<DDim> m_reference;

    public:
        Impl() = default;

        /** @brief Initialise the possible Lagrange bases.
         *
         * Initialise the class such that Lagrange bases can be evaluated on domains
         * derived from the break point domain.
         */
        template <class Grid1D>
        explicit Impl(IdxRange<Grid1D> break_point_domain)
            : m_reference(ddc::create_reference_discrete_element<DDim>())
        {
            static_assert(std::is_same_v<typename Grid1D::continuous_dimension_type, Dim>);
            static_assert(ddc::is_uniform_point_sampling<Grid1D>());
            assert(ddc::is_discrete_space_initialized<Grid1D>());
            assert(break_point_domain.size() >= D + 1);

            coord_type bp_rmin = ddc::coordinate(break_point_domain.front());
            DataType step = ddc::discrete_space<Grid1D>().step();

            ddc::init_discrete_space<knot_grid>(bp_rmin, step);
            // Initialise knot grid
            if constexpr (is_periodic()) {
                m_knot_domain = IdxRange<knot_grid>(
                        ddc::discrete_space<knot_grid>().front(),
                        IdxStep<knot_grid>(break_point_domain.size() + D));
                m_break_point_domain = m_knot_domain.remove_last(IdxStep<knot_grid>(D));
            } else {
                m_knot_domain = IdxRange<knot_grid>(
                        ddc::discrete_space<knot_grid>().front(),
                        IdxStep<knot_grid>(break_point_domain.size()));
                m_break_point_domain = m_knot_domain;
            }

            // Calculate weights
            DataType dx = ddc::discrete_space<knot_grid>().step();
            for (int i(0); i < D + 1; ++i) {
                DataType numerator = dx;
                for (int j(0); j < D + 1; ++j) {
                    if (i == j)
                        continue;
                    numerator *= (i - j);
                }
                m_weights[i] = 1.0 / numerator;
            }
        }

        /** @brief Copy-constructs from another Impl with a different Kokkos memory space.
         *
         * @param impl A reference to the other Impl.
         */
        template <class OriginMemorySpace>
        explicit Impl(Impl<DDim, OriginMemorySpace> const& impl)
            : m_weights(impl.m_weights)
            , m_knot_domain(impl.m_knot_domain)
            , m_break_point_domain(impl.m_break_point_domain)
            , m_reference(impl.m_reference)
        {
        }

        /** @brief Returns the coordinate of the lower bound of the domain on which the B-splines are defined.
         *
         * @return Coordinate of the lower bound of the domain.
         */
        KOKKOS_INLINE_FUNCTION coord_type rmin() const noexcept
        {
            return ddc::coordinate(m_break_point_domain.front());
        }

        /** @brief Returns the coordinate of the upper bound of the domain on which the B-splines are defined.
         *
         * @return Coordinate of the upper bound of the domain.
         */
        KOKKOS_INLINE_FUNCTION coord_type rmax() const noexcept
        {
            return ddc::coordinate(m_break_point_domain.back());
        }

        /** @brief Returns the length of the domain.
         *
         * @return The length of the domain.
         */
        KOKKOS_INLINE_FUNCTION DataType length() const noexcept
        {
            return static_cast<DataType>(rmax()) - static_cast<DataType>(rmin());
        }

        /** @brief Returns the index range of the break points.
         *
         * @return The index range describing the break points.
         */
        KOKKOS_INLINE_FUNCTION IdxRange<knot_grid> break_point_domain() const
        {
            return m_break_point_domain;
        }

        /** @brief Returns the index range including eventual duplicate values for the periodic case.
         *
         * @return The index range including eventual duplicate values.
         */
        KOKKOS_INLINE_FUNCTION IdxRange<knot_grid> full_domain() const
        {
            return m_knot_domain;
        }

        /** @brief Evaluate the selected set of bases at the coordinate.
         *
         * Evaluate all d+1 bases which span the domain
         * [coordinate(poly_start), coordinate(poly_start+d)]
         * at the coordinate x.
         *
         * @param[out] values The values of each basis at the coordinate x.
         * @param[in] x The coordinate where the bases are evaluated.
         * @param[in] poly_start The index of the first of the d+1 knots describing
         *                       the set of bases to be evaluated.
         */
        KOKKOS_INLINE_FUNCTION void eval_basis(
                Span1D<DataType> values,
                coord_type x,
                Idx<knot_grid> poly_start) const
        {
            KOKKOS_ASSERT(values.size() == degree() + 1);
            if constexpr (is_periodic()) {
                if (x < ddc::coordinate(poly_start)) {
                    x = x + length();
                }
            }
            KOKKOS_ASSERT(x >= ddc::coordinate(poly_start));
            KOKKOS_ASSERT(x <= ddc::coordinate(poly_start + degree()));

            DataType dx = ddc::discrete_space<knot_grid>().step();
            double inv_dx = 1. / dx;
            DataType offset = (x - ddc::coordinate(poly_start)) * inv_dx;
            DataType eps = Kokkos::Experimental::epsilon_v<DataType> * 4;
            int icell = static_cast<int>(offset);
            DataType local_offset = offset - icell;
            if (local_offset < eps) {
                for (int j(0); j < D + 1; ++j) {
                    values(j) = static_cast<int>(j == icell);
                }
            } else if (local_offset > (1 - eps)) {
                for (int j(0); j < D + 1; ++j) {
                    values(j) = static_cast<int>(j == (icell + 1));
                }
            } else {
                for (int i(0); i < D + 1; ++i) {
                    DataType numerator = m_weights[i] / (dx * (offset - i));
                    DataType denominator(0.0);
                    for (int j(0); j < D + 1; ++j) {
                        denominator += m_weights[j] / (dx * (offset - j));
                    }
                    values(i) = numerator / denominator;
                }
            }
        }

        /**
         * @brief Evaluate the selected set of bases and their derivatives at the coordinate.
         *
         * Evaluate all d+1 bases which span the domain
         * [coordinate(poly_start), coordinate(poly_start+d)]
         * at the coordinate x.
         * The derivatives are also calculated. The calculations are derived using
         * the method described in:
         *   Barycentric Lagrange Interpolation
         *   Jean-Paul Berrut and Lloyd N. Trefethen
         *   SIAM Review 2004 46:3, 501-517
         *
         * @param[out] derivs The values and derivatives of each basis at the coordinate x.
         * @param[in] x The coordinate where the bases are evaluated.
         * @param[in] poly_start The index of the first of the d+1 knots describing
         *                       the set of bases to be evaluated.
         * @param[in] n_derivs The number of derivatives to be calculated.
         */
        KOKKOS_INLINE_FUNCTION void eval_basis_and_n_derivs(
                Span2D<DataType> derivs,
                coord_type const& x,
                Idx<knot_grid> poly_start,
                std::size_t n_derivs) const
        {
            KOKKOS_ASSERT(x >= ddc::coordinate(poly_start));
            KOKKOS_ASSERT(x <= ddc::coordinate(poly_start + degree()));
            KOKKOS_ASSERT(n_derivs <= degree())
            KOKKOS_ASSERT(derivs.extent(0) == degree() + 1)
            KOKKOS_ASSERT(derivs.extent(1) == n_derivs + 1)

            constexpr std::size_t n_basis = degree() + 1;
            DataType eps = Kokkos::Experimental::epsilon_v<DataType> * 4;

            int node(-1);
            for (std::size_t i(0); i < n_basis; ++i) {
                if (Kokkos::fabs(x - ddc::coordinate(poly_start + i)) < eps) {
                    node = i;
                    for (int j(0); j < n_basis; ++j) {
                        derivs(j, 0) = static_cast<int>(j == i);
                    }
                }
            }

            // If coordinate not found at a node
            if (node == -1) {
                std::array<DataType, n_basis> vals_ptr;
                Span1D<DataType> values(vals_ptr.data(), n_basis);
                eval_basis(values, x, poly_start);

                std::array<int, n_basis> combinations;

                // For each basis element
                for (std::size_t j(0); j < n_basis; ++j) {
                    derivs(j, 0) = values[j];
                    for (std::size_t n(1); n < n_derivs + 1; ++n) {
                        DataType factor = 0;
                        combinations[0] = -1;
                        for (std::size_t i(1); i < n; ++i) {
                            combinations[i] = 0;
                        }
                        int i = 0;
                        while (i > -1) {
                            combinations[i] += 1;
                            combinations[i] += static_cast<int>(combinations[i] == j);
                            for (std::size_t k(i + 1); k < n; ++k) {
                                combinations[k] = combinations[k - 1] + 1;
                                combinations[k] += static_cast<int>(combinations[k] == j);
                            }
                            i = n - 1;
                            double divisor(1);
                            for (std::size_t k(0); k < n; ++k) {
                                divisor *= (x - ddc::coordinate(poly_start + combinations[k]));
                            }
                            factor += 1.0 / divisor;
                            int max_val_along_dim = i + n_basis - n;
                            max_val_along_dim -= static_cast<int>(max_val_along_dim <= j);
                            while (i > -1 and combinations[i] == max_val_along_dim) {
                                i--;
                                max_val_along_dim = i + n_basis - n;
                                max_val_along_dim -= static_cast<int>(max_val_along_dim <= j);
                            }
                        }
                        derivs(j, n) = values(j) * factor * factorial(n);
                    }
                }
            } else {
                DataType xi = ddc::coordinate(poly_start + node);
                for (std::size_t j(0); j < n_basis; ++j) {
                    DataType xj = ddc::coordinate(poly_start + j);
                    if (j == node)
                        continue;
                    for (std::size_t n(1); n < n_derivs + 1; ++n) {
                        // \partial^nl_j(x_i) = n!/wi [ (-1)^(n+1) wj/(xi-xj)^n
                        //                      - \sum_k=1^{n-1}\sum_{p\ne i} (-1)^(n+1-k)
                        //                                       wp \partial^k l_j(x_i) / k! / (xi-xp)^n ]
                        derivs(j, n) = neg_1_pow(n + 1) * m_weights[j] / Kokkos::pow(xi - xj, n);
                        for (std::size_t p(0); p < n_basis; ++p) {
                            if (p == node)
                                continue;
                            DataType xp = ddc::coordinate(poly_start + p);
                            for (std::size_t k(1); k < n; ++k) {
                                derivs(j, n) -= neg_1_pow(n + 1 - k) * m_weights[p] * derivs(j, k)
                                                / factorial(k) / ipow(xi - xp, n - k);
                            }
                        }
                        derivs(j, n) *= factorial(n) / m_weights[node];
                    }
                }
                for (std::size_t n(1); n < n_derivs + 1; ++n) {
                    derivs(node, n) = 0;
                    for (std::size_t j(0); j < n_basis; ++j) {
                        if (j == node)
                            continue;
                        derivs(node, n) -= derivs(j, n);
                    }
                }
            }
        }
    };
};

template <class DDim>
struct is_uniform_lagrange_basis
    : public std::is_base_of<detail::UniformLagrangeBasisBase, DDim>::type
{
};

template <class DDim>
constexpr bool is_uniform_lagrange_basis_v = is_uniform_lagrange_basis<DDim>::value;
