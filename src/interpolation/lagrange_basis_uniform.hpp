// SPDX-License-Identifier: MIT
#pragma once
#include <array>

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
    };
};

template <class DDim>
struct is_uniform_lagrange_basis
    : public std::is_base_of<detail::UniformLagrangeBasisBase, DDim>::type
{
};
