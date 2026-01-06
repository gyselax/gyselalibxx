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
 *
 * @tparam Grid1D The grid on which the Lagrange polynomials
 *                  are defined.
 * @tparam D The degree of the polynomials, equal to the number of
 *                  cells over which the Lagrange polynomials are
 *                  defined.
 * @tparam DataType The data type used for the calculations.
 *                  Double by default.
 */
template <class Grid1D, std::size_t D, class DataType = double>
class UniformLagrangeBasis : detail::UniformLagrangeBasisBase
{
    static_assert(D > 0, "Parameter `D` must be positive");
    static_assert(std::is_floating_point_v<DataType>);

public:
    /// @brief The tag identifying the continuous dimension on which the Lagrange polynomials are defined.
    using continuous_dimension_type = typename Grid1D::continuous_dimension_type;

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
        using knot_discrete_dimension_type = UniformLagrangeKnots<DDim>;

    private:
        std::array<DataType, D + 1> m_weights;

        IdxRange<knot_discrete_dimension_type> m_knot_domain;
        IdxRange<knot_discrete_dimension_type> m_break_point_domain;

        Idx<DDim> m_reference;

    public:
        Impl() = default;

        explicit Impl(IdxRange<Grid1D> break_point_domain);

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


        //discrete_element_type get_icell_and_offset(
        //        int& icell,
        //        DataType& offset,
        //        ddc::Coordinate<CDim> const& x) const;

        KOKKOS_INLINE_FUNCTION void eval_basis(
                Span1D<DataType> values,
                coord_type const& x,
                Idx<knot_discrete_dimension_type> poly_start) const
        {
            std::cout << values.size() << " " << degree() + 1 << std::endl;
            KOKKOS_ASSERT(values.size() == degree() + 1);
            //discrete_element_type poly_start = get_icell_and_offset(first_elem, icell, offset, x);
            DataType dx = ddc::discrete_space<knot_discrete_dimension_type>().step();
            double inv_dx = 1. / dx;
            DataType offset = (x - ddc::coordinate(poly_start)) * inv_dx;
            DataType eps = std::numeric_limits<DataType>::epsilon() * 4;
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

template <class Grid1D, std::size_t D, class DataType>
template <class DDim, class MemorySpace>
UniformLagrangeBasis<Grid1D, D, DataType>::Impl<DDim, MemorySpace>::Impl(
        IdxRange<Grid1D> break_point_domain)
    : m_reference(ddc::create_reference_discrete_element<DDim>())
{
    std::size_t ncells = break_point_domain.size() - 1;
    assert(ncells >= D);

    coord_type bp_rmin = ddc::coordinate(break_point_domain.front());
    DataType step = ddc::discrete_space<Grid1D>().step();

    // Initialise knot grid
    if constexpr (is_periodic()) {
        coord_type rmin_local = bp_rmin - step * (D - 1);
        ddc::init_discrete_space<knot_discrete_dimension_type>(rmin_local, step);
        m_knot_domain = IdxRange<knot_discrete_dimension_type>(
                ddc::discrete_space<knot_discrete_dimension_type>().front(),
                IdxStep<knot_discrete_dimension_type>(break_point_domain.size() + 2 * (D - 1)));
        m_break_point_domain = m_knot_domain
                                       .remove(IdxStep<knot_discrete_dimension_type>(D - 1),
                                               IdxStep<knot_discrete_dimension_type>(D - 1));
    } else {
        ddc::init_discrete_space<knot_discrete_dimension_type>(bp_rmin, step);
        m_knot_domain = IdxRange<knot_discrete_dimension_type>(
                ddc::discrete_space<knot_discrete_dimension_type>().front(),
                IdxStep<knot_discrete_dimension_type>(break_point_domain.size()));
        m_break_point_domain = m_knot_domain;
    }

    // Calculate weights
    DataType dx = ddc::discrete_space<knot_discrete_dimension_type>().step();
    for (std::size_t i(0); i < D + 1; ++i) {
        DataType numerator = dx;
        for (std::size_t j(0); j < D + 1; ++j) {
            if (i == j)
                continue;
            numerator *= (i - j);
        }
        m_weights[i] = 1.0 / numerator;
    }
}
