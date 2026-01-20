// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "tensor.hpp"

template <class StartDim, class EndDim>
class LinearCoordTransform
{
public:
    using CoordArg = Coord<StartDim>;
    using CoordResult = Coord<EndDim>;
    using CoordJacobian = CoordArg;

private:
    Coord<StartDim> m_start_reference;
    Coord<EndDim> m_end_reference;
    double m_scaling_factor;

public:
    explicit KOKKOS_FUNCTION LinearCoordTransform(
            Coord<StartDim> start_reference,
            Coord<EndDim> end_reference,
            double scaling_factor)
        : m_start_reference(start_reference)
        , m_end_reference(end_reference)
        , m_scaling_factor(scaling_factor)
    {
    }

    KOKKOS_DEFAULTED_FUNCTION LinearCoordTransform(LinearCoordTransform const&) = default;

    KOKKOS_DEFAULTED_FUNCTION ~LinearCoordTransform() = default;

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        double arg_dist = coord - m_start_reference;
        double res_dist = arg_dist * m_scaling_factor;
        return m_end_reference + res_dist;
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(CoordArg const& coord) const
    {
        return m_scaling_factor;
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the function jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<EndDim>, VectorIndexSet<StartDim>> jacobian_matrix(
            CoordArg const& coord) const
    {
        DTensor<VectorIndexSet<EndDim>, VectorIndexSet<StartDim>> jacobian_matrix;
        ddcHelper::get<EndDim, StartDim>(jacobian_matrix) = m_scaling_factor;
        return jacobian_matrix;
    }

    /**
     * @brief Compute the (i,j) coefficient of the Jacobian matrix.
     * 
     * @param[in] coord
     *              The coordinate where we evaluate the Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(std::is_same_v<IndexTag1, EndDim>);
        static_assert(std::is_same_v<IndexTag2, StartDim>);
        return m_scaling_factor;
    }


    /**
     * @brief Compute full inverse Jacobian matrix.
     *
     * For some computations, we need the complete inverse Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given independently with the function inv_jacobian_component.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The inverse Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<VectorIndexSet<StartDim>, VectorIndexSet<EndDim>> inv_jacobian_matrix(
            CoordArg const& coord) const
    {
        DTensor<VectorIndexSet<StartDim>, VectorIndexSet<EndDim>> matrix;
        ddcHelper::get<StartDim, EndDim>(matrix) = 1.0 / m_scaling_factor;
        return matrix;
    }


    /**
     * @brief Compute the (i,j) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
     *
     * @param[in] coord
     *              The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (i,j) coefficient of the inverse Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double inv_jacobian_component(CoordArg const& coord) const
    {
        static_assert(std::is_same_v<IndexTag1, StartDim>);
        static_assert(std::is_same_v<IndexTag2, EndDim>);
        return 1.0 / m_scaling_factor;
    }

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    KOKKOS_INLINE_FUNCTION LinearCoordTransform<EndDim, StartDim> get_inverse_mapping() const
    {
        return LinearCoordTransform<
                EndDim,
                StartDim>(m_end_reference, m_start_reference, 1.0 / m_scaling_factor);
    }
};
