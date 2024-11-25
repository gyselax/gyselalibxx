// SPDX-License-Identifier: MIT
#pragma once

#include <sll/view.hpp>

#include "mapping_tools.hpp"

/**
 * A class to calculate the inverse of the Jacobian matrix.
 * If specialised methods are available then these are used by the class.
 * Otherwise the inverse is calculated from the Jacobian matrix.
 *
 * @tparam Mapping The mapping whose inverse we are interested in.
 * @tparam PositionCoordinate The coordinate system in which the inverse should be calculated.
 */
template <class Mapping, class PositionCoordinate = typename Mapping::CoordArg>
class InverseJacobianMatrix
{
private:
    Mapping m_mapping;

public:
    /**
     * @brief A constructor for the InverseJacobianMatrix.
     * @param[in] mapping The mapping whose inverse we are interested in.
     */
    KOKKOS_FUNCTION explicit InverseJacobianMatrix(Mapping const& mapping) : m_mapping(mapping) {}

    /**
     * @brief Compute full inverse Jacobian matrix.
     *
     * For some computations, we need the complete inverse Jacobian matrix or just the
     * coefficients.
     * The coefficients can be given indendently with the functions
     * inv_jacobian_11, inv_jacobian_12, inv_jacobian_21 and inv_jacobian_22.
     *
     * @param[in] coord The coordinate where we evaluate the Jacobian matrix.
     * @returns The inverse Jacobian matrix returned.
     *
     * @see inv_jacobian_11
     * @see inv_jacobian_12
     * @see inv_jacobian_21
     * @see inv_jacobian_22
     */
    KOKKOS_INLINE_FUNCTION Matrix_2x2 operator()(PositionCoordinate const& coord) const
    {
        Matrix_2x2 matrix;
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            m_mapping.inv_jacobian_matrix(coord, matrix);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            matrix[0][0] = m_mapping.jacobian_22(coord) / jacob;
            matrix[0][1] = -m_mapping.jacobian_12(coord) / jacob;
            matrix[1][0] = -m_mapping.jacobian_21(coord) / jacob;
            matrix[1][1] = m_mapping.jacobian_11(coord) / jacob;
        }
        return matrix;
    }

    /**
     * @brief Compute the (1,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (1,1) coefficient of the inverse Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION double inv_jacobian_11(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_11(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            return m_mapping.jacobian_22(coord) / jacob;
        }
    }

    /**
     * @brief Compute the (1,2) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (1,2) coefficient of the inverse Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION double inv_jacobian_12(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_12(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            return -m_mapping.jacobian_12(coord) / jacob;
        }
    }
    /**
     * @brief Compute the (2,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (2,1) coefficient of the inverse Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION double inv_jacobian_21(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_21(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            return -m_mapping.jacobian_21(coord) / jacob;
        }
    }

    /**
     * @brief Compute the (2,2) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the center point.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the inverse Jacobian matrix.
     *
     * @return A double with the value of the (2,2) coefficient of the inverse Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION double inv_jacobian_22(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_22(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            return m_mapping.jacobian_11(coord) / jacob;
        }
    }
};
