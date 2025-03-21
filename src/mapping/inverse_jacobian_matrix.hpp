// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "math_tools.hpp"
#include "view.hpp"

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
    using ValidArgIndices = ddc::to_type_seq_t<typename Mapping::CoordArg>;
    using ValidResultIndices = ddc::to_type_seq_t<typename Mapping::CoordResult>;
    using DimArg0 = ddc::type_seq_element_t<0, ValidArgIndices>;
    using DimArg1 = ddc::type_seq_element_t<1, ValidArgIndices>;
    using DimRes0_cov = typename ddc::type_seq_element_t<0, ValidResultIndices>::Dual;
    using DimRes1_cov = typename ddc::type_seq_element_t<1, ValidResultIndices>::Dual;

public:
    /// The type of the tensor representing the inverse Jacobian.
    using InverseJacobianTensor
            = DTensor<ValidArgIndices, vector_index_set_dual_t<ValidResultIndices>>;

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
     * The coefficients can be given independently with the functions
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
    KOKKOS_INLINE_FUNCTION InverseJacobianTensor operator()(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_matrix(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            DTensor<ValidResultIndices, vector_index_set_dual_t<ValidArgIndices>> jacobian
                    = m_mapping.jacobian_matrix(coord);
            assert(fabs(determinant(jacobian)) > 1e-15);
            return inverse(jacobian);
        }
    }

    /**
     * @brief Compute the (1,1) coefficient of the inverse Jacobian matrix.
     *
     * Be careful because not all mappings are invertible, especially at the centre point.
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
     * Be careful because not all mappings are invertible, especially at the centre point.
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
     * Be careful because not all mappings are invertible, especially at the centre point.
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
     * Be careful because not all mappings are invertible, especially at the centre point.
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
