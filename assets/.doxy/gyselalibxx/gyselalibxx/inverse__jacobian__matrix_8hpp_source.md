

# File inverse\_jacobian\_matrix.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**inverse\_jacobian\_matrix.hpp**](inverse__jacobian__matrix_8hpp.md)

[Go to the documentation of this file](inverse__jacobian__matrix_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "math_tools.hpp"
#include "view.hpp"

template <class Mapping, class PositionCoordinate = typename Mapping::CoordArg>
class InverseJacobianMatrix
{
private:
    using ValidArgIndices = ddc::to_type_seq_t<typename Mapping::CoordArg>;
    using ValidResultIndices = ddc::to_type_seq_t<typename Mapping::CoordResult>;
    using DimArg0 = ddc::type_seq_element_t<0, ValidArgIndices>;
    using DimArg1 = ddc::type_seq_element_t<1, ValidArgIndices>;
    using DimArg0_cov = typename DimArg0::Dual;
    using DimArg1_cov = typename ddc::type_seq_element_t<1, ValidArgIndices>::Dual;
    using DimRes0 = typename ddc::type_seq_element_t<0, ValidResultIndices>;
    using DimRes1 = typename ddc::type_seq_element_t<1, ValidResultIndices>;
    using DimRes0_cov = typename ddc::type_seq_element_t<0, ValidResultIndices>::Dual;
    using DimRes1_cov = typename ddc::type_seq_element_t<1, ValidResultIndices>::Dual;

public:
    using InverseJacobianTensor
            = DTensor<ValidArgIndices, vector_index_set_dual_t<ValidResultIndices>>;

private:
    Mapping m_mapping;

public:
    KOKKOS_FUNCTION explicit InverseJacobianMatrix(Mapping const& mapping) : m_mapping(mapping) {}

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

    KOKKOS_INLINE_FUNCTION double inv_jacobian_11(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_11(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            // J^{-1}(1,1) = J(2,2) / det(J)
            return m_mapping.template jacobian_component<DimRes1, DimArg1_cov>(coord) / jacob;
        }
    }

    KOKKOS_INLINE_FUNCTION double inv_jacobian_12(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_12(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            // J^{-1}(1,2) = -J(1,2) / det(J)
            return -m_mapping.template jacobian_component<DimRes0, DimArg1_cov>(coord) / jacob;
        }
    }
    KOKKOS_INLINE_FUNCTION double inv_jacobian_21(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_21(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            // J^{-1}(2,1) = -J(2,1) / det(J)
            return -m_mapping.template jacobian_component<DimRes1, DimArg0_cov>(coord) / jacob;
        }
    }

    KOKKOS_INLINE_FUNCTION double inv_jacobian_22(PositionCoordinate const& coord) const
    {
        if constexpr (has_2d_inv_jacobian_v<Mapping, PositionCoordinate>) {
            return m_mapping.inv_jacobian_22(coord);
        } else {
            static_assert(has_2d_jacobian_v<Mapping, PositionCoordinate>);
            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            // J^{-1}(2,2) = J(1,1) / det(J)
            return m_mapping.template jacobian_component<DimRes0, DimArg0_cov>(coord) / jacob;
        }
    }
};
```


