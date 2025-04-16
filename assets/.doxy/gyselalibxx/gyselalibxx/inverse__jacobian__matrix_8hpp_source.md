

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

public:
    using InverseJacobianTensor
            = DTensor<ValidArgIndices, vector_index_set_dual_t<ValidResultIndices>>;

private:
    Mapping m_mapping;

public:
    KOKKOS_FUNCTION explicit InverseJacobianMatrix(Mapping const& mapping) : m_mapping(mapping) {}

    KOKKOS_INLINE_FUNCTION InverseJacobianTensor operator()(PositionCoordinate const& coord) const
    {
        if constexpr (has_inv_jacobian_v<Mapping, PositionCoordinate, false>) {
            return m_mapping.inv_jacobian_matrix(coord);
        } else {
            static_assert(has_jacobian_v<Mapping, PositionCoordinate>);
            DTensor<ValidResultIndices, vector_index_set_dual_t<ValidArgIndices>> jacobian
                    = m_mapping.jacobian_matrix(coord);
            assert(fabs(determinant(jacobian)) > 1e-15);
            return inverse(jacobian);
        }
    }


    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double inv_jacobian_component(PositionCoordinate const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ValidArgIndices>);
        static_assert(ddc::in_tags_v<IndexTag2, get_covariant_dims_t<ValidResultIndices>>);

        if constexpr (has_inv_jacobian_v<Mapping, PositionCoordinate, false>) {
            return m_mapping.template inv_jacobian_component<IndexTag1, IndexTag2>(coord);
        } else {
            static_assert(has_jacobian_v<Mapping, PositionCoordinate>);
            static_assert(Mapping::CoordArg::size() == 2);

            using DimArg0 = ddc::type_seq_element_t<0, ValidArgIndices>;
            using DimArg1 = ddc::type_seq_element_t<1, ValidArgIndices>;
            using DimArg0_cov = typename DimArg0::Dual;
            using DimArg1_cov = typename ddc::type_seq_element_t<1, ValidArgIndices>::Dual;
            using DimRes0 = typename ddc::type_seq_element_t<0, ValidResultIndices>;
            using DimRes1 = typename ddc::type_seq_element_t<1, ValidResultIndices>;
            using DimRes0_cov = typename ddc::type_seq_element_t<0, ValidResultIndices>::Dual;
            using DimRes1_cov = typename ddc::type_seq_element_t<1, ValidResultIndices>::Dual;

            double jacob = m_mapping.jacobian(coord);
            assert(fabs(jacob) > 1e-15);
            if constexpr (
                    std::is_same_v<IndexTag1, DimArg0> && std::is_same_v<IndexTag2, DimRes0_cov>) {
                //Compute the (1,1) coefficient of the inverse Jacobian matrix.
                // J^{-1}(1,1) = J(2,2) / det(J)
                return m_mapping.template jacobian_component<DimRes1, DimArg1_cov>(coord) / jacob;
            } else if constexpr (
                    std::is_same_v<IndexTag1, DimArg0> && std::is_same_v<IndexTag2, DimRes1_cov>) {
                //Compute the (1,2) coefficient of the inverse Jacobian matrix.
                // J^{-1}(1,2) = -J(1,2) / det(J)
                return -m_mapping.template jacobian_component<DimRes0, DimArg1_cov>(coord) / jacob;
            } else if constexpr (
                    std::is_same_v<IndexTag1, DimArg1> && std::is_same_v<IndexTag2, DimRes0_cov>) {
                //Compute the (2,1) coefficient of the inverse Jacobian matrix.
                // J^{-1}(2,1) = -J(2,1) / det(J)
                return -m_mapping.template jacobian_component<DimRes1, DimArg0_cov>(coord) / jacob;
            } else {
                //Compute the (2,1) coefficient of the inverse Jacobian matrix.
                // J^{-1}(2,2) = J(1,1) / det(J)
                return m_mapping.template jacobian_component<DimRes0, DimArg0_cov>(coord) / jacob;
            }
        }
    }
};
```


