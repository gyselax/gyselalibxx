

# File combined\_mapping.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**combined\_mapping.hpp**](combined__mapping_8hpp.md)

[Go to the documentation of this file](combined__mapping_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "inv_jacobian_o_point.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "mapping_tools.hpp"
#include "math_tools.hpp"

template <class Mapping1, class Mapping2>
class CombinedMapping
{
    static_assert(is_mapping_v<Mapping1>);
    static_assert(is_mapping_v<Mapping2>);
    static_assert(std::is_same_v<typename Mapping2::CoordResult, typename Mapping1::CoordArg>);
    static_assert(
            is_analytical_mapping_v<Mapping2>,
            "The second mapping must be analytical to evaluate the jacobian");

public:
    using CoordArg = typename Mapping2::CoordArg;
    using CoordResult = typename Mapping1::CoordResult;
    using CoordJacobian = typename Mapping2::CoordResult;
    using JacobianMatrixType = DTensor<
            ddc::to_type_seq_t<CoordResult>,
            vector_index_set_dual_t<ddc::to_type_seq_t<CoordArg>>>;
    using InvJacobianMatrixType = DTensor<
            ddc::to_type_seq_t<CoordArg>,
            vector_index_set_dual_t<ddc::to_type_seq_t<CoordResult>>>;

private:
    using DimArg1 = ddc::type_seq_element_t<0, ddc::to_type_seq_t<CoordArg>>;
    using DimArg2 = ddc::type_seq_element_t<1, ddc::to_type_seq_t<CoordArg>>;
    using DimResult1 = ddc::type_seq_element_t<0, ddc::to_type_seq_t<CoordResult>>;
    using DimResult2 = ddc::type_seq_element_t<1, ddc::to_type_seq_t<CoordResult>>;

private:
    using InverseMapping2 = inverse_mapping_t<Mapping2>;

    static_assert(has_jacobian_v<Mapping1, CoordJacobian>);
    static_assert(has_jacobian_v<InverseMapping2, CoordJacobian>);

private:
    Mapping1 m_mapping_1;
    Mapping2 m_mapping_2;
    InverseMapping2 m_inv_mapping_2;
    // The Jacobian defined on CoordJacobian is the inverse of the inverse mapping
    InverseJacobianMatrix<InverseMapping2, CoordJacobian> m_jacobian_mapping_2;
    double m_epsilon;

public:
    template <
            class Map1,
            std::enable_if_t<
                    (has_singular_o_point_inv_jacobian_v<Map1>)
                            || (has_singular_o_point_inv_jacobian_v<InverseMapping2>),
                    bool> = true>
    CombinedMapping(Map1 mapping_1, Mapping2 mapping_2, double epsilon)
        : m_mapping_1(mapping_1)
        , m_mapping_2(mapping_2)
        , m_inv_mapping_2(mapping_2.get_inverse_mapping())
        , m_jacobian_mapping_2(mapping_2.get_inverse_mapping())
        , m_epsilon(epsilon)
    {
        static_assert(std::is_same_v<Mapping1, Map1>);
    }

    template <
            class Map1,
            std::enable_if_t<
                    !((has_singular_o_point_inv_jacobian_v<Map1>)
                      || (has_singular_o_point_inv_jacobian_v<InverseMapping2>)),
                    bool> = true>
    CombinedMapping(Map1 mapping_1, Mapping2 mapping_2)
        : m_mapping_1(mapping_1)
        , m_mapping_2(mapping_2)
        , m_inv_mapping_2(mapping_2.get_inverse_mapping())
        , m_jacobian_mapping_2(mapping_2.get_inverse_mapping())
        , m_epsilon(0.0)
    {
        static_assert(std::is_same_v<Mapping1, Map1>);
    }

    CoordResult operator()(CoordArg coord)
    {
        return m_mapping_1(m_mapping_2(coord));
    }

    KOKKOS_INLINE_FUNCTION JacobianMatrixType jacobian_matrix(CoordJacobian const& coord) const
    {
        return tensor_mul(
                index<'i', 'j'>(m_mapping_1.jacobian_matrix(coord)),
                index<'j', 'k'>(m_jacobian_mapping_2(coord)));
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(CoordJacobian const& coord_rtheta) const
    {
        JacobianMatrixType J = jacobian_matrix(coord_rtheta);
        return ddcHelper::get<IndexTag1, IndexTag2>(J);
    }

    KOKKOS_INLINE_FUNCTION double jacobian(CoordJacobian const& coord_rtheta) const
    {
        return determinant(jacobian_matrix(coord_rtheta));
    }

    KOKKOS_FUNCTION InvJacobianMatrixType inv_jacobian_matrix(CoordJacobian const& coord) const
    {
        if constexpr (
                (has_singular_o_point_inv_jacobian_v<Mapping1>)
                || (has_singular_o_point_inv_jacobian_v<InverseMapping2>)) {
            using R = ddc::type_seq_element_t<0, ddc::to_type_seq_t<CoordJacobian>>;
            using Theta = ddc::type_seq_element_t<1, ddc::to_type_seq_t<CoordJacobian>>;
            double r = ddc::get<R>(coord);
            if (r < m_epsilon) {
                InvJacobianOPoint<CombinedMapping<Mapping1, Mapping2>, CoordJacobian> o_point_val(
                        *this);
                CoordJacobian coord_eps(m_epsilon, ddc::get<Theta>(coord));
                Tensor J_0 = o_point_val();
                Tensor J_eps = non_singular_inverse_jacobian_matrix(coord_eps);
                return (1 - r / m_epsilon) * J_0 + r / m_epsilon * J_eps;
            } else {
                return non_singular_inverse_jacobian_matrix(coord);
            }
        } else {
            return non_singular_inverse_jacobian_matrix(coord);
        }
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double inv_jacobian_component(CoordJacobian const& coord_rtheta) const
    {
        InvJacobianMatrixType J = inv_jacobian_matrix(coord_rtheta);
        return ddcHelper::get<IndexTag1, IndexTag2>(J);
    }


    KOKKOS_INLINE_FUNCTION double inv_jacobian(CoordJacobian const& coord_rtheta) const
    {
        return determinant(inv_jacobian_matrix(coord_rtheta));
    }

    template <class Mapping>
    KOKKOS_INLINE_FUNCTION Mapping const& get() const
    {
        static_assert(std::is_same_v<Mapping, Mapping1> || std::is_same_v<Mapping, Mapping2>);
        if constexpr (std::is_same_v<Mapping, Mapping1>) {
            return m_mapping_1;
        } else {
            return m_mapping_2;
        }
    }

private:
    KOKKOS_INLINE_FUNCTION InvJacobianMatrixType
    non_singular_inverse_jacobian_matrix(CoordJacobian const& coord) const
    {
        InverseJacobianMatrix<Mapping1, CoordJacobian> inv_jacobian_matrix_1(m_mapping_1);
        return tensor_mul(
                index<'i', 'j'>(m_inv_mapping_2.jacobian_matrix(coord)),
                index<'j', 'k'>(inv_jacobian_matrix_1(coord)));
    }
};


namespace mapping_detail {
template <class Mapping1, class Mapping2, class ExecSpace>
struct MappingAccessibility<ExecSpace, CombinedMapping<Mapping1, Mapping2>>
{
    static constexpr bool value = MappingAccessibility<ExecSpace, Mapping1>::value
                                  && MappingAccessibility<ExecSpace, Mapping2>::value;
};
} // namespace mapping_detail
```


