// SPDX-License-Identifier: MIT
#pragma once

#include "inv_jacobian_o_point.hpp"
#include "inverse_jacobian_matrix.hpp"
#include "mapping_tools.hpp"
#include "math_tools.hpp"

/**
 * @brief A class which describes a mapping which is constructed by
 * combining two mappings.
 * Let us denote Mapping1 as @f$ \mathcal{F} @f$ and Mapping2 as @f$ \mathcal{G} @f$
 * then this mapping represents:
 * @f$ \mathcal{F} \circ \mathcal{G} @f$
 *
 * There are therefore 3 domains in this calculation @f$ \Omega_{start} @f$, @f$ \Omega_{mid} @f$
 * and @f$ \Omega_{end} @f$ with @f$ \mathcal{F}(\Omega_{mid})\rightarrow\Omega_{end} @f$
 * and @f$ \mathcal{G}(\Omega_{start})\rightarrow\Omega_{mid} @f$
 *
 * The functions in this mapping are defined on the coordinate system associated
 * with the domain @f$ \Omega_{mid} @f$.
 */
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
    /// The type of the argument of the function described by this mapping.
    using CoordArg = typename Mapping2::CoordArg;
    /// The type of the result of the function described by this mapping.
    using CoordResult = typename Mapping1::CoordResult;
    /// The coordinate system on which the Jacobian is described.
    using CoordJacobian = typename Mapping2::CoordResult;
    /// The type of the Jacobian matrix.
    using JacobianMatrixType = DTensor<
            ddc::to_type_seq_t<CoordResult>,
            vector_index_set_dual_t<ddc::to_type_seq_t<CoordArg>>>;
    /// The type of the inverse Jacobian matrix.
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

    static_assert(has_2d_jacobian_v<Mapping1, CoordJacobian>);
    static_assert(has_2d_jacobian_v<InverseMapping2, CoordJacobian>);

private:
    Mapping1 m_mapping_1;
    Mapping2 m_mapping_2;
    InverseMapping2 m_inv_mapping_2;
    // The Jacobian defined on CoordJacobian is the inverse of the inverse mapping
    InverseJacobianMatrix<InverseMapping2, CoordJacobian> m_jacobian_mapping_2;
    double m_epsilon;

public:
    /**
     * @brief Build a CombinedMapping from the component mappings.
     * This constructor should be used if the inverse jacobian of the first mapping, or the
     * jacobian of the second mapping cannot be evaluated at the O-point.
     * @param[in] mapping_1 The first mapping.
     * @param[in] mapping_2 The second mapping.
     * @param[in] epsilon The parameter @f$ \varepsilon @f$ which determines when a point is
     *          close enough to the central O-point for linearisation to be required when
     *          calculating the inverse of the Jacobian. The Jacobian is linearised on
     *          @f$ r \in [0, \varepsilon] @f$.
     */
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

    /**
     * @brief Build a CombinedMapping from the component mappings.
     * This constructor should be used if both mappings can be safely evaluated at all points in
     * space.
     * @param[in] mapping_1 The first mapping.
     * @param[in] mapping_2 The second mapping.
     */
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

    /**
     * @brief Convert the argument coordinate to the equivalent result coordinate.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
    CoordResult operator()(CoordArg coord)
    {
        return m_mapping_1(m_mapping_2(coord));
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * For some computations, we need the complete Jacobian matrix or just the
     * coefficients.
     * This is calculated by combining the Jacobian matrices of the 2 curvilinear mappings.
     *
     * @f$ J_{(x_{in}, y_{in})->(x_{out},y_{out})} = J_{F\circ G} @f$
     * @f$ = J_{(x_{mid}, y_{mid})->(x_{out},y_{out})}J_{(x_{in}, y_{in})->(x_{mid}, y_{mid})} = J_F J_G @f$
     * @f$ = J_{(x_{mid}, y_{mid})->(x_{out},y_{out})} [J_{(x_{mid}, y_{mid})->(x_{in}, y_{in})}]^{-1} = J_F [J_{G^{-1}}]^{-1} @f$
     *
     * The Jacobians that are used for the calculation must be mappings from @f$ (x_{mid}, y_{mid}) @f$
     * so they can be calculated on the correct coordinate system.
     *
     * @param[in] coord The coordinate where we evaluate the Jacobian matrix.
     * @return The calculated Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION JacobianMatrixType jacobian_matrix(CoordJacobian const& coord) const
    {
        return tensor_mul(
                index<'i', 'j'>(m_mapping_1.jacobian_matrix(coord)),
                index<'j', 'k'>(m_jacobian_mapping_2(coord)));
    }

    /**
     * @brief Compute the (i,j) component of the Jacobian matrix.
     *
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @return The (i,j) component of the Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(CoordJacobian const& coord_rtheta) const
    {
        static_assert(
                std::is_same_v<IndexTag1, DimResult1> && std::is_same_v<IndexTag1, DimResult2>);
        static_assert((std::is_same_v<IndexTag2, typename DimArg1::Dual>)&&(
                std::is_same_v<IndexTag2, typename DimArg2::Dual>));
        //VG// TODO replace by static_assert(ddc::in_tags_v<IndexTag1, ResultTags>);
        JacobianMatrixType J = jacobian_matrix(coord_rtheta);
        return ddcHelper::get<IndexTag1, IndexTag2>(J);
    }

    /**
     * @brief Compute the determinant of the Jacobian matrix.
     * @see jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @returns The determinant of the Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION double jacobian(CoordJacobian const& coord_rtheta) const
    {
        return determinant(jacobian_matrix(coord_rtheta));
    }

    /**
     * @brief Compute the full inverse Jacobian matrix.
     *
     * If one of the mappings is singular then this function linearises the inverse Jacobian
     * between the O-point and @f$ r = \varepsilon @f$.
     * The inverse Jacobian at the O-point is calculated using the class InvJacobianOPoint in
     * this case.
     * The inverse Jacobian at a point where the matrices are not singular is calculated
     * using non_singular_inverse_jacobian_matrix
     *
     * @param[in] coord The coordinate where we evaluate the inverse Jacobian matrix.
     * @return The calculated inverse Jacobian matrix.
     *
     * @see non_singular_inverse_jacobian_matrix
     */
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

    /**
     * @brief Compute the (i,j) coefficient of the inverse Jacobian matrix.
     * @see inv_jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @return The (i,j) coefficient of the Jacobian matrix.
     */
    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double inv_jacobian_component(CoordJacobian const& coord_rtheta) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<DimArg1, DimArg2>>);
        static_assert(ddc::in_tags_v<
                      IndexTag2,
                      VectorIndexSet<typename DimResult1::Dual, typename DimResult2::Dual>>);

        InvJacobianMatrixType J = inv_jacobian_matrix(coord_rtheta);
        return ddcHelper::get<IndexTag1, IndexTag2>(J);
    }


    /**
     * @brief Compute the determinant of the Jacobian matrix.
     * @see inv_jacobian_matrix
     * @param[in] coord_rtheta The coordinate where we evaluate the Jacobian matrix.
     * @returns The determinant of the Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION double inv_jacobian(CoordJacobian const& coord_rtheta) const
    {
        return determinant(inv_jacobian_matrix(coord_rtheta));
    }

    /**
     * @brief Get access to one of the internal mappings.
     * @tparam Mapping The mapping that we want to access.
     * @return A constant reference to the internal mapping.
     */
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
    /**
     * Compute the inverse Jacobian matrix at a point where the matrices are not singular.
     *
     * This is calculated by combining the Jacobian matrices of the 2 curvilinear mappings.
     *
     * @f$ [J_{(x_{in}, y_{in})->(x_{out},y_{out})}]^{-1} = [J_{F\circ G}]^{-1} @f$
     * @f$ = [J_{(x_{mid}, y_{mid})->(x_{out},y_{out})}J_{(x_{in}, y_{in})->(x_{mid}, y_{mid})}]^{-1} = [J_F J_G]^{-1} @f$
     * @f$ = J_{(x_{in}, y_{in})->(x_{mid}, y_{mid})}^{-1}J_{(x_{mid}, y_{mid})->(x_{out},y_{out})}^{-1} = J_G^{-1} J_F^{-1} @f$
     * @f$ = J_{(x_{mid}, y_{mid})->(x_{in}, y_{in})} J_{(x_{mid}, y_{mid})->(x_{out},y_{out})}^{-1} = J_{G^{-1}} J_F^{-1} @f$
     *
     * The Jacobians that are used for the calculation must be mappings from @f$ (x_{mid}, y_{mid}) @f$
     * so they can be calculated on the correct coordinate system.
     *
     * @param[in] coord The coordinate where we evaluate the Jacobian matrix.
     */
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
