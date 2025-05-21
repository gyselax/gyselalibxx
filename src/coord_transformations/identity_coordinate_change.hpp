// SPDX-License-Identifier: MIT
#pragma once
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"

/**
 * @brief A class describing an identity transformation.
 *
 * It is not expected that this class appear in a final simulation,
 * but it may be useful when debugging vector calculations on
 * general coordinates.
 *
 * @tparam ArgBasis A VectorIndexSet containing the continuous
 *      dimensions on which the argument is described.
 * @tparam ResultBasis A VectorIndexSet containing the continuous
 *      dimensions on which the result is described.
 */
template <class ArgBasis, class ResultBasis>
class IdentityCoordinateChange
{
public:
    /// The type of the argument of the function described by this mapping
    using CoordArg = ddcHelper::to_coord_t<ArgBasis>;
    /// The type of the result of the function described by this mapping
    using CoordResult = ddcHelper::to_coord_t<ResultBasis>;
    /// The type of the coordinate that can be used to evaluate the Jacobian of this mapping
    using CoordJacobian = ddcHelper::to_coord_t<ArgBasis>;

private:
    using ArgBasisCov = vector_index_set_dual_t<ArgBasis>;

public:
    /**
     * @brief Convert the coordinate in the argument basis to the equivalent coordinate in the result basis.
     *
     * @param[in] coord The coordinate to be converted.
     *
     * @return The equivalent coordinate.
     */
    template <class... ArgDims>
    KOKKOS_FUNCTION CoordResult operator()(Coord<ArgDims...> const& coord) const
    {
        static_assert(std::is_same_v<CoordArg, Coord<ArgDims...>>);
        return CoordResult((ddc::get<ArgDims>(coord))...);
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION double jacobian(CoordArg const& coord) const
    {
        return 1;
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<ResultBasis, ArgBasisCov> jacobian_matrix(CoordArg const& coord) const
    {
        DTensor<ResultBasis, ArgBasisCov> jacobian_matrix(0.0);
        fill_diagonal_elements(jacobian_matrix);
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
    KOKKOS_INLINE_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ResultBasis>);
        static_assert(ddc::in_tags_v<IndexTag2, ArgBasisCov>);
        if constexpr (
                (ddc::type_seq_rank_v<IndexTag1, ResultBasis>)
                == (ddc::type_seq_rank_v<IndexTag2, ArgBasisCov>)) {
            return 1;
        } else {
            return 0;
        }
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_INLINE_FUNCTION double inv_jacobian(CoordArg const& coord) const
    {
        return 1;
    }

    /**
     * @brief Compute full Jacobian matrix.
     *
     * @param[in] coord
     * 				The coordinate where we evaluate the Jacobian matrix.
     * @return The Jacobian matrix.
     */
    KOKKOS_FUNCTION DTensor<ResultBasis, ArgBasisCov> inv_jacobian_matrix(
            CoordArg const& coord) const
    {
        DTensor<ResultBasis, ArgBasisCov> inv_jacobian_matrix(0.0);
        fill_diagonal_elements(inv_jacobian_matrix);
        return inv_jacobian_matrix;
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
    KOKKOS_INLINE_FUNCTION double inv_jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ResultBasis>);
        static_assert(ddc::in_tags_v<IndexTag2, ArgBasisCov>);
        if constexpr (
                (ddc::type_seq_rank_v<IndexTag1, ResultBasis>)
                == (ddc::type_seq_rank_v<IndexTag2, ArgBasisCov>)) {
            return 1;
        } else {
            return 0;
        }
    }

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    KOKKOS_INLINE_FUNCTION IdentityCoordinateChange<ResultBasis, ArgBasis> get_inverse_mapping()
            const
    {
        return IdentityCoordinateChange<ResultBasis, ArgBasis>();
    }

private:
    template <class... RowDims, class... ColDims>
    void fill_diagonal_elements(
            DTensor<VectorIndexSet<RowDims...>, VectorIndexSet<ColDims...>>& matrix) const
    {
        ((ddcHelper::get<RowDims, ColDims>(matrix) = 1.0), ...);
    }
};

namespace mapping_detail {
template <class ArgBasis, class ResultBasis, class ExecSpace>
struct MappingAccessibility<ExecSpace, IdentityCoordinateChange<ArgBasis, ResultBasis>>
    : std::true_type
{
};
} // namespace mapping_detail
