// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_helper.hpp"
#include "tensor.hpp"
#include "type_seq_tools.hpp"

namespace detail {

template <class CoordArg, class TypeSeqGrid>
struct FindTransform;

template <class CoordArg, class HeadTransform, class... Transforms>
struct FindTransform<CoordArg, std::tuple<HeadTransform, Transforms...>>
{
    using type = std::conditional_t<
            std::is_same_v<typename HeadTransform::CoordArg, CoordArg>,
            HeadTransform,
            typename FindTransform<CoordArg, std::tuple<Transforms...>>::type>;
};

template <class CoordArg>
struct FindTransform<CoordArg, std::tuple<>>
{
    static_assert(std::is_same_v<CoordArg, CoordArg>, "Transform not found");
    using type = void;
};

template <class ContraDim, class CovDimsSet, class TensorType, class CoordTransformType>
struct FillJacobianSubsetRow;

template <class ContraDim, class... CovDims, class TensorType, class CoordTransformType>
struct FillJacobianSubsetRow<ContraDim, VectorIndexSet<CovDims...>, TensorType, CoordTransformType>
{
    static void fill(
            TensorType& mat,
            CoordTransformType const& transform,
            typename CoordTransformType::CoordArg const& coord)
    {
        typename CoordTransformType::CoordArg coord_arg(coord);
        ((ddcHelper::get<ContraDim, CovDims>(mat)
          = transform.template jacobian_component<ContraDim, CovDims>(coord_arg)),
         ...);
    }
};

template <class ContraDimsSet, class CovDimsSet, class TensorType, class CoordTransformType>
struct FillJacobianSubset;
template <class... ContraDims, class CovDimsSet, class TensorType, class CoordTransformType>
struct FillJacobianSubset<VectorIndexSet<ContraDims...>, CovDimsSet, TensorType, CoordTransformType>
{
    static void fill(
            TensorType& mat,
            CoordTransformType const& transform,
            typename CoordTransformType::CoordArg const& coord)
    {
        ((FillJacobianSubsetRow<ContraDims, CovDimsSet, TensorType, CoordTransformType>::
                  fill(mat, transform, coord)),
         ...);
    }
};

} // namespace detail

template <class ArgCoord, class ResultCoord, class... CoordTransform>
class OrthogonalCoordTransforms
{
    static_assert(sizeof...(CoordTransform) > 1);

public:
    using CoordArg = ArgCoord;
    using CoordResult = ResultCoord;
    using CoordJacobian = ddcHelper::to_coord_t<
            type_seq_cat_t<ddc::to_type_seq_t<typename CoordTransform::CoordJacobian>...>>;

private:
    std::tuple<CoordTransform...> m_transforms;

public:
    explicit KOKKOS_FUNCTION OrthogonalCoordTransforms(CoordTransform const&... transform)
        : m_transforms(transform...)
    {
    }

    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        return CoordResult(std::get<CoordTransform>(m_transforms)(
                typename CoordTransform::CoordArg(coord))...);
    }

    template <class CoordType>
    KOKKOS_INLINE_FUNCTION auto operator()(CoordType const& coord) const
    {
        using SelectedCoordTransform =
                typename detail::FindTransform<CoordType, std::tuple<CoordTransform...>>::type;
        return std::get<SelectedCoordTransform>(m_transforms)(coord);
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
    KOKKOS_FUNCTION DTensor<
            ddc::to_type_seq_t<CoordResult>,
            get_covariant_dims_t<ddc::to_type_seq_t<CoordArg>>>
    jacobian_matrix(CoordArg const& coord) const
    {
        using TensorType = DTensor<
                ddc::to_type_seq_t<CoordResult>,
                get_covariant_dims_t<ddc::to_type_seq_t<CoordArg>>>;
        TensorType jacobian_matrix;
        ((detail::FillJacobianSubset<
                 ddc::to_type_seq_t<typename CoordTransform::CoordResult>,
                 get_covariant_dims_t<ddc::to_type_seq_t<typename CoordTransform::CoordArg>>,
                 TensorType,
                 CoordTransform>::
                  fill(jacobian_matrix,
                       std::get<CoordTransform>(m_transforms),
                       typename CoordTransform::CoordArg(coord))),
         ...);
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
        static_assert(ddc::in_tags_v<IndexTag1, ddc::to_type_seq_t<CoordResult>>);
        static_assert(
                ddc::in_tags_v<IndexTag2, get_covariant_dims_t<ddc::to_type_seq_t<CoordArg>>>);
        double result
                = ((((ddc::in_tags_v<IndexTag1, ddc::to_type_seq_t<typename CoordTransform::CoordResult>>)&&(
                            ddc::in_tags_v<
                                    IndexTag2,
                                    get_covariant_dims_t<ddc::to_type_seq_t<
                                            typename CoordTransform::CoordArg>>>))
                            ? std::get<CoordTransform>(m_transforms)
                                      .template jacobian_component<IndexTag1, IndexTag2>(
                                              typename CoordTransform::CoordArg(coord))
                            : 0.0)
                   + ...);
        return result;
    }

    /**
     * @brief Compute the Jacobian, the determinant of the Jacobian matrix of the mapping.
     *
     * @param[in] coord
     *          The coordinate where we evaluate the Jacobian.
     *
     * @return A double with the value of the determinant of the Jacobian matrix.
     */
    KOKKOS_FUNCTION double jacobian(CoordJacobian const& coord) const
    {
        double result
                = ((std::get<CoordTransform>(m_transforms)
                            .jacobian(typename CoordTransform::CoordJacobian(coord)))
                   * ...);
        result *= type_seq_permutation_parity_v<
                ddc::to_type_seq_t<CoordResult>,
                type_seq_cat_t<ddc::to_type_seq_t<typename CoordTransform::CoordResult>...>>;
        return result;
    }

    /**
     * @brief Get the inverse mapping.
     *
     * @return The inverse mapping.
     */
    KOKKOS_INLINE_FUNCTION auto get_inverse_mapping() const
    {
        return OrthogonalCoordTransforms<
                ResultCoord,
                ArgCoord,
                inverse_mapping_t<CoordTransform>...>(
                (std::get<CoordTransform>(m_transforms).get_inverse_mapping())...);
    }
};
