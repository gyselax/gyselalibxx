// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_helper.hpp"
#include "tensor.hpp"
#include "type_seq_tools.hpp"

namespace detail {

/**
 * A structure used at compile time to find which of the transformations in a given
 * TypeSeq converts from coordinates of the type CoordArg.
 */
template <class CoordArg, class TypeSeqTransforms>
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

} // namespace detail

/**
 * @brief A multi-dimensional coordinate transformation which can be decomposed
 * into multiple orthogonal coordinate transformations.
 *
 * E.g. @f$ (x_1,x_2) = (y_1 + 3, 4*y_2 + 7) @f$
 * Here the coordinates @f$ y_1 @f$ and @f$ y_2 @f$ only appear in the description
 * of either @f$ x_1 @f$ or @f$ x_2 @f$.

 * E.g. a cylindrical transformation which can be decomposed into a 2D circular
 * transformation and a linear transformation.
 *
 * @tparam ArgCoord The type of the input coordinates.
 * @tparam ResultCoord The type of the output coordinates.
 * @tparam CoordTransform The coordinate transformations comprising this coordinate
 *                  transformation. Note that the order of these is unrelated to
 *                  the ordering chosen for the coordinates.
 */
template <class ArgCoord, class ResultCoord, class... CoordTransform>
class OrthogonalCoordTransforms
{
    static_assert(sizeof...(CoordTransform) > 1);

public:
    /// The type of the argument of the function described by this mapping
    using CoordArg = ArgCoord;
    /// The type of the result of the function described by this mapping
    using CoordResult = ResultCoord;
    /// The type of the coordinate that can be used to evaluate the Jacobian of this mapping
    using CoordJacobian = ddcHelper::to_coord_t<
            type_seq_cat_t<ddc::to_type_seq_t<typename CoordTransform::CoordJacobian>...>>;

private:
    std::tuple<CoordTransform...> m_transforms;

public:
    /**
     * @brief Construct a multi-dimensional coordinate transformation.
     *
     * @param[in] transform The coordinate transformations comprising this coordinate transformation.
     */
    explicit KOKKOS_FUNCTION OrthogonalCoordTransforms(CoordTransform const&... transform)
        : m_transforms(transform...)
    {
    }

    /**
     * @brief Convert the coordinate to the output coordinate system.
     *
     * @param[in] coord The coordinate to be converted expressed on the input coordinate system.
     *
     * @return The coordinate expressed on the output coordinate system.
     */
    KOKKOS_FUNCTION CoordResult operator()(CoordArg const& coord) const
    {
        return CoordResult(std::get<CoordTransform>(m_transforms)(
                typename CoordTransform::CoordArg(coord))...);
    }

    /**
     * @brief Convert a subset of coordinates to the corresponding subset of the output coordinate system.
     *
     * Neglected elements of the coordinate system must be orthogonal to the provided coordinates.
     *
     * @param[in] coord The coordinate to be converted expressed on the input coordinate system.
     *
     * @return The coordinate expressed on the output coordinate system.
     */
    template <class CoordType, std::enable_if_t<!std::is_same_v<CoordType, CoordArg>, bool> = true>
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
        TensorType jacobian_matrix(0);
        ((ddcHelper::assign_elements(
                 jacobian_matrix,
                 std::get<CoordTransform>(m_transforms)
                         .jacobian_matrix(typename CoordTransform::CoordJacobian(coord)))),
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
