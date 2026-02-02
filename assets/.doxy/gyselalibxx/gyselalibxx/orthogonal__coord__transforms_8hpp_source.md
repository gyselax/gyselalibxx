

# File orthogonal\_coord\_transforms.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**orthogonal\_coord\_transforms.hpp**](orthogonal__coord__transforms_8hpp.md)

[Go to the documentation of this file](orthogonal__coord__transforms_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_helper.hpp"
#include "tensor.hpp"
#include "type_seq_tools.hpp"

namespace detail {

template <class CoordArg, class TypeSeqTransforms>
struct FindTransform;

template <class CoordArg, class HeadTransform, class... Transforms>
struct FindTransform<CoordArg, std::tuple<HeadTransform, Transforms...>>
{
    using type = std::conditional_t<
            ((ddc::detail::is_tagged_vector_v<CoordArg>)&&(
                    std::is_same_v<typename HeadTransform::CoordArg, CoordArg>))
                    || (ddc::in_tags_v<
                            CoordArg,
                            ddc::to_type_seq_t<typename HeadTransform::CoordArg>>),
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

template <class ArgCoord, class ResultCoord, class JacobianCoord, class... CoordTransform>
class OrthogonalCoordTransforms
{
    static_assert(sizeof...(CoordTransform) > 1);

public:
    using CoordArg = ArgCoord;
    using CoordResult = ResultCoord;
    using CoordJacobian = JacobianCoord;
    static_assert(ddc::type_seq_same_v<
                  type_seq_cat_t<ddc::to_type_seq_t<typename CoordTransform::CoordJacobian>...>,
                  ddc::to_type_seq_t<JacobianCoord>>);

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

    template <class CoordType, std::enable_if_t<!std::is_same_v<CoordType, CoordArg>, bool> = true>
    KOKKOS_INLINE_FUNCTION auto operator()(CoordType const& coord) const
    {
        using SelectedCoordTransform =
                typename detail::FindTransform<CoordType, std::tuple<CoordTransform...>>::type;
        return std::get<SelectedCoordTransform>(m_transforms)(coord);
    }

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

    template <class IndexTag1, class IndexTag2>
    KOKKOS_FUNCTION double jacobian_component(CoordArg const& coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, ddc::to_type_seq_t<CoordResult>>);
        static_assert(
                ddc::in_tags_v<IndexTag2, get_covariant_dims_t<ddc::to_type_seq_t<CoordArg>>>);
        using RelevantMapping =
                typename detail::FindTransform<IndexTag2, std::tuple<CoordTransform...>>::type;
        if constexpr (ddc::in_tags_v<
                              IndexTag1,
                              ddc::to_type_seq_t<typename RelevantMapping::CoordResult>>) {
            return std::get<RelevantMapping>(m_transforms)
                    .template jacobian_component<IndexTag1, IndexTag2>(
                            typename RelevantMapping::CoordJacobian(coord));
        } else {
            return 0;
        }
    }

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

    KOKKOS_INLINE_FUNCTION auto get_inverse_mapping() const
    {
        return OrthogonalCoordTransforms<
                ResultCoord,
                ArgCoord,
                ddcHelper::to_coord_t<type_seq_cat_t<ddc::to_type_seq_t<
                        typename inverse_mapping_t<CoordTransform>::CoordJacobian>...>>,
                inverse_mapping_t<CoordTransform>...>(
                (std::get<CoordTransform>(m_transforms).get_inverse_mapping())...);
    }
};
```


