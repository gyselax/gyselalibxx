// SPDX-License-Identifier: MIT

#pragma once

#include <utility>

#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/for_each.hpp"

namespace detail {

/** A serial reduction over a nD domain
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] init the initial value of the generalized sum
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and init. 
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <
        class... DDims,
        class T,
        class BinaryReductionOp,
        class UnaryTransformOp,
        class... DCoords>
inline T transform_reduce_serial(
        DiscreteDomain<DDims...> const& domain,
        T init,
        BinaryReductionOp const& reduce,
        UnaryTransformOp const& transform,
        DCoords const&... dcoords) noexcept
{
    if constexpr (sizeof...(DCoords) == sizeof...(DDims)) {
        return reduce(init, transform(DiscreteCoordinate<DDims...>(dcoords...)));
    } else {
        using CurrentDDim = type_seq_element_t<sizeof...(DCoords), detail::TypeSeq<DDims...>>;
        for (DiscreteCoordinate<CurrentDDim> const ii : select<CurrentDDim>(domain)) {
            init = transform_reduce_serial(domain, init, reduce, transform, dcoords..., ii);
        }
        return init;
    }
}

} // namespace detail

/** A reduction over a nD domain using the default execution policy
 * @param[in] policy the execution policy to use
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] init the initial value of the generalized sum
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and init.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class... DDims, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce(
        [[maybe_unused]] serial_policy policy,
        DiscreteDomain<DDims...> const& domain,
        T init,
        BinaryReductionOp&& reduce,
        UnaryTransformOp&& transform) noexcept
{
    return detail::transform_reduce_serial(
            domain,
            init,
            std::forward<BinaryReductionOp>(reduce),
            std::forward<UnaryTransformOp>(transform));
}

/** A reduction over a nD domain using the default execution policy
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] init the initial value of the generalized sum
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and init.
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class... DDims, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce(
        DiscreteDomain<DDims...> const& domain,
        T init,
        BinaryReductionOp&& reduce,
        UnaryTransformOp&& transform) noexcept
{
    return transform_reduce(
            default_policy(),
            domain,
            init,
            std::forward<BinaryReductionOp>(reduce),
            std::forward<UnaryTransformOp>(transform));
}
