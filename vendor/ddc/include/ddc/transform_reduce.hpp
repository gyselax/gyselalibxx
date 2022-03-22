// SPDX-License-Identifier: MIT

#pragma once

#include <utility>

#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_domain.hpp"
#include "ddc/for_each.hpp"

/** A sequential reduction over a 1D domain
 * @param[in] policy the execution policy to use
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] init the initial value of the generalized sum
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and init. 
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class DDim, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce(
        serial_policy policy,
        DiscreteDomain<DDim> const& domain,
        T init,
        BinaryReductionOp&& reduce,
        UnaryTransformOp&& transform) noexcept
{
    DiscreteDomainIterator<DDim> const it_b = domain.begin();
    DiscreteDomainIterator<DDim> const it_e = domain.end();
    for (DiscreteDomainIterator<DDim> it = it_b; it != it_e; ++it) {
        init = reduce(init, transform(*it));
    }
    return init;
}

/** A sequential reduction over a 2D domain
 * @param[in] policy the execution policy to use
 * @param[in] domain the range over which to apply the algorithm
 * @param[in] init the initial value of the generalized sum
 * @param[in] reduce a binary FunctionObject that will be applied in unspecified order to the
 *            results of transform, the results of other reduce and init. 
 * @param[in] transform a unary FunctionObject that will be applied to each element of the input
 *            range. The return type must be acceptable as input to reduce
 */
template <class DDim1, class DDim2, class T, class BinaryReductionOp, class UnaryTransformOp>
inline T transform_reduce(
        [[maybe_unused]] serial_policy policy,
        DiscreteDomain<DDim1, DDim2> const& domain,
        T init,
        BinaryReductionOp&& reduce,
        UnaryTransformOp&& transform) noexcept
{
    DiscreteDomainIterator<DDim1> const it1_b = select<DDim1>(domain).begin();
    DiscreteDomainIterator<DDim1> const it1_e = select<DDim1>(domain).end();
    DiscreteDomainIterator<DDim2> const it2_b = select<DDim2>(domain).begin();
    DiscreteDomainIterator<DDim2> const it2_e = select<DDim2>(domain).end();
    for (DiscreteDomainIterator<DDim1> it1 = it1_b; it1 != it1_e; ++it1) {
        for (DiscreteDomainIterator<DDim2> it2 = it2_b; it2 != it2_e; ++it2) {
            init = reduce(init, transform(DiscreteCoordinate<DDim1, DDim2>(*it1, *it2)));
        }
    }
    return init;
}

/** A reduction over a 2D domain using the default execution policy
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
            serial_policy(),
            domain,
            init,
            std::forward<BinaryReductionOp>(reduce),
            std::forward<UnaryTransformOp>(transform));
}
