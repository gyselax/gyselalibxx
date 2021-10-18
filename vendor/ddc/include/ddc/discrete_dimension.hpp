#pragma once

#include <type_traits>

/** `DiscreteDimension` should be used as a base class of tag classes for discrete dimensions.
 *
 * It should implement the method `to_real` if this is a discretization of a continuous (real)
 * dimension.
 */
struct DiscreteDimension
{
    void operator=(DiscreteDimension const&) = delete;
    void operator=(DiscreteDimension&&) = delete;
};

template <class T>
inline constexpr bool is_discrete_dimension_v = std::is_base_of_v<DiscreteDimension, T>;
