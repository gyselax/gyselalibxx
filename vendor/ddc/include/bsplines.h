#pragma once

#include <cstdint>
#include <type_traits>

template <class Domain, std::size_t D>
class BSplines;

template <class Domain, std::size_t D>
constexpr BSplines<Domain, D> make_bsplines(
        const Domain& domain,
        std::integral_constant<std::size_t, D>) noexcept
{
    return BSplines<Domain, D>(domain);
}
