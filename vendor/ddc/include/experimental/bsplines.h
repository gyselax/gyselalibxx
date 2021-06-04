#pragma once

#include <cstdint>
#include <type_traits>

namespace experimental {

template <class Domain, std::size_t D>
class BSplines;

// Helper for template deduction
struct BSplinesHelper
{
    template <class Domain, std::size_t D>
    constexpr BSplines<Domain, D> operator()(
            const Domain& domain,
            std::integral_constant<std::size_t, D>) const noexcept
    {
        return BSplines<Domain, D>(domain);
    }
};

inline constexpr BSplinesHelper bsplines_helper;

} // namespace experimental
