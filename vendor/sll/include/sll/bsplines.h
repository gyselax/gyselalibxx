#pragma once

#include <cstdint>
#include <type_traits>

template <class Mesh, std::size_t D>
class BSplines;

template <class Domain, std::size_t D>
constexpr BSplines<typename Domain::mesh_type, D> make_bsplines_new(
        Domain const& domain,
        std::integral_constant<std::size_t, D>) noexcept
{
    return BSplines<typename Domain::mesh_type, D>(domain);
}
