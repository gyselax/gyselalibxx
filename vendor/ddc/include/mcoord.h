#pragma once

#include <cstddef>
#include <utility>

#include "taggedarray.h"

#include <experimental/mdspan>

using MCoordElement = std::size_t;

using MLengthElement = std::ptrdiff_t;

template <class... Tags>
using MCoord = TaggedArray<MCoordElement, Tags...>;


namespace detail {

template <class, class>
struct ExtentToMCoord;

template <class MCoordType, std::size_t... Ints>
struct ExtentToMCoord<MCoordType, std::index_sequence<Ints...>>
{
    static_assert(MCoordType::size() == sizeof...(Ints));

    template <std::size_t... Extents>
    static inline constexpr MCoordType mcoord(
            const std::experimental::extents<Extents...>& extent) noexcept
    {
        return MCoordType(extent.extent(Ints)...);
    }
};

} // namespace detail

template <class MCoordType>
struct ExtentToMCoordEnd
    : detail::ExtentToMCoord<MCoordType, std::make_index_sequence<MCoordType::size()>>
{
};

template <class... Tags, class Extents>
inline constexpr MCoord<Tags...> mcoord_end(Extents extent) noexcept
{
    return ExtentToMCoordEnd<MCoord<Tags...>>::mcoord(extent);
}
