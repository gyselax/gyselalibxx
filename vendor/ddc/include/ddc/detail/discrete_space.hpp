#pragma once

#include <tuple>
#include <type_traits>

#include "ddc/coordinate.hpp"
#include "ddc/discrete_coordinate.hpp"
#include "ddc/discrete_dimension.hpp"

namespace detail {
template <class... DDims>
class DiscreteSpace
{
    template <class DDim>
    using rdim_t = typename DDim::rdim_type;

    template <class DDim>
    using storage_t = DDim const*;

private:
    // static_assert((... && is_mesh_v<DDims>), "A template parameter is not a mesh");

    static_assert(sizeof...(DDims) > 0, "At least 1 mesh must be provided");

    static_assert((... && (DDims::rank() <= 1)), "Only ddims of rank <= 1 are allowed");

    std::tuple<storage_t<DDims>...> m_meshes;

public:
    using rcoord_type = Coordinate<rdim_t<DDims>...>;

    using mcoord_type = DiscreteCoordinate<DDims...>;

public:
    static constexpr std::size_t rank() noexcept
    {
        return (0 + ... + DDims::rank());
    }

public:
    DiscreteSpace() = default;

    constexpr explicit DiscreteSpace(DDims const&... ddims) : m_meshes(&ddims...) {}

    DiscreteSpace(DiscreteSpace const& x) = default;

    DiscreteSpace(DiscreteSpace&& x) = default;

    ~DiscreteSpace() = default;

    DiscreteSpace& operator=(DiscreteSpace const& x) = default;

    DiscreteSpace& operator=(DiscreteSpace&& x) = default;

    template <class QueryDDim>
    QueryDDim const& get() const noexcept
    {
        return *std::get<storage_t<QueryDDim>>(m_meshes);
    }

    template <class... QueryDDims>
    Coordinate<rdim_t<QueryDDims>...> to_real(
            DiscreteCoordinate<QueryDDims...> const& mcoord) const noexcept
    {
        return Coordinate<rdim_t<QueryDDims>...>(
                std::get<storage_t<QueryDDims>>(m_meshes)->to_real(select<QueryDDims>(mcoord))...);
    }

    friend constexpr bool operator==(DiscreteSpace const& lhs, DiscreteSpace const& rhs)
    {
        return (...
                && (*std::get<storage_t<DDims>>(lhs.m_meshes)
                    == *std::get<storage_t<DDims>>(rhs.m_meshes)));
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    friend constexpr bool operator!=(DiscreteSpace const& lhs, DiscreteSpace const& rhs)
    {
        return !(lhs == rhs);
    }
#endif
};

template <class QueryDDim, class... DDims>
constexpr QueryDDim const& get(DiscreteSpace<DDims...> const& mesh)
{
    return mesh.template get<QueryDDim>();
}

template <class... QueryDDims, class... DDims>
constexpr DiscreteSpace<QueryDDims...> select(DiscreteSpace<DDims...> const& mesh)
{
    return DiscreteSpace(get<QueryDDims>(mesh)...);
}

} // namespace detail
