#pragma once

#include <type_traits>

#include "ddc/mcoord.h"
#include "ddc/mesh.h"
#include "ddc/rcoord.h"
#include "ddc/taggedtuple.h"

template <class... Meshes>
class ProductMesh
{
    template <class Mesh>
    using rdim_t = typename Mesh::rdim_type;

private:
    static_assert((... && is_mesh_v<Meshes>), "A template parameter is not a mesh");

    static_assert(sizeof...(Meshes) > 0, "At least 1 mesh must be provided");

    static_assert((... && (Meshes::rank() <= 1)), "Only meshes of rank <= 1 are allowed");

    TaggedTuple<detail::TypeSeq<Meshes const&...>, detail::TypeSeq<Meshes...>> m_meshes;

public:
    using rcoord_type = RCoord<rdim_t<Meshes>...>;

    using mcoord_type = MCoord<Meshes...>;

public:
    static constexpr std::size_t rank() noexcept
    {
        return (0 + ... + Meshes::rank());
    }

public:
    ProductMesh() = default;

    constexpr explicit ProductMesh(Meshes const&... meshes) : m_meshes(meshes...) {}

    // template <class... OMeshes>
    // constexpr ProductMesh(ProductMesh<OMeshes...> const& mesh)
    //     : m_meshes(mesh.template get<Meshes>()...)
    // {
    // }

    // template <class... OMeshes>
    // constexpr ProductMesh(ProductMesh<OMeshes...>&& mesh)
    //     : m_meshes(std::move(mesh.template get<Meshes>())...)
    // {
    // }

    ProductMesh(ProductMesh const& x) = default;

    ProductMesh(ProductMesh&& x) = default;

    ~ProductMesh() = default;

    ProductMesh& operator=(ProductMesh const& x) = default;

    ProductMesh& operator=(ProductMesh&& x) = default;

    template <
            std::size_t N = sizeof...(Meshes),
            class Mesh0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<Meshes...>>>>
    constexpr operator Mesh0 const &() const
    {
        return ::get<Mesh0>(m_meshes);
    }

    template <
            std::size_t N = sizeof...(Meshes),
            class Mesh0 = std::enable_if_t<N == 1, std::tuple_element_t<0, std::tuple<Meshes...>>>>
    constexpr operator Mesh0&()
    {
        return ::get<Mesh0>(m_meshes);
    }

    template <class Mesh>
    auto const& get() const noexcept
    {
        return ::get<Mesh>(m_meshes);
    }

    template <class Mesh>
    auto& get() noexcept
    {
        return ::get<Mesh>(m_meshes);
    }

    template <class Mesh>
    friend auto const& get(ProductMesh const& pmesh) noexcept
    {
        return pmesh.get<Mesh>();
    }

    template <class Mesh>
    friend auto& get(ProductMesh& pmesh) noexcept
    {
        return pmesh.get<Mesh>();
    }

    template <class... QueryMeshes>
    RCoord<rdim_t<QueryMeshes>...> to_real(MCoord<QueryMeshes...> const& mcoord) const noexcept
    {
        return RCoord<rdim_t<QueryMeshes>...>(
                ::get<QueryMeshes>(m_meshes).to_real(::get<QueryMeshes>(mcoord))...);
    }

    friend constexpr bool operator==(ProductMesh const& lhs, ProductMesh const& rhs)
    {
        return (... && (::get<Meshes>(lhs.m_meshes) == ::get<Meshes>(rhs.m_meshes)));
    }

#if __cplusplus <= 201703L
    // Shall not be necessary anymore in C++20
    // `a!=b` shall be translated by the compiler to `!(a==b)`
    friend constexpr bool operator!=(ProductMesh const& lhs, ProductMesh const& rhs)
    {
        return !(lhs == rhs);
    }
#endif
};


template <class QueryMesh, class... Meshes>
constexpr auto const& get(ProductMesh<Meshes...> const& mesh)
{
    return mesh.template get<QueryMesh>();
}

template <class QueryMesh, class... Meshes>
constexpr auto& get(ProductMesh<Meshes...>& mesh)
{
    return mesh.template get<QueryMesh>();
}

template <class... QueryMeshes, class... Meshes>
constexpr ProductMesh<QueryMeshes...> select(ProductMesh<Meshes...> const& mesh)
{
    return ProductMesh(::get<QueryMeshes>(mesh)...);
}
