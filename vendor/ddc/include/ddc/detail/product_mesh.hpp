#pragma once

#include <tuple>
#include <type_traits>

#include "ddc/mcoord.hpp"
#include "ddc/mesh.hpp"
#include "ddc/rcoord.hpp"

namespace detail {
template <class... Meshes>
class ProductMesh
{
    template <class Mesh>
    using rdim_t = typename Mesh::rdim_type;

    template <class Mesh>
    using storage_t = Mesh const&;

private:
    // static_assert((... && is_mesh_v<Meshes>), "A template parameter is not a mesh");

    static_assert(sizeof...(Meshes) > 0, "At least 1 mesh must be provided");

    static_assert((... && (Meshes::rank() <= 1)), "Only meshes of rank <= 1 are allowed");

    std::tuple<storage_t<Meshes>...> m_meshes;

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

    ProductMesh(ProductMesh const& x) = default;

    ProductMesh(ProductMesh&& x) = default;

    ~ProductMesh() = default;

    ProductMesh& operator=(ProductMesh const& x) = default;

    ProductMesh& operator=(ProductMesh&& x) = default;

    template <class QueryMesh>
    QueryMesh const& get() const noexcept
    {
        return std::get<storage_t<QueryMesh>>(m_meshes);
    }

    template <class... QueryMeshes>
    RCoord<rdim_t<QueryMeshes>...> to_real(MCoord<QueryMeshes...> const& mcoord) const noexcept
    {
        return RCoord<rdim_t<QueryMeshes>...>(
                std::get<storage_t<QueryMeshes>>(m_meshes).to_real(select<QueryMeshes>(mcoord))...);
    }

    friend constexpr bool operator==(ProductMesh const& lhs, ProductMesh const& rhs)
    {
        return (...
                && (std::get<storage_t<Meshes>>(lhs.m_meshes)
                    == std::get<storage_t<Meshes>>(rhs.m_meshes)));
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
constexpr QueryMesh const& get(ProductMesh<Meshes...> const& mesh)
{
    return mesh.template get<QueryMesh>();
}

template <class... QueryMeshes, class... Meshes>
constexpr ProductMesh<QueryMeshes...> select(ProductMesh<Meshes...> const& mesh)
{
    return ProductMesh(get<QueryMeshes>(mesh)...);
}

}
