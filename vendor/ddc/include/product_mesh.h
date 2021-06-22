#pragma once

#include <type_traits>

#include "mcoord.h"
#include "rcoord.h"
#include "single_mesh.h"
#include "taggedtuple.h"

template <class Mesh>
struct SubmeshImpl;

template <class... Meshes>
class ProductMesh
{
    template <class Mesh>
    using tag_t = typename Mesh::tag_type;

private:
    static_assert(sizeof...(Meshes) > 0, "At least 1 mesh must be provided");

    static_assert((... && (Meshes::rank() <= 1)), "Only meshes of rank <= 1 are allowed");

    TaggedTuple<detail::TypeSeq<Meshes...>, detail::TypeSeq<tag_t<Meshes>...>> m_meshes;

public:
    using rcoord_type = RCoord<tag_t<Meshes>...>;

    using mcoord_type = MCoord<tag_t<Meshes>...>;

public:
    static constexpr std::size_t rank() noexcept
    {
        return (0 + ... + Meshes::rank());
    }

public:
    ProductMesh() = delete;

    constexpr ProductMesh(Meshes const&... meshes) : m_meshes(meshes...) {}

    constexpr ProductMesh(Meshes&&... meshes) : m_meshes(std::move(meshes)...) {}

    template <class... OMeshes>
    constexpr ProductMesh(ProductMesh<OMeshes...> const& mesh)
        : m_meshes(mesh.template get<tag_t<Meshes>>()...)
    {
    }

    template <class... OMeshes>
    constexpr ProductMesh(ProductMesh<OMeshes...>&& mesh)
        : m_meshes(std::move(mesh.template get<tag_t<Meshes>>())...)
    {
    }

    ProductMesh(ProductMesh const& x) = default;

    ProductMesh(ProductMesh&& x) = default;

    ~ProductMesh() = default;

    ProductMesh& operator=(ProductMesh const& x) = default;

    ProductMesh& operator=(ProductMesh&& x) = default;

    template <class Tag>
    auto const& get() const noexcept
    {
        return ::get<Tag>(m_meshes);
    }

    template <class Tag>
    auto& get() noexcept
    {
        return ::get<Tag>(m_meshes);
    }

    template <class... QueryTags>
    RCoord<QueryTags...> to_real(MCoord<QueryTags...> const& mcoord) const noexcept
    {
        return RCoord<QueryTags...>(
                ::get<QueryTags>(m_meshes).to_real(::get<QueryTags>(mcoord))...);
    }

    template <class... Slicespecs>
    auto submesh(Slicespecs&&... slicespecs) const noexcept
    {
        return SubmeshImpl<ProductMesh>::submesh(*this, std::forward<Slicespecs>(slicespecs)...);
    }

    friend constexpr bool operator==(ProductMesh const& lhs, ProductMesh const& rhs)
    {
        return (... && (::get<tag_t<Meshes>>(lhs.m_meshes) == ::get<tag_t<Meshes>>(rhs.m_meshes)));
    }

    friend constexpr bool operator!=(ProductMesh const& lhs, ProductMesh const& rhs)
    {
        return !(lhs == rhs);
    }
};

template <class... Meshes>
struct SubmeshImpl<ProductMesh<Meshes...>>
{
    template <class Mesh, class Slicespec>
    static auto submesh_rank_1(Mesh const& mesh, Slicespec&& slicespec)
    {
        static_assert(Mesh::rank() <= 1);
        using slicespec_type = std::remove_cv_t<std::remove_reference_t<Slicespec>>;
        if constexpr (std::is_same_v<slicespec_type, std::experimental::all_type>) {
            return mesh;
        } else if constexpr (std::is_integral_v<slicespec_type>) {
            return SingleMesh(mesh.to_real(slicespec));
        }
    }

    template <class... Slicespecs>
    static auto submesh(ProductMesh<Meshes...> const& mesh, Slicespecs&&... slicespecs)
    {
        static_assert(sizeof...(Meshes) == sizeof...(Slicespecs));
        return ProductMesh(submesh_rank_1(
                mesh.template get<typename Meshes::tag_type>(),
                std::forward<Slicespecs>(slicespecs))...);
    }
};
