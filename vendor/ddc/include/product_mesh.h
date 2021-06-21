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
private:
    static_assert((... && (Meshes::rank() <= 1)), "Only meshes of rank <= 1 are allowed");

    TaggedTuple<detail::TypeSeq<Meshes...>, detail::TypeSeq<typename Meshes::tag_type...>> m_meshes;

public:
    using rcoord_type = RCoord<typename Meshes::tag_type...>;

    using mcoord_type = MCoord<typename Meshes::tag_type...>;

public:
    static constexpr std::size_t rank() noexcept
    {
        return (0 + ... + Meshes::rank());
    }

public:
    constexpr ProductMesh(Meshes const&... meshes) : m_meshes(meshes...) {}

    constexpr ProductMesh(Meshes&&... meshes) : m_meshes(std::move(meshes)...) {}

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
