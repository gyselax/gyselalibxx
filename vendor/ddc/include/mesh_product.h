#pragma once

#include "mcoord.h"
#include "rcoord.h"
#include "taggedtuple.h"

template <class... Meshes>
class MeshProduct
{
private:
    static_assert((... && (Meshes::rank() == 1)), "Only meshes of rank 1 allowed");

    TaggedTuple<detail::TypeSeq<Meshes...>, detail::TypeSeq<typename Meshes::Tag_...>> m_meshes;

public:
    using RCoord_ = RCoord<typename Meshes::Tag_...>;

    using MCoord_ = MCoord<typename Meshes::Tag_...>;

public:
    static constexpr std::size_t rank() noexcept
    {
        return (0 + ... + Meshes::rank());
    }

public:
    constexpr MeshProduct(Meshes const&... meshes) : m_meshes(meshes...) {}

    constexpr MeshProduct(Meshes&&... meshes) : m_meshes(std::move(meshes)...) {}

    MeshProduct(MeshProduct const& x) = default;

    MeshProduct(MeshProduct&& x) = default;

    ~MeshProduct() = default;

    MeshProduct& operator=(MeshProduct const& x) = default;

    MeshProduct& operator=(MeshProduct&& x) = default;

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
};
