#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include "mdomain.h"
#include "taggedarray.h"
#include "view.h"

template <class, class>
class Block;

template <class, class, bool = true>
class BlockView;

template <class OElementType, bool O_CONTIGUOUS, class Mesh>
static BlockView<MDomainImpl<Mesh>, OElementType, O_CONTIGUOUS> make_view(
        Mesh const& mesh,
        SpanND<Mesh::rank(), OElementType, O_CONTIGUOUS> const& raw_view);

template <class Mesh, class ElementType, bool CONTIGUOUS>
class BlockView<MDomainImpl<Mesh>, ElementType, CONTIGUOUS>
{
public:
    /// ND memory view
    using raw_view_type = SpanND<Mesh::rank(), ElementType, CONTIGUOUS>;

    using mdomain_type = MDomainImpl<Mesh>;

    using mesh_type = Mesh;

    using mcoord_type = typename mdomain_type::mcoord_type;

    using extents_type = typename raw_view_type::extents_type;

    using layout_type = typename raw_view_type::layout_type;

    using accessor_type = typename raw_view_type::accessor_type;

    using mapping_type = typename raw_view_type::mapping_type;

    using element_type = typename raw_view_type::element_type;

    using value_type = typename raw_view_type::value_type;

    using index_type = typename raw_view_type::index_type;

    using difference_type = typename raw_view_type::difference_type;

    using pointer = typename raw_view_type::pointer;

    using reference = typename raw_view_type::reference;

    template <class, class, bool>
    friend class BlockView;

protected:
    template <class QTag, class... CTags>
    static auto get_slicer_for(const MCoord<CTags...>& c)
    {
        if constexpr (has_tag_v<QTag, MCoord<CTags...>>) {
            return c.template get<QTag>();
        } else {
            return std::experimental::all;
        }
    }

    template <class... SliceSpecs>
    struct Slicer
    {
        template <class... OSliceSpecs>
        static inline constexpr auto slice(const BlockView& block, OSliceSpecs&&... slices)
        {
            auto view = subspan(block.raw_view(), std::forward<OSliceSpecs>(slices)...);
            auto mesh = submesh(block.mesh(), std::forward<OSliceSpecs>(slices)...);
            return make_view<element_type, ::is_contiguous_v<decltype(view)>>(mesh, view);
        }
    };

    template <class... STags>
    struct Slicer<MCoord<STags...>>
    {
        static inline constexpr auto slice(const BlockView& block, const MCoord<STags...>& slices)
        {
            return Slicer<MCoord<STags...>, mdomain_type>::slice(block, std::move(slices));
        }
    };

    template <class... STags>
    struct Slicer<MCoord<STags...>, MDomainImpl<typename Mesh::template Mesh<>>>
    {
        template <class... SliceSpecs>
        static inline constexpr auto slice(
                const BlockView& block,
                const MCoord<STags...>&,
                SliceSpecs&&... oslices)
        {
            return Slicer<SliceSpecs...>::slice(block, oslices...);
        }
    };

    template <class... STags, class OTag0, class... OTags>
    struct Slicer<MCoord<STags...>, MDomainImpl<typename Mesh::template Mesh<OTag0, OTags...>>>
    {
        template <class... SliceSpecs>
        static inline constexpr auto slice(
                const BlockView& block,
                const MCoord<STags...>& slices,
                SliceSpecs&&... oslices)
        {
            return Slicer<MCoord<STags...>, MDomainImpl<typename Mesh::template Mesh<OTags...>>>::
                    slice(block, slices, oslices..., get_slicer_for<OTag0>(slices));
        }
    };

    /// The raw view of the data
    raw_view_type m_raw;

    /// The mesh on which this block is defined
    Mesh m_mesh;

public:
    /** Constructs a new BlockView by copy, yields a new view to the same data
     * @param other the BlockView to copy
     */
    inline constexpr BlockView(const BlockView& other) = default;

    /** Constructs a new BlockView by move
     * @param other the BlockView to move
     */
    inline constexpr BlockView(BlockView&& other) = default;

    /** Constructs a new BlockView by copy of a block, yields a new view to the same data
     * @param other the BlockView to move
     */
    template <class OElementType>
    inline constexpr BlockView(const Block<mdomain_type, OElementType>& other) noexcept
        : m_raw(other.raw_view())
        , m_mesh(other.mesh())
    {
    }

    /** Constructs a new BlockView by copy of a block, yields a new view to the same data
     * @param other the BlockView to move
     */
    template <class OElementType>
    inline constexpr BlockView(
            const BlockView<mdomain_type, OElementType, CONTIGUOUS>& other) noexcept
        : m_raw(other.raw_view())
        , m_mesh(other.mesh())
    {
    }

    /** Constructs a new BlockView from scratch
     * @param mesh the mesh that sustains the view
     * @param raw_view the raw view to the data
     */
    inline constexpr BlockView(const Mesh& mesh, raw_view_type raw_view)
        : m_raw(raw_view)
        , m_mesh(mesh)
    {
    }

    /** Copy-assigns a new value to this BlockView, yields a new view to the same data
     * @param other the BlockView to copy
     * @return *this
     */
    inline constexpr BlockView& operator=(const BlockView& other) = default;

    /** Move-assigns a new value to this BlockView
     * @param other the BlockView to move
     * @return *this
     */
    inline constexpr BlockView& operator=(BlockView&& other) = default;

    template <class QueryTag>
    static constexpr std::size_t tag_rank()
    {
        return Mesh::template tag_rank<QueryTag>();
    }

    /** Slice out some dimensions
     * @param slices the coordinates to 
     */
    template <class SliceSpec>
    inline constexpr auto operator[](SliceSpec&& slice) const
    {
        return Slicer<std::remove_cv_t<std::remove_reference_t<SliceSpec>>>::
                slice(*this, std::forward<SliceSpec>(slice));
    }

    template <class... IndexType>
    inline constexpr reference operator()(IndexType&&... indices) const noexcept
    {
        return m_raw(std::forward<IndexType>(indices)...);
    }

    template <class... OTags>
    inline constexpr reference operator()(const MCoord<OTags...>& indices) const noexcept
    {
        return m_raw(mcoord_type(indices).array());
    }

    inline constexpr reference operator()(const mcoord_type& indices) const noexcept
    {
        return m_raw(indices.array());
    }

    inline accessor_type accessor() const
    {
        return m_raw.accessor();
    }

    static inline constexpr int rank() noexcept
    {
        return extents_type::rank();
    }

    static inline constexpr int rank_dynamic() noexcept
    {
        return extents_type::rank_dynamic();
    }

    static inline constexpr index_type static_extent(size_t r) noexcept
    {
        return extents_type::static_extent(r);
    }

    inline constexpr extents_type extents() const noexcept
    {
        return m_raw.extents();
    }

    inline constexpr index_type extent(size_t dim) const noexcept
    {
        return m_raw.extent(dim);
    }

    inline constexpr index_type size() const noexcept
    {
        return m_raw.size();
    }

    inline constexpr index_type unique_size() const noexcept
    {
        return m_raw.unique_size();
    }

    static inline constexpr bool is_always_unique() noexcept
    {
        return mapping_type::is_always_unique();
    }

    static inline constexpr bool is_always_contiguous() noexcept
    {
        return mapping_type::is_always_contiguous();
    }

    static inline constexpr bool is_always_strided() noexcept
    {
        return mapping_type::is_always_strided();
    }

    inline constexpr mapping_type mapping() const noexcept
    {
        return m_raw.mapping();
    }

    inline constexpr bool is_unique() const noexcept
    {
        return m_raw.is_unique();
    }

    inline constexpr bool is_contiguous() const noexcept
    {
        return m_raw.is_contiguous();
    }

    inline constexpr bool is_strided() const noexcept
    {
        return m_raw.is_strided();
    }

    inline constexpr index_type stride(size_t r) const
    {
        return m_raw.stride(r);
    }

    /** Swaps this field with another
     * @param other the Block to swap with this one
     */
    inline constexpr void swap(BlockView& other)
    {
        BlockView tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

    /** Provide access to the mesh on which this block is defined
     * @return the mesh on which this block is defined
     */
    inline constexpr Mesh mesh() const noexcept
    {
        return m_mesh;
    }

    /** Provide access to the domain on which this block is defined
     * @return the domain on which this block is defined
     */
    inline constexpr mdomain_type domain() const noexcept
    {
        return mdomain_type(mesh(), ExtentToMCoordEnd<mcoord_type>::mcoord(raw_view().extents()));
    }

    /** Provide access to the domain on which this block is defined
     * @return the domain on which this block is defined
     */
    template <class... OTags>
    inline constexpr MDomainImpl<typename Mesh::template Mesh<OTags...>> domain() const noexcept
    {
        return MDomainImpl<typename Mesh::template Mesh<OTags...>>(domain());
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr raw_view_type raw_view()
    {
        return m_raw;
    }

    /** Provide a constant view of the data
     * @return a constant view of the data
     */
    inline constexpr const raw_view_type raw_view() const
    {
        return m_raw;
    }

    /** Slice out some dimensions
     * @param slices the coordinates to 
     */
    template <class... SliceSpecs>
    inline constexpr auto subblockview(SliceSpecs&&... slices) const
    {
        return Slicer<std::remove_cv_t<std::remove_reference_t<SliceSpecs>>...>::
                slice(*this, std::forward<SliceSpecs>(slices)...);
    }
};

/** Construct a new BlockView from scratch
 * @param[in] mesh      the mesh that sustains the view
 * @param[in] raw_view  the raw view to the data
 * @return the newly constructed view
 */
template <class OElementType, bool O_CONTIGUOUS, class Mesh>
static BlockView<MDomainImpl<Mesh>, OElementType, O_CONTIGUOUS> make_view(
        Mesh const& mesh,
        SpanND<Mesh::rank(), OElementType, O_CONTIGUOUS> const& raw_view)
{
    return BlockView<MDomainImpl<Mesh>, OElementType, O_CONTIGUOUS>(mesh, raw_view);
}

/** Access the domain (or subdomain) of a view
 * @param[out] view  the view whose domain to iterate
 * @param[in]  f     a functor taking the list of indices as parameter
 * @return the domain of view in the queried dimensions
 */
template <class... QueryTags, class Mesh, class ElementType, bool CONTIGUOUS>
MDomainImpl<typename Mesh::template Mesh<QueryTags...>> get_domain(
        const BlockView<MDomainImpl<Mesh>, ElementType, CONTIGUOUS>& v)
{
    return v.template domain<QueryTags...>();
}
