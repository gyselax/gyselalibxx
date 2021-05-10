#pragma once

#include <cassert>
#include <memory>
#include <type_traits>
#include <vector>

#include "mdomain.h"
#include "view.h"

template <class, class>
class Block;


template <class, class, bool = true>
class BlockView;

template <class... Tags, class ElementType, bool CONTIGUOUS>
class BlockView<MDomain<Tags...>, ElementType, CONTIGUOUS>
{
public:
    /// ND memory view
    using RawView = ViewND<sizeof...(Tags), ElementType, CONTIGUOUS>;

    using MDomain_ = MDomain<Tags...>;

    using Mesh = typename MDomain_::Mesh;

    using MCoord_ = typename MDomain_::MCoord_;

    using extents_type = typename RawView::extents_type;

    using layout_type = typename RawView::layout_type;

    using accessor_type = typename RawView::accessor_type;

    using mapping_type = typename RawView::mapping_type;

    using element_type = typename RawView::element_type;

    using value_type = typename RawView::value_type;

    using index_type = typename RawView::index_type;

    using difference_type = typename RawView::difference_type;

    using pointer = typename RawView::pointer;

    using reference = typename RawView::reference;

    template <class, class, bool>
    friend class BlockView;

private:
    /// The raw view of the data
    RawView m_raw;

    /// The mesh on which this block is defined
    Mesh m_mesh;

public:
    /** Constructs a new Block by copy
     * @param other the Block to copy
     */
    inline constexpr BlockView(const BlockView& other) noexcept = default;

    /** Constructs a new Block by move
     * @param other the Block to move
     */
    inline constexpr BlockView(BlockView&& other) noexcept = default;

    /** Copy-assigns a new value to this field
     * @param other the Block to copy
     * @return *this
     */
    inline constexpr BlockView& operator=(const BlockView& other) noexcept = default;

    /** Move-assigns a new value to this field
     * @param other the Block to move
     * @return *this
     */
    inline constexpr BlockView& operator=(BlockView&& other) noexcept = default;

    template <class... IndexType>
    inline constexpr reference operator()(IndexType... indices) const noexcept
    {
        return m_raw(indices...);
    }

    inline constexpr reference operator()(const MCoord_& indices) const noexcept
    {
        return m_raw(indices);
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
        return m_raw.stride();
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
    inline constexpr MDomain_ domain() const noexcept
    {
        return MDomain_(mesh(), mcoord_end<Tags...>(raw_view().extents()));
    }

    /** Provide access to the domain on which this block is defined
     * @return the domain on which this block is defined
     */
    template <class... OTags>
    inline constexpr MDomain<OTags...> domain() const noexcept
    {
        return MDomain<OTags...>(domain());
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    inline constexpr RawView raw_view()
    {
        return m_raw;
    }

    /** Provide a constant view of the data
     * @return a constant view of the data
     */
    inline constexpr const RawView raw_view() const
    {
        return m_raw;
    }

    template <class>
    struct Slicer;

    template <class... STags>
    struct Slicer<RegularMesh<STags...>>
    {
        template <class... SliceSpecs>
        static inline constexpr auto mesh_for(
                const RegularMesh<STags...>& mesh,
                std::experimental::all_type,
                SliceSpecs... slices)
        {
        }

        template <class... SliceSpecs>
        static inline constexpr auto mesh_for(
                const RegularMesh<STags...>& mesh,
                SliceSpecs... slices)
        {
        }
    };

    /** Slice out some dimensions
     * @param slices the coordinates to 
     */
    template <class... SliceSpecs>
    inline constexpr auto slice(SliceSpecs&&... slices) const
    {
        auto view = subspan(raw_view(), std::forward<SliceSpecs>(slices)...);
        auto mesh = submesh(m_mesh, std::forward<SliceSpecs>(slices)...);
        return (BlockView(mesh, view));
    }

    /** Duplicate the data of this view
     * @return a copy of the data of this view
     */
    inline constexpr Block<MDomain<Tags...>, ElementType> duplicate() const
    {
        Block<MDomain<Tags...>, ElementType> result(this->domain());
        deepcopy(result, *this);
        return result;
    }

protected:
    inline constexpr BlockView(const Mesh& mesh, RawView raw_view) : m_raw(raw_view), m_mesh(mesh)
    {
    }
};

using DBlockViewX = BlockView<MDomain<Dim::X>, double>;

using DBlockViewVx = BlockView<MDomain<Dim::Vx>, double>;

using DBlockViewXVx = BlockView<MDomain<Dim::X, Dim::Vx>, double>;

template <class... Tags, class View>
RegularMDomain<Tags...> get_domain(const View& v)
{
    return v.template domain<Tags...>();
}


template <class... Tags, class ElementType, bool CONTIGUOUS, bool OCONTIGUOUS, class... Indices>
inline BlockView<MDomain<Tags...>, ElementType, CONTIGUOUS>& deepcopy(
        BlockView<MDomain<Tags...>, ElementType, CONTIGUOUS>& to,
        BlockView<MDomain<Tags...>, ElementType, OCONTIGUOUS> const& from,
        Indices... idxs) noexcept
{
    assert(to.extents() == from.extents());
    //TODO !
    //     if constexpr (sizeof...(Indices) == to.rank()) {
    //         to(idxs...) = from(idxs...);
    //     } else {
    //         for (size_t ii = 0; ii < to.extent(sizeof...(Indices)); ++ii) {
    //             copy(to, from, idxs..., ii);
    //         }
    //     }
    return to;
}


template <class... Tags, class ElementType>
class Block<MDomain<Tags...>, ElementType> : public BlockView<MDomain<Tags...>, ElementType>
{
public:
    /// ND view on this block
    using BlockView_ = BlockView<MDomain<Tags...>, ElementType>;

    /// ND memory view
    using RawView = ViewND<sizeof...(Tags), ElementType>;

    using MDomain_ = MDomain<Tags...>;

    using Mesh = typename MDomain_::Mesh;

    using MCoord_ = typename MDomain_::MCoord_;

    using extents_type = typename BlockView_::extents_type;

    using layout_type = typename BlockView_::layout_type;

    using accessor_type = typename BlockView_::accessor_type;

    using mapping_type = typename BlockView_::mapping_type;

    using element_type = typename BlockView_::element_type;

    using value_type = typename BlockView_::value_type;

    using index_type = typename BlockView_::index_type;

    using difference_type = typename BlockView_::difference_type;

    using pointer = typename BlockView_::pointer;

    using reference = typename BlockView_::reference;

public:
    /** Construct a Block on a domain with uninitialized values
     */
    template <class... OTags>
    explicit inline constexpr Block(const MDomain<OTags...>& domain)
        : BlockView_(
                domain,
                RawView(new double[domain.size()],
                        ExtentsND<sizeof...(Tags)>(domain.template extent<Tags>()...)))
    {
    }

    /** Constructs a new Block by copy
     * 
     * This is deleted, one should use deepcopy
     * @param other the Block to copy
     */
    inline constexpr Block(const Block& other) = delete;

    /** Constructs a new Block by move
     * @param other the Block to move
     */
    inline constexpr Block(Block&& other) = default;

    inline ~Block()
    {
        delete[] this->raw_view().data();
    }

    /** Copy-assigns a new value to this field
     * @param other the Block to copy
     * @return *this
     */
    inline constexpr Block& operator=(const Block& other) = default;

    /** Move-assigns a new value to this field
     * @param other the Block to move
     * @return *this
     */
    inline constexpr Block& operator=(Block&& other) = default;

    /** Copy-assigns a new value to this field
     * @param other the Block to copy
     * @return *this
     */
    template <class... OTags, class OElementType>
    inline Block& operator=(Block<MDomain<OTags...>, OElementType>&& other)
    {
        copy(*this, other);
        return *this;
    }

    /** Swaps this field with another
     * @param other the Block to swap with this one
     */
    inline constexpr void swap(Block& other)
    {
        Block tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }
};

template <class ElementType>
using BlockX = Block<MDomain<Dim::X>, ElementType>;

using DBlockX = BlockX<double>;

template <class ElementType>
using BlockVx = Block<MDomain<Dim::Vx>, ElementType>;

using DBlockVx = BlockVx<double>;

template <class ElementType>
using BlockXVx = Block<MDomain<Dim::X, Dim::Vx>, ElementType>;

using DBlockXVx = BlockXVx<double>;
