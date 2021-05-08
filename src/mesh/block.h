#pragma once

#include <memory>
#include <type_traits>
#include <vector>

#include "mesh.h"

class BlockBase
{
public:
    class AllDim
    {
        constexpr AllDim() noexcept = default;

    public:
        constexpr AllDim(const AllDim&) noexcept = default;
        constexpr AllDim(AllDim&&) noexcept = default;
        friend class BlockBase;
    };

    static constexpr AllDim ALL_DIM = AllDim();
};

template <int NDIMS, class ElementType>
class BlockND;


template <int NDIMS, class ElementType, bool CONTIGUOUS = true>
class BlockViewND : public BlockBase
{
public:
    /// ND memory view
    using RawView = ViewND<NDIMS, ElementType, CONTIGUOUS>;

    ///
    using Domain = MDomainND<NDIMS>;

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

    template <int ONDIMS, class OElementType, bool OCONTIGUOUS>
    friend class BlockViewND;

private:
    /// The raw view of the data
    RawView m_raw;

    /// The mesh on which this block is defined
    MeshND<NDIMS> m_mesh;

public:
    /** Constructs a new BlockND by copy
     * @param other the BlockND to copy
     */
    inline constexpr BlockViewND(const BlockViewND& other) noexcept = default;

    /** Constructs a new BlockND by move
     * @param other the BlockND to move
     */
    inline constexpr BlockViewND(BlockViewND&& other) noexcept = default;

    /** Copy-assigns a new value to this field
     * @param other the BlockND to copy
     * @return *this
     */
    inline constexpr BlockViewND& operator=(const BlockViewND& other) noexcept = default;

    /** Move-assigns a new value to this field
     * @param other the BlockND to move
     * @return *this
     */
    inline constexpr BlockViewND& operator=(BlockViewND&& other) noexcept = default;

    template <class... IndexType>
    inline constexpr reference operator()(IndexType... indices) const noexcept
    {
        return m_raw(indices...);
    }

    inline constexpr reference operator()(const MCoordND<NDIMS>& indices) const noexcept
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

    inline constexpr index_type extent(size_t r) const noexcept
    {
        return extents().extent(r);
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
     * @param other the BlockND to swap with this one
     */
    inline constexpr void swap(BlockViewND& other)
    {
        BlockViewND tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

    /** Provide access to the domain on which this field block is defined
     * @return the domain on which this field block is defined
     */
    inline constexpr MDomain domain(size_t dim) const noexcept
    {
        return MDomain(m_mesh[dim], m_raw.extent(dim));
    }

    /** Provide access to the domain on which this field block is defined
     * @return the domain on which this field block is defined
     */
    inline constexpr MDomainND<NDIMS> domain() const noexcept
    {
        return FullDom<std::make_index_sequence<NDIMS>>::eval(*this);
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

    /** Slice out some dimensions of
     */
    template <class... SliceSpecs>
    inline constexpr auto slice(SliceSpecs... slices) const
    {
        constexpr size_t RESULT_DIM = ((std::is_integral_v<SliceSpecs> ? 0 : 1) + ... + 0);
        auto&& new_view
                = std::experimental::subspan(raw_view(), std::forward<SliceSpecs>(slices)...);
        auto&& new_mesh = submesh(m_mesh, std::forward<SliceSpecs>(slices)...);
        return BlockViewND<
                RESULT_DIM,
                ElementType,
                new_view.is_always_contiguous()>(new_mesh, new_view);
    }

    /** Duplicate the data of this view
     * @return a copy of the data of this view
     */
    inline constexpr BlockND<NDIMS, ElementType> duplicate() const
    {
        return BlockND<NDIMS, ElementType>(*this);
    }

protected:
    inline constexpr BlockViewND(MeshND<NDIMS> mesh, RawView raw_view)
        : m_raw(raw_view)
        , m_mesh(mesh)
    {
    }

private:
    template <class>
    struct FullDom;
    template <size_t... IDXS>
    struct FullDom<std::index_sequence<IDXS...>>
    {
        static MDomainND<sizeof...(IDXS)> eval(
                const BlockViewND<sizeof...(IDXS), ElementType, CONTIGUOUS>& self)
        {
            return MDomainND<sizeof...(IDXS)> {self.domain(IDXS)...};
        }
    };
};

using DBlockView1D = BlockViewND<1, double>;

using DBlockView2D = BlockViewND<2, double>;


template <int NDIMS, class ElementType, bool CONTIGUOUS, bool OCONTIGUOUS, class... Indices>
inline BlockViewND<NDIMS, ElementType, CONTIGUOUS>& copy(
        BlockViewND<NDIMS, ElementType, CONTIGUOUS>& to,
        BlockViewND<NDIMS, ElementType, OCONTIGUOUS> const& from,
        Indices... idxs) noexcept
{
    assert(to.extents() == from.extents());
    if constexpr (sizeof...(Indices) == to.rank()) {
        to(idxs...) = from(idxs...);
    } else {
        for (size_t ii = 0; ii < to.extent(sizeof...(Indices)); ++ii) {
            copy(to, from, idxs..., ii);
        }
    }
    return to;
}

template <int NDIMS, class ElementType>
class BlockND : public BlockViewND<NDIMS, ElementType>
{
public:
    /// ND view on this block
    using BlockView = BlockViewND<NDIMS, ElementType>;

    /// ND memory view
    using RawView = typename BlockView::RawView;

    using Domain = typename BlockView::Domain;

    using extents_type = typename BlockView::extents_type;

    using layout_type = typename BlockView::layout_type;

    using accessor_type = typename BlockView::accessor_type;

    using mapping_type = typename BlockView::mapping_type;

    using element_type = typename BlockView::element_type;

    using value_type = typename BlockView::value_type;

    using index_type = typename BlockView::index_type;

    using difference_type = typename BlockView::difference_type;

    using pointer = typename BlockView::pointer;

    using reference = typename BlockView::reference;

public:
    /** Construct a BlockND on a domain with uninitialized values
     */
    inline constexpr BlockND(const MDomainND<NDIMS>& domain)
        : BlockView(
                array_apply([](const MDomain& d) { return Mesh(d.mesh()); }, domain),
                RawView(new double[domain.size()],
                        ExtentsND<NDIMS>(
                                array_apply([](const MDomain& d) { return d.size(); }, domain))))
    {
    }

    /** Constructs a new BlockND by copy
     * @param other the BlockND to copy
     */
    inline constexpr BlockND(const BlockND& other) = default;

    /** Constructs a new BlockND by move
     * @param other the BlockND to move
     */
    inline constexpr BlockND(BlockND&& other) = default;

    /** Constructs a new BlockND by copy from a view
     * @param the view to copy
     */
    template <int ONDIMS, class OElementType, bool OCONTIGUOUS>
    inline constexpr BlockND(const BlockViewND<ONDIMS, OElementType, OCONTIGUOUS>& other)
        : BlockND(other.domain())
    {
        *this = other;
    }

    inline ~BlockND()
    {
        delete[] this->raw_view().data();
    }

    /** Copy-assigns a new value to this field
     * @param other the BlockND to copy
     * @return *this
     */
    inline constexpr BlockND& operator=(const BlockND& other) = default;

    /** Move-assigns a new value to this field
     * @param other the BlockND to move
     * @return *this
     */
    inline constexpr BlockND& operator=(BlockND&& other) = default;

    /** Copy-assigns a new value to this field
     * @param other the BlockND to copy
     * @return *this
     */
    template <int ONDIMS, class OElementType, bool OCONTIGUOUS>
    inline BlockND& operator=(const BlockViewND<ONDIMS, OElementType, OCONTIGUOUS>& other)
    {
        copy(*this, other);
        return *this;
    }

    /** Swaps this field with another
     * @param other the BlockND to swap with this one
     */
    inline constexpr void swap(BlockND& other)
    {
        BlockND tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }
};

using DBlock1D = BlockND<1, double>;

using DBlock2D = BlockND<2, double>;
