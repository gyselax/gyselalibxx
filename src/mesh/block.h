#pragma once

#include <memory>
#include <span>
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
class BlockViewND
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

private:
    /// The raw view of the data
    RawView m_raw;

    /// The mesh on which this block is defined
    MeshND<NDIMS> m_mesh;

public:
    /** Constructs a new BlockND by copy
     * @param other the BlockND to copy
     */
    constexpr BlockViewND(const BlockViewND& other) noexcept = default;

    /** Constructs a new BlockND by move
     * @param other the BlockND to move
     */
    constexpr BlockViewND(BlockViewND&& other) noexcept = default;

    /** Copy-assigns a new value to this field
     * @param other the BlockND to copy
     * @return *this
     */
    constexpr BlockViewND& operator=(const BlockViewND& other) noexcept = default;

    /** Move-assigns a new value to this field
     * @param other the BlockND to move
     * @return *this
     */
    constexpr BlockViewND& operator=(BlockViewND&& other) noexcept = default;

    template <class... IndexType>
    constexpr reference operator()(IndexType... indices) const noexcept
    {
        return m_raw(indices...);
    }

    constexpr reference operator()(const MCoordND<NDIMS>& indices) const noexcept
    {
        return m_raw(indices);
    }

    accessor_type accessor() const
    {
        return m_raw.accessor();
    }

    static constexpr int rank() noexcept
    {
        return extents_type::rank();
    }

    static constexpr int rank_dynamic() noexcept
    {
        return extents_type::rank_dynamic();
    }

    static constexpr index_type static_extent(size_t r) noexcept
    {
        return extents_type::static_extent(r);
    }

    constexpr extents_type extents() const noexcept
    {
        return m_raw.extents();
    }

    constexpr index_type extent(size_t r) const noexcept
    {
        return extents().extent(r);
    }

    constexpr index_type size() const noexcept
    {
        return m_raw.size();
    }

    constexpr index_type unique_size() const noexcept
    {
        return m_raw.unique_size();
    }

    static constexpr bool is_always_unique() noexcept
    {
        return mapping_type::is_always_unique();
    }

    static constexpr bool is_always_contiguous() noexcept
    {
        return mapping_type::is_always_contiguous();
    }

    static constexpr bool is_always_strided() noexcept
    {
        return mapping_type::is_always_strided();
    }

    constexpr mapping_type mapping() const noexcept
    {
        return m_raw.mapping();
    }

    constexpr bool is_unique() const noexcept
    {
        return m_raw.is_unique();
    }

    constexpr bool is_contiguous() const noexcept
    {
        return m_raw.is_contiguous();
    }

    constexpr bool is_strided() const noexcept
    {
        return m_raw.is_strided();
    }

    constexpr index_type stride(size_t r) const
    {
        return m_raw.stride();
    }

    /** Swaps this field with another
     * @param other the BlockND to swap with this one
     */
    void swap(BlockViewND& other)
    {
        BlockViewND tmp = std::move(other);
        other = std::move(*this);
        *this = std::move(tmp);
    }

    /** Provide access to the domain on which this field block is defined
     * @return the domain on which this field block is defined
     */
    constexpr MDomain domain(size_t dim) const noexcept
    {
        return MDomain(m_mesh[dim], m_raw.extent(dim));
    }

    /** Provide access to the domain on which this field block is defined
     * @return the domain on which this field block is defined
     */
    constexpr MDomainND<NDIMS> domain() const noexcept
    {
        return FullDom<std::make_index_sequence<NDIMS>>::eval(*this);
    }

    /** Provide a modifiable view of the data
     * @return a modifiable view of the data
     */
    RawView raw_view()
    {
        return m_raw;
    }

    /** Provide a constant view of the data
     * @return a constant view of the data
     */
    const RawView raw_view() const
    {
        return m_raw;
    }

    /** Slice out some dimensions of
     */
    template <size_t... dims>
    const BlockViewND<sizeof...(dims), ElementType, false> slice() const
    {
        return m_raw;
    }

    /** Duplicate the data of this view
     * @return a copy of the data of this view
     */
    BlockND<NDIMS, ElementType> duplicate() const;

protected:
    /** Construct a BlockND on a domain with default initialized values
     */
    BlockViewND(MeshND<NDIMS> mesh, RawView raw_view) : m_raw(raw_view), m_mesh(mesh) { }

    static const BlockViewND make_const(const Domain& domain, const RawView raw_view);

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

    template <class>
    struct Slicer;
    template <size_t... IDXS>
    struct Slicer<std::index_sequence<IDXS...>>
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
    BlockND(const MDomainND<NDIMS>& domain)
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
    BlockND(const BlockND& other) = default;

    /** Constructs a new BlockND by copy
     * @param other the BlockND to copy
     */
    BlockND(const BlockView& other);

    /** Constructs a new BlockND by move
     * @param other the BlockND to move
     */
    BlockND(BlockND&& other) = default;

    ~BlockND()
    {
        delete[] this->raw_view().data();
    }

    /** Copy-assigns a new value to this field
     * @param other the BlockND to copy
     * @return *this
     */
    BlockND& operator=(const BlockND& other) = default;

    /** Copy-assigns a new value to this field
     * @param other the BlockND to copy
     * @return *this
     */
    BlockND& operator=(const BlockView& other);

    /** Move-assigns a new value to this field
     * @param other the BlockND to move
     * @return *this
     */
    BlockND& operator=(BlockND&& other) = default;

    /** Swaps this field with another
     * @param other the BlockND to swap with this one
     */
    void swap(BlockND& other);
};

using DBlock1D = BlockND<1, double>;

using DBlock2D = BlockND<2, double>;
