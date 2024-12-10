// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "derivative_field.hpp"
#include "derivative_field_common.hpp"

/**
 * See @ref DerivFieldMemImplementation
 */
template <class ElementType, class Domain, int NDerivs, class MemSpace = Kokkos::HostSpace>
class DerivFieldMem;

template <class ElementType, class SupportType, int NDerivs, class MemSpace>
inline constexpr bool
        enable_deriv_field<DerivFieldMem<ElementType, SupportType, NDerivs, MemSpace>> = true;

template <class ElementType, class SupportType, int NDerivs, class Allocator>
inline constexpr bool enable_data_access_methods<
        DerivFieldMem<ElementType, SupportType, NDerivs, Allocator>> = true;

template <class ElementType, class SupportType, int NDerivs, class Allocator>
inline constexpr bool
        enable_mem_type<DerivFieldMem<ElementType, SupportType, NDerivs, Allocator>> = true;


/**
 * @brief A class which holds a chunk of memory describing a field and its derivatives.
 *
 * The values of the field and the derivatives may be defined on different index ranges, but
 * the underlying mesh must be the same for both.
 *
 * @anchor DerivFieldMemImplementation
 *
 * @tparam ElementType The type of the elements inside the chunks.
 * @tparam IdxRange<DDims...> The index range on which the internal fields are defined.
 *          This index range is the physical index range on which the values are defined combined with
 *          the index range of the derivatives of interest (e.g. IdxRange<Deriv<IDimX>, IDimX, IDimY>).
 * @tparam NDerivs The number of derivatives which are defined in the dimensions where derivatives
 *          appear.
 * @tparam MemSpace The memory space where the data will be saved.
 */
template <class ElementType, class... DDims, int NDerivs, class MemSpace>
class DerivFieldMem<ElementType, IdxRange<DDims...>, NDerivs, MemSpace>
    : public DerivFieldCommon<
              FieldMem<ElementType, IdxRange<DDims...>, MemSpace>,
              IdxRange<DDims...>>
{
private:
    /// @brief The class from which this object inherits.
    using base_type = DerivFieldCommon<
            FieldMem<ElementType, IdxRange<DDims...>, MemSpace>,
            IdxRange<DDims...>>;

public:
    /// @brief The type of the elements in the chunks.
    using element_type = typename base_type::element_type;

    /**
     * @brief The IdxRange on which the chunks in this object are defined.
     *
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC IdxRange types are referred to as DiscreteDomain types.
     */
    using discrete_domain_type = typename base_type::discrete_domain_type;

    /**
     * @brief The Idx which can be used to index this object.
     *
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC Idx types are referred to as DiscreteElement types.
     */
    using discrete_element_type = typename base_type::discrete_element_type;

    /// @brief A type sequence containing all derivatives present in this object.
    using deriv_tags = typename base_type::deriv_tags;

    /// @brief A type sequence containing all dimensions for which derivatives are present in this object.
    using physical_deriv_grids = typename base_type::physical_deriv_grids;

    /// @brief A type sequence containing all the physical dimensions on which the chunks are defined.
    using physical_grids = typename base_type::physical_grids;

    /// @brief The physical index range on which the field is defined.
    using physical_idx_range_type = typename base_type::physical_idx_range_type;

    /**
     * @brief The type of the field stored in the array
     *
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC FieldMem types are referred to as Chunk types and Field types are
     * referred to as ChunkSpan/ChunkView.
     */
    using chunk_type = typename base_type::chunk_type;

    template <class, class, class, class>
    friend class DerivField;

    /// @brief The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field.
    using span_type = DerivField<
            ElementType,
            IdxRange<DDims...>,
            typename chunk_type::memory_space,
            typename chunk_type::layout_type>;

    /// @brief The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field.
    using view_type = DerivField<
            ElementType const,
            IdxRange<DDims...>,
            typename chunk_type::memory_space,
            typename chunk_type::layout_type>;

private:
    /**
     * @brief A index range describing all the index ranges where derivatives are defined.
     *
     * A index range describing all the index ranges where derivatives are defined. If this index range is
     * sliced into a 1D object then it describes the index range on which the derivative in that
     * direction is defined. This object is used internally to initialise the chunks.
     */
    using physical_deriv_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<physical_deriv_grids>;

    /// @brief The IdxRange which describes the derivatives present on each chunk.
    using discrete_deriv_idx_range_type = typename base_type::discrete_deriv_idx_range_type;

    /** @brief The Idx which describes the order of the derivatives in each dimension.
     * (e.g. second-order derivative).
     */
    using discrete_deriv_index_type = typename base_type::discrete_deriv_index_type;

    using mapping_type = typename chunk_type::mapping_type;

    using extents_type = typename chunk_type::extents_type;

    using allocator_type = ddc::KokkosAllocator<element_type, typename chunk_type::memory_space>;

    using internal_mdspan_type = typename base_type::internal_mdspan_type;

    /// @brief The number of chunks which must be created to describe this object.
    using base_type::n_fields;

private:
    /// @brief A function to get the index range along direction Tag for the ArrayIndex-th element of internal_fields.
    template <std::size_t ArrayIndex, class Tag>
    IdxRange<Tag> get_idx_range(
            physical_idx_range_type val_idx_range,
            physical_deriv_idx_range_type deriv_idx_ranges)
    {
        // If the Tag describes a dimension for which a derivative is defined
        if constexpr (ddc::in_tags_v<Tag, physical_deriv_grids>) {
            // If the FieldMem at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, physical_deriv_grids>)) {
                return ddc::select<Tag>(deriv_idx_ranges);
            }
            // If the FieldMem at this index doesn't contain derivatives of this dimension
            else {
                return ddc::select<Tag>(val_idx_range);
            }
        }
        // If the Tag describes a derivative
        else if constexpr (ddc::in_tags_v<Tag, deriv_tags>) {
            // If the FieldMem at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, deriv_tags>)) {
                return IdxRange<Tag>(Idx<Tag>(1), IdxStep<Tag>(NDerivs));
            }
            // If the FieldMem at this index doesn't contain derivatives of this dimension
            else {
                return IdxRange<Tag>(Idx<Tag>(0), IdxStep<Tag>(1));
            }
        }
        // If the Tag describes a dimension for which derivatives are not defined.
        else {
            return ddc::select<Tag>(val_idx_range);
        }
    }

    /// @brief Get the index range of the FieldMem at the index ArrayIndex in internal_fields.
    template <std::size_t ArrayIndex>
    IdxRange<DDims...> get_chunk_idx_range(
            physical_idx_range_type val_idx_range,
            physical_deriv_idx_range_type deriv_idx_ranges)
    {
        return IdxRange<DDims...>(
                get_idx_range<ArrayIndex, DDims>(val_idx_range, deriv_idx_ranges)...);
    }

    /// @brief Get the size of the mdspan in dimension QueryDDim, at the index ArrayIndex.
    template <class QueryDDim, std::size_t ArrayIndex>
    std::size_t get_mdspan_size()
    {
        if constexpr (ddc::in_tags_v<QueryDDim, deriv_tags>) {
            if (ArrayIndex & (1 << ddc::type_seq_rank_v<QueryDDim, deriv_tags>)) {
                return NDerivs;
            } else {
                return 1;
            }
        } else if constexpr (ddc::in_tags_v<QueryDDim, physical_deriv_grids>) {
            if constexpr (
                    ArrayIndex & (1 << ddc::type_seq_rank_v<ddc::Deriv<QueryDDim>, deriv_tags>)) {
                IdxRangeSlice<QueryDDim> idx_range_local(base_type::m_cross_derivative_idx_range);
                return idx_range_local.extents().value();
            } else {
                IdxRange<QueryDDim> idx_range_local(base_type::m_physical_idx_range);
                return idx_range_local.extents().value();
            }
        } else {
            IdxRange<QueryDDim> idx_range_local(base_type::m_physical_idx_range);
            return idx_range_local.extents().value();
        }
    }

    /// @brief Make the internal mdspan that will be saved in internal_fields at the index ArrayIndex.
    template <std::size_t ArrayIndex>
    std::enable_if_t<std::is_constructible_v<mapping_type, extents_type>, internal_mdspan_type>
    make_internal_mdspan(allocator_type allocator)
    {
        std::size_t alloc_size(((get_mdspan_size<DDims, ArrayIndex>()) * ...));
        element_type* ptr = allocator.allocate("", alloc_size);

        extents_type extents_r(get_mdspan_size<DDims, ArrayIndex>()...);
        mapping_type mapping_r(extents_r);

        std::array<std::size_t, sizeof...(DDims)> strides_s {
                mapping_r.stride(ddc::type_seq_rank_v<DDims, ddc::detail::TypeSeq<DDims...>>)...};
        Kokkos::layout_stride::mapping<extents_type> mapping_s(extents_r, strides_s);
        return internal_mdspan_type(ptr, mapping_s);
    }

    /// @brief Initialise the chunks inside internal_fields.
    template <std::size_t... ArrayIndex>
    void initialise_chunks(allocator_type allocator, std::index_sequence<ArrayIndex...>)
    {
        ((base_type::internal_fields[ArrayIndex] = make_internal_mdspan<ArrayIndex>(allocator)),
         ...);
    }

public:
    /**
     * @brief The constructor for DerivFieldMem. The constructor initialises the chunks using
     * the provided index ranges.
     *
     * @param val_idx_range The index range on which the values of the field are defined.
     * @param m_deriv_idx_range The 1D sub-index ranges on which the derivatives of the field are defined.
     */
    template <class... DerivDoms>
    DerivFieldMem(
            physical_idx_range_type val_idx_range,
            IdxRangeSlice<DerivDoms>... m_deriv_idx_range)
        : DerivFieldMem(allocator_type(), val_idx_range, m_deriv_idx_range...)
    {
    }

    /**
     * @brief The constructor for DerivFieldMem. The constructor initialises the chunks using
     * the provided index ranges.
     *
     * @param allocator The object which allocates the memory on the CPU or GPU.
     * @param val_idx_range The index range on which the values of the field are defined.
     * @param m_deriv_idx_range The 1D sub-index ranges on which the derivatives of the field are defined.
     */
    template <class... DerivDoms>
    DerivFieldMem(
            allocator_type allocator,
            physical_idx_range_type val_idx_range,
            IdxRangeSlice<DerivDoms>... m_deriv_idx_range)
        : base_type(
                val_idx_range,
                discrete_deriv_idx_range_type(IdxRange<ddc::Deriv<DerivDoms>>(
                        Idx<ddc::Deriv<DerivDoms>>(1),
                        IdxStep<ddc::Deriv<DerivDoms>>(NDerivs))...),
                to_subidx_range_collection<physical_deriv_grids>(m_deriv_idx_range...))
    {
        static_assert(
                ddc::type_seq_same_v<ddc::detail::TypeSeq<DerivDoms...>, physical_deriv_grids>);
        initialise_chunks(allocator, std::make_integer_sequence<std::size_t, n_fields> {});
    }

    /// Defaulted destructor
    ~DerivFieldMem() = default;

    /// Deleted copy operator
    DerivFieldMem& operator=(DerivFieldMem const& other) = delete;

    /** Move-assigns a new value to this field
     * @param other the FieldMem to move
     * @return *this
     */
    DerivFieldMem& operator=(DerivFieldMem&& other) = default;

    /**
     * @brief Get a modifiable reference to an element from a constant field.
     * A Idx describes the element of interest. If information about the derivatives is
     * missing then it is assumed that the 0-th order derivative is requested.
     *
     * @param elems The element of interest.
     *
     * @returns The requested element.
     */
    template <class... DElem>
    element_type& operator()(DElem... elems) noexcept
    {
        static_assert((ddc::is_discrete_element_v<DElem> && ...));
        using full_index_type = detail::combine_t<DElem...>;
        full_index_type elem(elems...);
        return base_type::get_internal_field(elem)();
    }

    /**
     * @brief Get an element from a constant field.
     * A Idx describes the element of interest. If information about the derivatives is
     * missing then it is assumed that the 0-th order derivative is requested.
     *
     * @param elems The element of interest.
     *
     * @returns The requested element.
     */
    template <class... DElem>
    element_type const& operator()(DElem... elems) const noexcept
    {
        static_assert((ddc::is_discrete_element_v<DElem> && ...));
        using full_index_type = detail::combine_t<DElem...>;
        full_index_type elem(elems...);
        return base_type::get_internal_field(elem)();
    }

    /**
     * @brief Get a constant DerivField of this field.
     *
     * This function is designed to match the equivalent function in DDC. In Gysela it should
     * not be called directly. Instead the global function get_const_field should be used.
     *
     * @returns A constant span of this field.
     */
    view_type span_cview() const
    {
        return view_type(*this);
    }

    /**
     * @brief Get a constant DerivField of this field.
     *
     * This function is designed to match the equivalent function in DDC. In Gysela it should
     * not be called directly. Instead the global function get_field should be used.
     *
     * @returns A constant span of this field.
     */
    view_type span_view() const
    {
        return view_type(*this);
    }

    /**
     * @brief Get a modifiable DerivField of this field.
     *
     * This function is designed to match the equivalent function in DDC. In Gysela it should
     * not be called directly. Instead the global function get_field should be used.
     *
     * @returns A span of this field.
     */
    span_type span_view()
    {
        return span_type(*this);
    }
};

namespace detail {
template <class NewMemorySpace, class ElementType, class SupportType, int NDerivs, class MemSpace>
struct OnMemorySpace<NewMemorySpace, DerivFieldMem<ElementType, SupportType, NDerivs, MemSpace>>
{
    using type = DerivFieldMem<ElementType, SupportType, NDerivs, NewMemorySpace>;
};
} // namespace detail
