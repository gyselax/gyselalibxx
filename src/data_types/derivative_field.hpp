// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "derivative_field_common.hpp"

template <class, class, int, class>
class DerivFieldMem;

/**
 * See @ref DerivFieldImplementation
 */
template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right,
        class MemorySpace = Kokkos::HostSpace>
class DerivField;

template <class ElementType, class SupportType, class LayoutStridedPolicy, class MemorySpace>
inline constexpr bool enable_deriv_field<
        DerivField<ElementType, SupportType, LayoutStridedPolicy, MemorySpace>> = true;

template <class ElementType, class SupportType, class LayoutStridedPolicy, class MemorySpace>
inline constexpr bool enable_borrowed_deriv_field<
        DerivField<ElementType, SupportType, LayoutStridedPolicy, MemorySpace>> = true;

template <class ElementType, class SupportType, class LayoutStridedPolicy, class MemorySpace>
inline constexpr bool enable_data_access_methods<
        DerivField<ElementType, SupportType, LayoutStridedPolicy, MemorySpace>> = true;

namespace ddcHelper {

/**
 * @brief Copy the contents of one DerivField into another.
 *
 * @param dst The DerivField where the data will be saved.
 * @param src The DerivField whose data will be copied.
 */
template <
        class FieldDst,
        class FieldSrc,
        std::enable_if_t<
                is_borrowed_deriv_field_v<FieldDst> && is_borrowed_deriv_field_v<FieldSrc>,
                bool> = true>
auto deepcopy(FieldDst&& dst, FieldSrc&& src)
{
    assert(dst.get_values_field().domain().extents() == src.get_values_field().domain().extents());

    DerivField dst_field = dst.span_view();
    DerivField src_field = src.span_view();

    dst_field.deepcopy(src_field);

    return dst_field;
}

/**
 * @brief Copy the contents of one DerivField into another.
 *
 * @param execution_space The Kokkos execution space where the copy will be carried out.
 * @param dst The DerivField where the data will be saved.
 * @param src The DerivField whose data will be copied.
 */
template <class ExecSpace, class FieldDst, class FieldSrc>
auto deepcopy(ExecSpace const& execution_space, FieldDst&& dst, FieldSrc&& src)
{
    static_assert(is_borrowed_deriv_field_v<FieldDst>);
    static_assert(is_borrowed_deriv_field_v<FieldSrc>);

    assert(dst.get_values_field().idx_range().extents()
           == src.get_values_field().idx_range().extents());

    DerivField dst_field = get_field(dst);
    DerivField src_field = get_field(src);

    dst_field.deepcopy(execution_space, src_field);

    return dst_field;
}

} // namespace ddcHelper

/**
 * @brief A class which holds references to chunks of memory describing a field and its derivatives.
 *
 * The values of the field and the derivatives may be defined on different index ranges, but
 * the underlying mesh must be the same for both.
 *
 * @anchor DerivFieldImplementation
 *
 * @tparam ElementType The type of the elements inside the chunks.
 * @tparam IdxRange<DDims...> The index range on which the internal fields are defined.
 *          This index range is the physical index range on which the values are defined combined with
 *          the index range of the derivatives of interest (e.g. IdxRange<Deriv<IDimX>, IDimX, IDimY>).
 * @tparam LayoutStridedPolicy The way in which the memory is laid out in memory (contiguous in
 *          the leading/trailing dimension, strided, etc).
 * @tparam MemorySpace The memory space where the data is saved (CPU/GPU).
 */
template <class ElementType, class... DDims, class LayoutStridedPolicy, class MemorySpace>
class DerivField<ElementType, IdxRange<DDims...>, LayoutStridedPolicy, MemorySpace>
    : public DerivFieldCommon<
              Field<ElementType, IdxRange<DDims...>, LayoutStridedPolicy, MemorySpace>,
              IdxRange<DDims...>>
{
private:
    /// @brief The class from which this object inherits.
    using base_type = DerivFieldCommon<
            Field<ElementType, IdxRange<DDims...>, LayoutStridedPolicy, MemorySpace>,
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
    /// @brief The IdxRange on which the fields in this object are defined.
    using index_range_type = typename base_type::index_range_type;

    /**
     * @brief The Idx which can be used to index this object.
     *
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC Idx types are referred to as DiscreteElement types.
     */
    using discrete_element_type = typename base_type::discrete_element_type;

    /// @brief A type sequence containing all derivatives present in this object.
    using deriv_tags = typename base_type::deriv_tags;

    /// @brief A type sequence containing all grid types for which derivatives are present in this object.
    using physical_deriv_grids = typename detail::strip_deriv_t<deriv_tags>;

    /// @brief A type sequence containing all the grids on which the fields are defined.
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

    /// @brief The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field.
    using span_type = DerivField<ElementType, IdxRange<DDims...>, LayoutStridedPolicy, MemorySpace>;

    /// @brief The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field.
    using view_type
            = DerivField<ElementType const, IdxRange<DDims...>, LayoutStridedPolicy, MemorySpace>;

private:
    /// @brief The IdxRange which describes the derivatives present on each chunk.
    using discrete_deriv_idx_range_type = typename base_type::discrete_deriv_idx_range_type;

    /** @brief The Idx which describes the order of the derivatives in each dimension.
     * (e.g. second-order derivative).
     */
    using discrete_deriv_index_type = typename base_type::discrete_deriv_index_type;

    /// @brief The type of a reference to an element of the mdspan
    using reference = typename chunk_type::reference;

    /// @brief The number of chunks which must be created to describe this object.
    using base_type::n_fields;

private:
    /** @brief Get the subindex range to be extracted from a DerivFieldMem to build the internal_chunk at position ArrayIndex
     * along dimension Tag.
     */
    template <std::size_t ArrayIndex, class Tag>
    KOKKOS_FUNCTION constexpr Idx<Tag> get_chunk_subidx_range_1d_idx()
    {
        // The Tag describes a dimension for which a derivative is defined
        if constexpr (ddc::in_tags_v<Tag, deriv_tags>) {
            // If the Field at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, deriv_tags>)) {
                return Idx<Tag>(1);
            }
            // If the Field at this index doesn't contain derivatives of this dimension
            else {
                return Idx<Tag>(0);
            }
        } else {
            // Empty IdxRange to be discarded
            return Idx<Tag>();
        }
    }

    /// @brief Get the subindex range to be extracted from a DerivFieldMem to build the internal_chunk at position ArrayIndex
    template <std::size_t ArrayIndex>
    KOKKOS_FUNCTION constexpr discrete_deriv_index_type get_chunk_subidx_range_idx()
    {
        return discrete_deriv_index_type(get_chunk_subidx_range_1d_idx<ArrayIndex, DDims>()...);
    }

    /// @brief Initialise fields inside internal_fields.
    template <class Field, std::size_t... ArrayIndex>
    KOKKOS_FUNCTION constexpr void initialise_fields(
            Field const& chunks,
            std::index_sequence<ArrayIndex...>)
    {
        ((base_type::internal_fields[ArrayIndex] = chunks.internal_fields[chunks.get_array_index(
                  get_chunk_subidx_range_idx<ArrayIndex>())]),
         ...);
    }

    auto get_kokkos_view_from_internal_chunk(int index)
    {
        typename base_type::internal_mdspan_type field = base_type::internal_fields[index];
        auto kokkos_layout = ddc::detail::build_kokkos_layout(
                field.extents(),
                field.mapping(),
                std::make_index_sequence<sizeof...(DDims)> {});
        return Kokkos::View<
                ddc::detail::mdspan_to_kokkos_element_t<ElementType, sizeof...(DDims)>,
                decltype(kokkos_layout),
                MemorySpace>(field.data_handle(), kokkos_layout);
    }

public:
    /**
     * @brief Copy-construct a DerivField
     *
     * @param other The DerivField being copied.
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivField(DerivField const& other) = default;

    /**
     * @brief Move-construct a DerivField
     *
     * @param other The DerivField being moved.
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivField(DerivField&& other) = default;

    /** 
     * @brief Constructs a new DerivField containing a modifiable view on the data in a DerivFieldMem.
     *
     * @param field The DerivFieldMem to view.
     */
    template <
            class OElementType,
            int NDerivs,
            class Allocator,
            class = std::enable_if_t<std::is_same_v<typename Allocator::memory_space, MemorySpace>>>
    explicit constexpr DerivField(
            DerivFieldMem<OElementType, index_range_type, NDerivs, Allocator>& field)
        : base_type(
                field.m_physical_idx_range,
                field.m_deriv_idx_range,
                field.m_cross_derivative_idx_range)
    {
        initialise_fields(field, std::make_integer_sequence<std::size_t, n_fields> {});
    }

    /** 
     * @brief Constructs a new DerivField containing a constant view on the data in a DerivFieldMem.
     *
     * @param field The DerivFieldMem to view.
     */
    // Disabled by SFINAE if `ElementType` is not `const` to avoid write access
    template <
            class OElementType,
            class SFINAEElementType = ElementType,
            class = std::enable_if_t<std::is_const_v<SFINAEElementType>>,
            int NDerivs,
            class Allocator,
            class = std::enable_if_t<std::is_same_v<typename Allocator::memory_space, MemorySpace>>>
    explicit constexpr DerivField(
            DerivFieldMem<OElementType, index_range_type, NDerivs, Allocator> const& field)
        : base_type(
                field.m_physical_idx_range,
                field.m_deriv_idx_range,
                field.m_cross_derivative_idx_range)
    {
        initialise_fields(field, std::make_integer_sequence<std::size_t, n_fields> {});
    }

    /**
     * @brief Copy construct a DerivField. The element type may be changed to a complatible type.
     * (e.g. double -> const double).
     *
     * @param field The DerivField to be copied.
     */
    template <class OElementType>
    KOKKOS_FUNCTION constexpr DerivField(
            DerivField<OElementType, index_range_type, LayoutStridedPolicy, MemorySpace> const&
                    field)
        : base_type(
                field.m_physical_idx_range,
                field.m_deriv_idx_range,
                field.m_cross_derivative_idx_range)
    {
        initialise_fields(field, std::make_integer_sequence<std::size_t, n_fields> {});
    }

    KOKKOS_DEFAULTED_FUNCTION ~DerivField() = default;

    /** Copy-assigns a new value to this DerivField, yields a new view to the same data
     * @param other the DerivField to copy
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivField& operator=(DerivField const& other) = default;

    /** Move-assigns a new value to this DerivField
     * @param other the DerivField to move
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivField& operator=(DerivField&& other) = default;

    /**
     * @brief Copy the source DerivField into this DerivField using Kokkos::deep_copy. 
     *
     * @param src The DerivField containing the data to be copied.
     */
    template <class OElementType, class OLayoutStridedPolicy, class OMemorySpace>
    void deepcopy(
            DerivField<OElementType, index_range_type, OLayoutStridedPolicy, OMemorySpace> src)
    {
        for (int i(0); i < n_fields; ++i) {
            auto kokkos_span = get_kokkos_view_from_internal_chunk(i);
            auto src_kokkos_span = src.get_kokkos_view_from_internal_chunk(i);
            Kokkos::deep_copy(kokkos_span, src_kokkos_span);
        }
    }

    /**
     * @brief Copy the source DerivField into this DerivField using Kokkos::deep_copy.
     *
     * @param execution_space The execution space on which the copy will be carried out.
     * @param src The DerivField containing the data to be copied.
     */
    template <class ExecSpace, class OElementType, class OLayoutStridedPolicy, class OMemorySpace>
    void deepcopy(
            ExecSpace const& execution_space,
            DerivField<OElementType, index_range_type, OLayoutStridedPolicy, OMemorySpace> src)
    {
        for (int i(0); i < n_fields; ++i) {
            auto kokkos_span = get_kokkos_view_from_internal_chunk(i);
            auto src_kokkos_span = src.get_kokkos_view_from_internal_chunk(i);
            Kokkos::deep_copy(execution_space, kokkos_span, src_kokkos_span);
        }
    }

    /**
     * @brief Get an element from a constant field.
     * An Idx describes the element of interest. If information about the derivatives is
     * missing then it is assumed that the 0-th order derivative is requested.
     *
     * @param elems The element of interest.
     *
     * @returns The requested element.
     */
    template <class... DElem>
    KOKKOS_FUNCTION constexpr reference operator()(DElem... elems) const noexcept
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
    KOKKOS_FUNCTION constexpr view_type span_cview() const
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
    KOKKOS_FUNCTION constexpr span_type span_view() const
    {
        return *this;
    }
};

template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right,
        class MemorySpace = Kokkos::HostSpace>
using DerivConstField
        = DerivField<ElementType const, SupportType, LayoutStridedPolicy, MemorySpace>;

namespace detail {
template <
        class NewMemorySpace,
        class ElementType,
        class SupportType,
        class Layout,
        class MemorySpace>
struct OnMemorySpace<NewMemorySpace, DerivField<ElementType, SupportType, Layout, MemorySpace>>
{
    using type = DerivField<ElementType, SupportType, Layout, NewMemorySpace>;
};
}; // namespace detail
