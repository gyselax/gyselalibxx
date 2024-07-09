// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "derivative_field_common.hpp"

template <class, class, int, class>
class DerivField;

template <
        class ElementType,
        class SupportType,
        class LayoutStridedPolicy = std::experimental::layout_right,
        class MemorySpace = Kokkos::DefaultHostExecutionSpace::memory_space>
class DerivFieldSpan;

template <class ElementType, class SupportType, class LayoutStridedPolicy, class MemorySpace>
inline constexpr bool enable_deriv_field<
        DerivFieldSpan<ElementType, SupportType, LayoutStridedPolicy, MemorySpace>> = true;

template <class ElementType, class SupportType, class LayoutStridedPolicy, class MemorySpace>
inline constexpr bool enable_borrowed_deriv_field<
        DerivFieldSpan<ElementType, SupportType, LayoutStridedPolicy, MemorySpace>> = true;

namespace ddcHelper {

/**
 * @brief Copy the contents of one DerivFieldSpan into another.
 *
 * @param dst The DerivFieldSpan where the data will be saved.
 * @param src The DerivFieldSpan whose data will be copied.
 */
template <
        class FieldDst,
        class FieldSrc,
        std::enable_if_t<
                is_borrowed_deriv_field_v<FieldDst> && is_borrowed_deriv_field_v<FieldSrc>,
                bool> = true>
auto deepcopy(FieldDst&& dst, FieldSrc&& src)
{
    assert(dst.get_values_span().domain().extents() == src.get_values_span().domain().extents());

    DerivFieldSpan dst_span = dst.span_view();
    DerivFieldSpan src_span = src.span_view();

    dst_span.deepcopy(src_span);

    return dst_span;
}

/**
 * @brief Copy the contents of one DerivFieldSpan into another.
 *
 * @param execution_space The Kokkos execution space where the copy will be carried out.
 * @param dst The DerivFieldSpan where the data will be saved.
 * @param src The DerivFieldSpan whose data will be copied.
 */
template <class ExecSpace, class FieldDst, class FieldSrc>
auto deepcopy(ExecSpace const& execution_space, FieldDst&& dst, FieldSrc&& src)
{
    static_assert(is_borrowed_deriv_field_v<FieldDst>);
    static_assert(is_borrowed_deriv_field_v<FieldSrc>);

    assert(dst.get_values_span().domain().extents() == src.get_values_span().domain().extents());

    DerivFieldSpan dst_span = dst.span_view();
    DerivFieldSpan src_span = src.span_view();

    dst_span.deepcopy(execution_space, src_span);

    return dst_span;
}

} // namespace ddcHelper

/**
 * @brief A class which holds references to chunks of memory describing a field and its derivatives.
 *
 * The values of the field and the derivatives may be defined on different domains, but
 * the underlying mesh must be the same for both.
 *
 * @tparam ElementType The type of the elements inside the chunks.
 * @tparam ddc::DiscreteDomain<DDims...> The domain on which the internal chunks are defined.
 *          This domain is the physical domain on which the values are defined combined with
 *          the domain of the derivatives of interest (e.g. ddc::DiscreteDomain<Deriv<IDimX>, IDimX, IDimY>).
 * @tparam LayoutStridedPolicy The way in which the memory is laid out in memory (contiguous in
 *          the leading/trailing dimension, strided, etc).
 * @tparam MemorySpace The memory space where the data is saved (CPU/GPU).
 */
template <class ElementType, class... DDims, class LayoutStridedPolicy, class MemorySpace>
class DerivFieldSpan<ElementType, ddc::DiscreteDomain<DDims...>, LayoutStridedPolicy, MemorySpace>
    : public DerivFieldCommon<
              ddc::ChunkSpan<
                      ElementType,
                      ddc::DiscreteDomain<DDims...>,
                      LayoutStridedPolicy,
                      MemorySpace>,
              ddc::DiscreteDomain<DDims...>>
{
private:
    /// @brief The class from which this object inherits.
    using base_type = DerivFieldCommon<
            ddc::ChunkSpan<
                    ElementType,
                    ddc::DiscreteDomain<DDims...>,
                    LayoutStridedPolicy,
                    MemorySpace>,
            ddc::DiscreteDomain<DDims...>>;

public:
    /// @brief The type of the elements in the chunks.
    using element_type = typename base_type::element_type;

    /// @brief The DiscreteDomain on which the chunks in this object are defined.
    using discrete_domain_type = typename base_type::discrete_domain_type;

    /// @brief The DiscreteElement which can be used to index this object.
    using discrete_element_type = typename base_type::discrete_element_type;

    /// @brief A type sequence containing all derivatives present in this object.
    using deriv_tags = typename base_type::deriv_tags;

    /// @brief A type sequence containing all dimensions for which derivatives are present in this object.
    using physical_deriv_dims = typename detail::strip_deriv_t<deriv_tags>;

    /// @brief A type sequence containing all the physical dimensions on which the chunks are defined.
    using physical_dims = typename base_type::physical_dims;

    /// @brief The physical domain on which the field is defined.
    using physical_domain_type = typename base_type::physical_domain_type;

    /// @brief The type of the chunk stored in the array
    using chunk_type = typename base_type::chunk_type;

    /// type of a modifiable span of this span
    using span_type = DerivFieldSpan<
            ElementType,
            ddc::DiscreteDomain<DDims...>,
            LayoutStridedPolicy,
            MemorySpace>;

    /// type of a constant view of this span
    using view_type = DerivFieldSpan<
            ElementType const,
            ddc::DiscreteDomain<DDims...>,
            LayoutStridedPolicy,
            MemorySpace>;

private:
    /// @brief The DiscreteDomain which describes the derivatives present on each chunk.
    using discrete_deriv_domain_type = typename base_type::discrete_deriv_domain_type;

    /** @brief The DiscreteElement which describes the order of the derivatives in each dimension.
     * (e.g. second-order derivative).
     */
    using discrete_deriv_element_type = typename base_type::discrete_deriv_element_type;

    /// @brief The type of a reference to an element of the mdspan
    using reference = typename chunk_type::reference;

    /// @brief The number of chunks which must be created to describe this object.
    static constexpr int n_chunks = base_type::n_chunks;

private:
    /** @brief Get the subdomain to be extracted from a DerivField to build the internal_chunk at position ArrayIndex
     * along dimension Tag.
     */
    template <std::size_t ArrayIndex, class Tag>
    KOKKOS_FUNCTION constexpr ddc::DiscreteElement<Tag> get_chunk_subdomain_1d_idx()
    {
        // The Tag describes a dimension for which a derivative is defined
        if constexpr (ddc::in_tags_v<Tag, deriv_tags>) {
            // If the Chunk at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, deriv_tags>)) {
                return ddc::DiscreteElement<Tag>(1);
            }
            // If the Chunk at this index doesn't contain derivatives of this dimension
            else {
                return ddc::DiscreteElement<Tag>(0);
            }
        } else {
            // Empty DiscreteDomain to be discarded
            return ddc::DiscreteElement<Tag>();
        }
    }

    /// @brief Get the subdomain to be extracted from a DerivField to build the internal_chunk at position ArrayIndex
    template <std::size_t ArrayIndex>
    KOKKOS_FUNCTION constexpr discrete_deriv_element_type get_chunk_subdomain_idx()
    {
        return discrete_deriv_element_type(get_chunk_subdomain_1d_idx<ArrayIndex, DDims>()...);
    }

    /// @brief Initialise chunk spans inside internal_chunks.
    template <class Field, std::size_t... ArrayIndex>
    KOKKOS_FUNCTION constexpr void initialise_chunk_spans(
            Field const& chunks,
            std::index_sequence<ArrayIndex...>)
    {
        ((base_type::internal_chunks[ArrayIndex]
          = chunks.internal_chunks[chunks.get_array_index(get_chunk_subdomain_idx<ArrayIndex>())]),
         ...);
    }

    auto get_kokkos_view_from_internal_chunk(int index)
    {
        typename base_type::internal_mdspan_type chunk_span = base_type::internal_chunks[index];
        auto kokkos_layout = ddc::detail::build_kokkos_layout(
                chunk_span.extents(),
                chunk_span.mapping(),
                std::make_index_sequence<sizeof...(DDims)> {});
        return Kokkos::View<
                ddc::detail::mdspan_to_kokkos_element_t<ElementType, sizeof...(DDims)>,
                decltype(kokkos_layout),
                MemorySpace>(chunk_span.data_handle(), kokkos_layout);
    }

public:
    /**
     * @brief Copy-construct a DerivFieldSpan
     *
     * @param other The DerivFieldSpan being copied.
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivFieldSpan(DerivFieldSpan const& other) = default;

    /**
     * @brief Move-construct a DerivFieldSpan
     *
     * @param other The DerivFieldSpan being moved.
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivFieldSpan(DerivFieldSpan&& other) = default;

    /** 
     * @brief Constructs a new DerivFieldSpan containing a modifiable view on the data in a DerivField.
     *
     * @param field The DerivField to view.
     */
    template <
            class OElementType,
            int NDerivs,
            class Allocator,
            class = std::enable_if_t<std::is_same_v<typename Allocator::memory_space, MemorySpace>>>
    constexpr DerivFieldSpan(
            DerivField<OElementType, discrete_domain_type, NDerivs, Allocator>& field)
        : base_type(field.m_physical_domain, field.m_deriv_domain, field.m_cross_derivative_domain)
    {
        initialise_chunk_spans(field, std::make_integer_sequence<std::size_t, n_chunks> {});
    }

    /** 
     * @brief Constructs a new DerivFieldSpan containing a constant view on the data in a DerivField.
     *
     * @param field The DerivField to view.
     */
    // Disabled by SFINAE if `ElementType` is not `const` to avoid write access
    template <
            class OElementType,
            class SFINAEElementType = ElementType,
            class = std::enable_if_t<std::is_const_v<SFINAEElementType>>,
            int NDerivs,
            class Allocator,
            class = std::enable_if_t<std::is_same_v<typename Allocator::memory_space, MemorySpace>>>
    constexpr DerivFieldSpan(
            DerivField<OElementType, discrete_domain_type, NDerivs, Allocator> const& field)
        : base_type(field.m_physical_domain, field.m_deriv_domain, field.m_cross_derivative_domain)
    {
        initialise_chunk_spans(field, std::make_integer_sequence<std::size_t, n_chunks> {});
    }

    /**
     * @brief Copy construct a DerivFieldSpan. The element type may be changed to a complatible type.
     * (e.g. double -> const double).
     *
     * @param field The DerivFieldSpan to be copied.
     */
    template <class OElementType>
    KOKKOS_FUNCTION constexpr DerivFieldSpan(DerivFieldSpan<
                                             OElementType,
                                             discrete_domain_type,
                                             LayoutStridedPolicy,
                                             MemorySpace> const& field)
        : base_type(field.m_physical_domain, field.m_deriv_domain, field.m_cross_derivative_domain)
    {
        initialise_chunk_spans(field, std::make_integer_sequence<std::size_t, n_chunks> {});
    }

    KOKKOS_DEFAULTED_FUNCTION ~DerivFieldSpan() = default;

    /** Copy-assigns a new value to this DerivFieldSpan, yields a new view to the same data
     * @param other the DerivFieldSpan to copy
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivFieldSpan& operator=(DerivFieldSpan const& other)
            = default;

    /** Move-assigns a new value to this DerivFieldSpan
     * @param other the DerivFieldSpan to move
     * @return *this
     */
    KOKKOS_DEFAULTED_FUNCTION constexpr DerivFieldSpan& operator=(DerivFieldSpan&& other) = default;

    /**
     * @brief Copy the source DerivFieldSpan into this DerivFieldSpan using ddc::parallel_deepcopy. 
     *
     * @param src The DerivFieldSpan containing the data to be copied.
     */
    template <class OElementType, class OLayoutStridedPolicy, class OMemorySpace>
    void deepcopy(
            DerivFieldSpan<OElementType, discrete_domain_type, OLayoutStridedPolicy, OMemorySpace>
                    src)
    {
        for (int i(0); i < n_chunks; ++i) {
            auto kokkos_span = get_kokkos_view_from_internal_chunk(i);
            auto src_kokkos_span = src.get_kokkos_view_from_internal_chunk(i);
            Kokkos::deep_copy(kokkos_span, src_kokkos_span);
        }
    }

    /**
     * @brief Copy the source DerivFieldSpan into this DerivFieldSpan using ddc::parallel_deepcopy.
     *
     * @param execution_space The execution space on which the copy will be carried out.
     * @param src The DerivFieldSpan containing the data to be copied.
     */
    template <class ExecSpace, class OElementType, class OLayoutStridedPolicy, class OMemorySpace>
    void deepcopy(
            ExecSpace const& execution_space,
            DerivFieldSpan<OElementType, discrete_domain_type, OLayoutStridedPolicy, OMemorySpace>
                    src)
    {
        for (int i(0); i < n_chunks; ++i) {
            auto kokkos_span = get_kokkos_view_from_internal_chunk(i);
            auto src_kokkos_span = src.get_kokkos_view_from_internal_chunk(i);
            Kokkos::deep_copy(execution_space, kokkos_span, src_kokkos_span);
        }
    }

    /**
     * @brief Get an element from a constant field.
     * A DiscreteElement describes the element of interest. If information about the derivatives is
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
        using full_element_type = ddcHelper::combine_t<DElem...>;
        full_element_type elem(elems...);
        return base_type::get_internal_chunk(elem)();
    }

    /**
     * @brief Get a constant DerivFieldSpan of this field.
     *
     * @returns A constant span of this field.
     */
    KOKKOS_FUNCTION constexpr view_type span_cview() const
    {
        return view_type(*this);
    }

    /**
     * @brief Get a modifiable DerivFieldSpan of this field.
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
using DerivFieldView
        = DerivFieldSpan<ElementType const, SupportType, LayoutStridedPolicy, MemorySpace>;

namespace detail {
template <
        class NewMemorySpace,
        class ElementType,
        class SupportType,
        class Layout,
        class MemorySpace>
struct OnMemorySpace<NewMemorySpace, DerivFieldSpan<ElementType, SupportType, Layout, MemorySpace>>
{
    using type = DerivFieldSpan<ElementType, SupportType, Layout, NewMemorySpace>;
};
}; // namespace detail
