// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

#include "derivative_field_common.hpp"
#include "derivative_field_span.hpp"

template <
        class ElementType,
        class Domain,
        int NDerivs,
        class Allocator = ddc::HostAllocator<ElementType>>
class DerivField;

template <class ElementType, class SupportType, int NDerivs, class Allocator>
inline constexpr bool
        enable_deriv_field<DerivField<ElementType, SupportType, NDerivs, Allocator>> = true;


/**
 * @brief A class which holds a chunk of memory describing a field and its derivatives.
 *
 * The values of the field and the derivatives may be defined on different domains, but
 * the underlying mesh must be the same for both.
 *
 * @tparam ElementType The type of the elements inside the chunks.
 * @tparam ddc::DiscreteDomain<DDims...> The domain on which the internal chunks are defined.
 *          This domain is the physical domain on which the values are defined combined with
 *          the domain of the derivatives of interest (e.g. ddc::DiscreteDomain<Deriv<IDimX>, IDimX, IDimY>).
 * @tparam NDerivs The number of derivatives which are defined in the dimensions where derivatives
 *          appear.
 * @tparam Allocator The Allocator which is used to build the Chunk. This provides information
 *          about the memory space where the data will be saved.
 */
template <class ElementType, class... DDims, int NDerivs, class Allocator>
class DerivField<ElementType, ddc::DiscreteDomain<DDims...>, NDerivs, Allocator>
    : public DerivFieldCommon<
              ddc::Chunk<ElementType, ddc::DiscreteDomain<DDims...>, Allocator>,
              ddc::DiscreteDomain<DDims...>>
{
private:
    /// @brief The class from which this object inherits.
    using base_type = DerivFieldCommon<
            ddc::Chunk<ElementType, ddc::DiscreteDomain<DDims...>, Allocator>,
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
    using physical_deriv_dims = typename base_type::physical_deriv_dims;

    /// @brief A type sequence containing all the physical dimensions on which the chunks are defined.
    using physical_dims = typename base_type::physical_dims;

    /// @brief The physical domain on which the field is defined.
    using physical_domain_type = typename base_type::physical_domain_type;

    /// @brief The type of the chunk stored in the array
    using chunk_type = typename base_type::chunk_type;

    template <class, class, class, class>
    friend class DerivSpan;

    /// type of a modifiable span of this span
    using span_type = DerivFieldSpan<
            ElementType,
            ddc::DiscreteDomain<DDims...>,
            typename chunk_type::layout_type,
            typename chunk_type::memory_space>;

    /// type of a constant view of this span
    using view_type = DerivFieldSpan<
            ElementType const,
            ddc::DiscreteDomain<DDims...>,
            typename chunk_type::layout_type,
            typename chunk_type::memory_space>;

private:
    /** @brief A domain describing all the domains where derivatives are defined.
     *
     * A domain describing all the domains where derivatives are defined. If this domain is
     * sliced into a 1D object then it describes the domain on which the derivative in that
     * direction is defined. This object is used internally to initialise the chunks.
     */
    using physical_deriv_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain<physical_deriv_dims>;

    /// @brief The DiscreteDomain which describes the derivatives present on each chunk.
    using discrete_deriv_domain_type = typename base_type::discrete_deriv_domain_type;

    /** @brief The DiscreteElement which describes the order of the derivatives in each dimension.
     * (e.g. second-order derivative).
     */
    using discrete_deriv_element_type = typename base_type::discrete_deriv_element_type;

    using mapping_type = typename chunk_type::mapping_type;

    using extents_type = typename chunk_type::extents_type;

    using allocator_type = ddc::KokkosAllocator<element_type, typename chunk_type::memory_space>;

    using internal_mdspan_type = typename base_type::internal_mdspan_type;

    /// @brief The number of chunks which must be created to describe this object.
    static constexpr int n_chunks = base_type::n_chunks;

private:
    /// @brief A function to get the domain along direction Tag for the ArrayIndex-th element of internal_chunks.
    template <std::size_t ArrayIndex, class Tag>
    ddc::DiscreteDomain<Tag> get_domain(
            physical_domain_type val_domain,
            physical_deriv_domain_type deriv_domains)
    {
        // If the Tag describes a dimension for which a derivative is defined
        if constexpr (ddc::in_tags_v<Tag, physical_deriv_dims>) {
            // If the Chunk at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, physical_deriv_dims>)) {
                return ddc::select<Tag>(deriv_domains);
            }
            // If the Chunk at this index doesn't contain derivatives of this dimension
            else {
                return ddc::select<Tag>(val_domain);
            }
        }
        // If the Tag describes a derivative
        else if constexpr (ddc::in_tags_v<Tag, deriv_tags>) {
            // If the Chunk at this index contains the derivatives of this dimension
            if constexpr (ArrayIndex & (1 << ddc::type_seq_rank_v<Tag, deriv_tags>)) {
                return ddc::DiscreteDomain<
                        Tag>(ddc::DiscreteElement<Tag>(1), ddc::DiscreteVector<Tag>(NDerivs));
            }
            // If the Chunk at this index doesn't contain derivatives of this dimension
            else {
                return ddc::DiscreteDomain<
                        Tag>(ddc::DiscreteElement<Tag>(0), ddc::DiscreteVector<Tag>(1));
            }
        }
        // If the Tag describes a dimension for which derivatives are not defined.
        else {
            return ddc::select<Tag>(val_domain);
        }
    }

    /// @brief Get the domain of the Chunk at the index ArrayIndex in internal_chunks.
    template <std::size_t ArrayIndex>
    ddc::DiscreteDomain<DDims...> get_chunk_domain(
            physical_domain_type val_domain,
            physical_deriv_domain_type deriv_domains)
    {
        return ddc::DiscreteDomain<DDims...>(
                get_domain<ArrayIndex, DDims>(val_domain, deriv_domains)...);
    }

    /// @brief Get the size of the mdspan in dimension QueryDDim, at the index ArrayIndex.
    template <class QueryDDim, std::size_t ArrayIndex>
    std::size_t get_span_size()
    {
        if constexpr (ddc::in_tags_v<QueryDDim, deriv_tags>) {
            if (ArrayIndex & (1 << ddc::type_seq_rank_v<QueryDDim, deriv_tags>)) {
                return NDerivs;
            } else {
                return 1;
            }
        } else if constexpr (ddc::in_tags_v<QueryDDim, physical_deriv_dims>) {
            if constexpr (
                    ArrayIndex & (1 << ddc::type_seq_rank_v<ddc::Deriv<QueryDDim>, deriv_tags>)) {
                DiscreteSubDomain<QueryDDim> local_dom(base_type::m_cross_derivative_domain);
                return local_dom.extents().value();
            } else {
                ddc::DiscreteDomain<QueryDDim> local_dom(base_type::m_physical_domain);
                return local_dom.extents().value();
            }
        } else {
            ddc::DiscreteDomain<QueryDDim> local_dom(base_type::m_physical_domain);
            return local_dom.extents().value();
        }
    }

    /// @brief Make the internal mdspan that will be saved in internal_chunks at the index ArrayIndex.
    template <std::size_t ArrayIndex>
    std::enable_if_t<std::is_constructible_v<mapping_type, extents_type>, internal_mdspan_type>
    make_internal_mdspan(
            ddc::KokkosAllocator<element_type, typename chunk_type::memory_space> allocator)
    {
        std::size_t alloc_size(((get_span_size<DDims, ArrayIndex>()) * ...));
        element_type* ptr = allocator.allocate("", alloc_size);

        extents_type extents_r(get_span_size<DDims, ArrayIndex>()...);
        mapping_type mapping_r(extents_r);

        std::array<std::size_t, sizeof...(DDims)> strides_s {
                mapping_r.stride(ddc::type_seq_rank_v<DDims, ddc::detail::TypeSeq<DDims...>>)...};
        std::experimental::layout_stride::mapping<extents_type> mapping_s(extents_r, strides_s);
        return internal_mdspan_type(ptr, mapping_s);
    }

    /// @brief Initialise the chunks inside internal_chunks.
    template <std::size_t... ArrayIndex>
    void initialise_chunks(allocator_type allocator, std::index_sequence<ArrayIndex...>)
    {
        ((base_type::internal_chunks[ArrayIndex] = make_internal_mdspan<ArrayIndex>(allocator)),
         ...);
    }

public:
    /**
     * @brief The constructor for DerivField. The constructor initialises the chunks using
     * the provided domains.
     *
     * @param val_domain The domain on which the values of the field are defined.
     * @param m_deriv_domain The 1D sub-domains on which the derivatives of the field are defined.
     */
    template <class... DerivDoms>
    DerivField(physical_domain_type val_domain, DiscreteSubDomain<DerivDoms>... m_deriv_domain)
        : DerivField(allocator_type(), val_domain, m_deriv_domain...)
    {
    }

    /**
     * @brief The constructor for DerivField. The constructor initialises the chunks using
     * the provided domains.
     *
     * @param allocator The object which allocates the memory on the CPU or GPU.
     * @param val_domain The domain on which the values of the field are defined.
     * @param m_deriv_domain The 1D sub-domains on which the derivatives of the field are defined.
     */
    template <class... DerivDoms>
    DerivField(
            allocator_type allocator,
            physical_domain_type val_domain,
            DiscreteSubDomain<DerivDoms>... m_deriv_domain)
        : base_type(
                val_domain,
                discrete_deriv_domain_type(ddc::DiscreteDomain<ddc::Deriv<DerivDoms>>(
                        ddc::DiscreteElement<ddc::Deriv<DerivDoms>>(1),
                        ddc::DiscreteVector<ddc::Deriv<DerivDoms>>(NDerivs))...),
                to_subdomain_collection<physical_deriv_dims>(m_deriv_domain...))
    {
        static_assert(
                ddc::type_seq_same_v<ddc::detail::TypeSeq<DerivDoms...>, physical_deriv_dims>);
        initialise_chunks(allocator, std::make_integer_sequence<std::size_t, n_chunks> {});
    }

    /// Defaulted destructor
    ~DerivField() = default;

    /// Deleted copy operator
    DerivField& operator=(DerivField const& other) = delete;

    /** Move-assigns a new value to this field
     * @param other the Chunk to move
     * @return *this
     */
    DerivField& operator=(DerivField&& other) = default;

    /**
     * @brief Get a modifiable reference to an element from a constant field.
     * A DiscreteElement describes the element of interest. If information about the derivatives is
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
        using full_element_type = ddcHelper::combine_t<DElem...>;
        full_element_type elem(elems...);
        return base_type::get_internal_chunk(elem)();
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
    element_type const& operator()(DElem... elems) const noexcept
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
    view_type span_cview() const
    {
        return view_type(*this);
    }

    /**
     * @brief Get a constant DerivFieldSpan of this field.
     *
     * @returns A constant span of this field.
     */
    view_type span_view() const
    {
        return view_type(*this);
    }

    /**
     * @brief Get a modifiable DerivFieldSpan of this field.
     *
     * @returns A span of this field.
     */
    span_type span_view()
    {
        return *this;
    }
};

namespace detail {
template <class NewMemorySpace, class ElementType, class SupportType, int NDerivs, class Allocator>
struct OnMemorySpace<NewMemorySpace, DerivField<ElementType, SupportType, NDerivs, Allocator>>
{
    using type = DerivField<
            ElementType,
            SupportType,
            NDerivs,
            ddc::KokkosAllocator<ElementType, NewMemorySpace>>;
};
} // namespace detail
