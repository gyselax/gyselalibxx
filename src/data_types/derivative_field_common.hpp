// SPDX-License-Identifier: MIT

#pragma once

#include <array>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp> // Needed for ddc::Deriv

#include <sll/math_tools.hpp>

#include <ddc_helper.hpp>

#include "deriv_details.hpp"
#include "discrete_subdomain.hpp"

template <class T>
inline constexpr bool enable_deriv_field = false;

template <class T>
inline constexpr bool enable_borrowed_deriv_field = false;

template <class T>
inline constexpr bool is_deriv_field_v
        = enable_deriv_field<std::remove_const_t<std::remove_reference_t<T>>>;

template <class T>
inline constexpr bool is_borrowed_deriv_field_v
        = (std::is_lvalue_reference_v<T>)
          || (enable_borrowed_deriv_field<std::remove_cv_t<std::remove_reference_t<T>>>);

template <class ChunkType, class SupportType>
class DerivFieldCommon;

/**
 * @brief An abstract class which holds a chunk of memory describing a field and its derivatives.
 * This is the superclass for DerivField and DerivSpan.
 *
 * @tparam ChunkType The type of the object stored in the internal_chunks array.
 * @tparam ddc::DiscreteDomain<DDims...> The domain on which the internal chunks are defined.
 *          This domain is the physical domain on which the values are defined combined with
 *          the domain of the derivatives of interest (e.g. ddc::DiscreteDomain<Deriv<IDimX>, IDimX, IDimY>).
 */
template <class ChunkType, class... DDims>
class DerivFieldCommon<ChunkType, ddc::DiscreteDomain<DDims...>>
{
public:
    /// @brief The type of the chunk stored in the array
    using chunk_type = ChunkType;

    /// @brief A type sequence containing all derivatives present in this object.
    using deriv_tags = detail::deriv_sub_set_t<ddc::detail::TypeSeq<DDims...>>;

    /// @brief A type sequence containing all physical dimensions for which derivatives are present in this object.
    using physical_deriv_dims = typename detail::strip_deriv_t<deriv_tags>;

    /// @brief A type sequence containing all the physical dimensions on which the chunks are defined.
    using physical_dims = ddc::type_seq_remove_t<ddc::detail::TypeSeq<DDims...>, deriv_tags>;

    /// @brief The type of the elements in the chunks.
    using element_type = typename chunk_type::element_type;

    /// @brief The DiscreteDomain on which the chunks in this object are defined.
    using discrete_domain_type = ddc::DiscreteDomain<DDims...>;

    /// @brief The DiscreteElement which can be used to index this object.
    using discrete_element_type = ddc::DiscreteElement<DDims...>;

protected:
    /// @brief The type of the memory block stored in the array internal_chunks
    using internal_mdspan_type = std::experimental::mdspan<
            element_type,
            std::experimental::dextents<std::size_t, sizeof...(DDims)>,
            std::experimental::layout_stride>;

    /// @brief The type of a constant view on the memory block stored in the array internal_chunks
    using internal_mdview_type = std::experimental::mdspan<
            const element_type,
            std::experimental::dextents<std::size_t, sizeof...(DDims)>,
            std::experimental::layout_stride>;

    /// @brief The type of the span of the chunk stored in the array
    using chunk_span = typename ChunkType::span_type;

    /// @brief The type of the view of the chunk stored in the array
    using chunk_view = typename ChunkType::view_type;

    /// @brief The domain for the chunk excluding derivatives
    using physical_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<physical_dims>;

    /// @brief The DiscreteElement which describes the physical position where values are defined.
    using physical_element_type = typename physical_domain_type::discrete_element_type;

    /// @brief The DiscreteDomain which describes the derivatives present on each chunk.
    using discrete_deriv_domain_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<deriv_tags>;

    /** @brief The DiscreteElement which describes the order of the derivatives in each dimension.
     * (e.g. second-order derivative).
     */
    using discrete_deriv_element_type = typename discrete_deriv_domain_type::discrete_element_type;

    /** @brief The DiscreteElement which describes the order of the derivatives in each dimension.
     * (e.g. second-order derivative).
     */
    using discrete_deriv_vector_type = typename discrete_deriv_domain_type::mlength_type;

    /// @brief The number of chunks which must be created to describe this object.
    static constexpr int n_chunks = 1 << ddc::type_seq_size_v<deriv_tags>;

    template <class, class, int, class>
    friend class DerivField;

    template <class, class, class, class>
    friend class DerivFieldSpan;

protected:
    /** @brief The internal chunks describing the values and derivatives.
     *
     * The chunks which contain the values have different domains to the chunks containing derivatives
     * so a DDC object cannot be used directly.
     * E.g. for a 2D field (X,Y) with derivatives provided in both directions the elements of internal_chunks
     * have the type : Chunk<double, DiscreteDomain<Deriv<IDimX>, Deriv<IDimY>, IDimX, IDimY>
     * The derivative domains are then defined such that the  elements of internal_chunks represent:
     * 0 : @f$f(x,y)@f$
     * 1 : @f$\partial_x^k f(x,y)  \quad \forall 1 \leq k \leq NDerivs@f$
     * 2 : @f$\partial_y^k f(x,y)  \quad \forall 1 \leq k \leq NDerivs@f$
     * 3 : @f$\partial_x^j \partial_y^k f(x,y)  \quad \forall 1 \leq j \leq NDerivs,  \forall 1 \leq k \leq NDerivs@f$
     */
    std::array<internal_mdspan_type, n_chunks> internal_chunks;

    /// @brief The physical domain on which the values are defined.
    physical_domain_type m_physical_domain;

    /// @brief The domain of available derivatives.
    discrete_deriv_domain_type m_deriv_domain;

    /// @brief The physical domains on which the derivatives are defined.
    to_subdomain_collection<physical_deriv_dims> m_cross_derivative_domain;

protected:
    /**
     * @brief An internal function which provides the index of an element inside the internal_chunks.
     * A DiscreteElement describes the element of interest. If information about the derivatives is
     * missing then it is assumed that the 0-th order derivative is requested.
     *
     * @param elem The element of interest.
     *
     * @returns int The index of the internal chunk inside the array internal_chunks.
     * @returns discrete_element_type The index of the element of interest inside the chunk of interest.
     */
    template <class DElem>
    KOKKOS_FUNCTION std::pair<int, discrete_element_type> get_index(DElem elem) const
    {
        discrete_deriv_element_type default_derivatives
                = detail::no_derivative_element<deriv_tags>();
        discrete_deriv_element_type deriv_index(detail::select_default(elem, default_derivatives));
        physical_element_type physical_index(elem);
        discrete_element_type index(physical_index, deriv_index);
        return std::pair<int, discrete_element_type>(get_array_index(deriv_index), index);
    }

    /**
     * @brief An internal function which provides the index of a chunk inside the internal_chunks array.
     * A DiscreteElement describes the derivatives of interest. n-th order derivatives are stored in the
     * same chunk for all n!=0 so it is sufficient to provide any valid element from the derivatives.
     *
     * The index is calculated using a bit mask.
     * Mathematically the equation which determines the index in `internal_chunks` is:
     * @f$ \sum_{\xi} 2^{i_\xi} \delta_{\text{using derivs in }\xi} @f$
     * where @f$i_\xi@f$ is the index of the dimension @f$\xi@f$ in the tags.
     *
     * @param idx The derivatives of interest.
     *
     * @returns int The index of the internal chunk inside the array internal_chunks.
     * @returns discrete_deriv_domain_type The domain of the derivatives at the chunk at the index.
     */
    template <class... Tag>
    KOKKOS_FUNCTION int get_array_index(ddc::DiscreteElement<Tag...> idx) const
    {
        static_assert(std::is_same_v<ddc::DiscreteElement<Tag...>, discrete_deriv_element_type>);
        return (0 + ...
                + (int(ddc::select<Tag>(idx) != ddc::DiscreteElement<Tag>(0))
                   << ddc::type_seq_rank_v<Tag, deriv_tags>));
    }

    /**
     * @brief Get an object which can be used to slice an mdspan.
     *
     * @param slice_idx The DDC element which should be used to slice the span.
     * @param array_idx The index of the mdspan in internal_chunks that will be sliced.
     *
     * @tparam QueryDDim The dimension along which we want to slice.
     *
     * @returns An index or a slice which can be used to slice an mdspan.
     */
    template <class QueryDDim, class... ODDims>
    KOKKOS_FUNCTION constexpr auto get_slicer_for(
            ddc::DiscreteElement<ODDims...> const& slice_idx,
            int array_idx) const
    {
        if constexpr (!ddc::in_tags_v<QueryDDim, ddc::detail::TypeSeq<ODDims...>>) {
            return std::experimental::full_extent;
        } else {
            if constexpr (ddc::in_tags_v<QueryDDim, physical_deriv_dims>) {
                // Physical dimension along which derivatives are known
                // If information is available about the physical domain
                if (array_idx & (1 << ddc::type_seq_rank_v<ddc::Deriv<QueryDDim>, deriv_tags>)) {
                    // If the derivative is being requested
                    return m_cross_derivative_domain.get_index(ddc::select<QueryDDim>(slice_idx));
                }
            }
            if constexpr (ddc::in_tags_v<QueryDDim, physical_dims>) {
                // Physical dimension along which derivatives are not known
                return std::size_t(
                        (ddc::select<QueryDDim>(slice_idx)
                         - ddc::select<QueryDDim>(m_physical_domain).front()));
            }
            if constexpr (ddc::in_tags_v<QueryDDim, deriv_tags>) {
                // Derivative dimension
                if (array_idx & (1 << ddc::type_seq_rank_v<QueryDDim, deriv_tags>)) {
                    // If array contains derivatives
                    return std::size_t(
                            (ddc::select<QueryDDim>(slice_idx)
                             - ddc::DiscreteElement<QueryDDim>(1)));
                } else {
                    // If array doesn't contain derivatives
                    assert(ddc::select<QueryDDim>(slice_idx) == ddc::DiscreteElement<QueryDDim>(0));
                    return std::size_t(0);
                }
            }
        }
    }

    /**
     * @brief Get an object which can be used to slice an mdspan.
     *
     * @param slice_domain The DDC domain which should be used to slice the span.
     * @param array_idx The index of the mdspan in internal_chunks that will be sliced.
     *
     * @tparam QueryDDim The dimension along which the we want to slice.
     *
     * @returns A slice (often in the form of a (start, end) pair) which can be used to slice an mdspan.
     */
    template <class QueryDDim, class... ODDims>
    KOKKOS_FUNCTION constexpr auto get_slicer_for(
            ddc::DiscreteDomain<ODDims...> const& slice_domain,
            int array_idx) const
    {
        if constexpr (!ddc::in_tags_v<QueryDDim, ddc::detail::TypeSeq<ODDims...>>) {
            return std::experimental::full_extent;
        } else {
            if constexpr (ddc::in_tags_v<QueryDDim, physical_deriv_dims>) {
                // Physical dimension along which derivatives are known
                // If information is available about the physical domain
                ddc::DiscreteDomain<QueryDDim> requested_dom(slice_domain);
                if (array_idx & (1 << ddc::type_seq_rank_v<ddc::Deriv<QueryDDim>, deriv_tags>)) {
                    // If the derivative is being requested
                    assert(m_cross_derivative_domain.contains(requested_dom));
                    return std::pair<std::size_t, std::size_t>(
                            m_cross_derivative_domain.get_index(requested_dom.front()),
                            m_cross_derivative_domain.get_index(requested_dom.back()) + 1);
                }
            }
            if constexpr (ddc::in_tags_v<QueryDDim, physical_dims>) {
                // Physical dimension along which derivatives are not known
                return std::pair<std::size_t, std::size_t>(
                        ddc::select<QueryDDim>(slice_domain).front()
                                - ddc::select<QueryDDim>(m_physical_domain).front(),
                        ddc::select<QueryDDim>(slice_domain).back() + 1
                                - ddc::select<QueryDDim>(m_physical_domain).front());
            }
            if constexpr (ddc::in_tags_v<QueryDDim, deriv_tags>) {
                // Derivative dimension
                if (array_idx & (1 << ddc::type_seq_rank_v<QueryDDim, deriv_tags>)) {
                    // If array contains derivatives
                    return std::pair<std::size_t, std::size_t>(
                            ddc::select<QueryDDim>(slice_domain).front()
                                    - ddc::DiscreteElement<QueryDDim>(1),
                            ddc::select<QueryDDim>(slice_domain).back() + 1
                                    - ddc::DiscreteElement<QueryDDim>(1));
                } else {
                    // If array doesn't contain derivatives
                    assert(ddc::select<QueryDDim>(slice_domain).front()
                           == ddc::DiscreteElement<QueryDDim>(0));
                    assert(ddc::select<QueryDDim>(slice_domain).back()
                           == ddc::select<QueryDDim>(slice_domain).front());
                    return std::pair<std::size_t, std::size_t>(0, 1);
                }
            }
        }
    }

    /**
     * Get a ddc::ChunkSpan from a subset of one of the mdspans in internal_chunks. The provided
     * DiscreteDomain is used to slice the mdspan in such a way that the resulting mdspan can be
     * saved in a ddc::ChunkSpan. This means that all information about derivatives must be provided.
     * If information about a derivative is missing then it is assumed that the 0-th order derivative
     * is requested.
     *
     * @param dom The domain used to slice the mdspan.
     *
     * @returns ChunkSpan The subset of the internal mdspan.
     */
    template <class... ODims>
    KOKKOS_FUNCTION auto get_internal_chunk(ddc::DiscreteDomain<ODims...> dom) const
    {
        // Get the types related to the provided information
        using provided_tags = ddc::detail::TypeSeq<ODims...>;
        using provided_deriv_tags = detail::deriv_sub_set_t<provided_tags>;

        // Get the types related to the implicit information
        using remaining_deriv_tags = ddc::type_seq_remove_t<deriv_tags, provided_deriv_tags>;
        using remaining_deriv_domain_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<remaining_deriv_tags>;

        // Find the domain of the derivatives (either provided or a domain containing only the 0-th derivative)
        remaining_deriv_domain_type no_deriv_domain = detail::get_domain_from_element(
                detail::no_derivative_element<remaining_deriv_tags>());
        discrete_deriv_domain_type deriv_domain(dom, no_deriv_domain);

        // Find the physical domain of the chunk
        physical_domain_type local_physical_domain(detail::select_default(dom, m_physical_domain));

        // Find the discrete domain of the chunk
        discrete_domain_type full_domain(local_physical_domain, deriv_domain);

        // Find the index of the internal chunk
        int const array_idx = get_array_index(deriv_domain.front());

        // Get the relevant internal chunk
        internal_mdspan_type internal_view = internal_chunks[array_idx];
        // Slice the relevant section of the internal chunk
        auto subview = std::experimental::
                submdspan(internal_view, get_slicer_for<DDims>(full_domain, array_idx)...);
        // Create a ddc::ChunkSpan with the expected domain
        ddc::ChunkSpan<
                element_type,
                discrete_domain_type,
                typename decltype(subview)::layout_type,
                typename chunk_type::memory_space>
                local_chunk(subview, full_domain);

        // If necessary, slice off the derivative dimensions deduced implicitly
        if constexpr (ddc::type_seq_size_v<remaining_deriv_tags> == 0) {
            return local_chunk;
        } else {
            return local_chunk[no_deriv_domain.front()];
        }
    }

    /**
     * Get a ddc::ChunkSpan from a subset of one of the mdspans in internal_chunks. The provided
     * DiscreteElement is used to slice the mdspan in such a way that the resulting mdspan can be
     * saved in a ddc::ChunkSpan. This means that all information about derivatives must be provided.
     * If information about a derivative is missing then it is assumed that the 0-th order derivative
     * is requested.
     *
     * @param elem The element used to slice the mdspan.
     *
     * @returns ChunkSpan The subset of the internal mdspan.
     */
    template <class... ODims>
    KOKKOS_FUNCTION auto get_internal_chunk(ddc::DiscreteElement<ODims...> elem) const
    {
        // Get the types related to the provided information
        using provided_tags = ddc::detail::TypeSeq<ODims...>;
        using provided_deriv_tags = detail::deriv_sub_set_t<provided_tags>;
        using provided_physical_tags = ddc::type_seq_remove_t<provided_tags, provided_deriv_tags>;
        using provided_deriv_domain_type
                = ddc::detail::convert_type_seq_to_discrete_domain_t<provided_deriv_tags>;
        using provided_deriv_element_type =
                typename provided_deriv_domain_type::discrete_element_type;

        // Get the types related to the implicit information
        using remaining_deriv_tags = ddc::type_seq_remove_t<deriv_tags, provided_deriv_tags>;
        using remaining_deriv_domain_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<remaining_deriv_tags>;
        using remaining_deriv_element_type =
                typename remaining_deriv_domain_type::discrete_element_type;

        // Get the types related to the final chunk type
        using sliced_tags = ddc::type_seq_merge_t<provided_physical_tags, deriv_tags>;
        using sliced_domain_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<sliced_tags>;
        using sliced_element_type = typename sliced_domain_type::discrete_element_type;
        using final_tags = ddc::type_seq_remove_t<ddc::detail::TypeSeq<DDims...>, sliced_tags>;
        using final_domain_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<final_tags>;

        // Get the index of the relevant derivatives
        provided_deriv_element_type requested_derivs(elem);
        remaining_deriv_element_type no_deriv_elements(
                detail::no_derivative_element<remaining_deriv_tags>());
        discrete_deriv_element_type deriv_index(requested_derivs, no_deriv_elements);

        // Get the element which will slice the mdspan
        sliced_element_type slice_idx(elem, no_deriv_elements);

        // Find the index of the internal chunk
        int const array_idx = get_array_index(deriv_index);

        // Get the final domain
        final_domain_type final_domain(m_physical_domain);

        // Get the relevant internal chunk
        internal_mdspan_type internal_view = internal_chunks[array_idx];
        // Slice the relevant section of the internal chunk
        auto subview = std::experimental::
                submdspan(internal_view, get_slicer_for<DDims>(slice_idx, array_idx)...);
        // Create a ddc::ChunkSpan with the expected domain
        ddc::ChunkSpan<
                element_type,
                final_domain_type,
                typename decltype(subview)::layout_type,
                typename chunk_type::memory_space>
                local_chunk(subview, final_domain);

        return local_chunk;
    }

    /**
     * @brief Protected constructor to be used by subclasses to initialise domains.
     *
     * @param physical_domain The domain on which the values of the function are defined.
     * @param deriv_domain The domain of the provided derivatives.
     * @param cross_derivative_domain The cross product of the domains on which the derivatives
     *              of the function are defined.
     */
    KOKKOS_FUNCTION DerivFieldCommon(
            physical_domain_type physical_domain,
            discrete_deriv_domain_type deriv_domain,
            to_subdomain_collection<physical_deriv_dims> cross_derivative_domain)
        : m_physical_domain(physical_domain)
        , m_deriv_domain(deriv_domain)
        , m_cross_derivative_domain(cross_derivative_domain)
    {
    }

public:
    KOKKOS_DEFAULTED_FUNCTION ~DerivFieldCommon() = default;

    /**
     * @brief Get a ChunkView describing a subset of the data.
     *
     * @param slice_spec A discrete element describing the position at which these dimensions should be
     *          indexed. If information about the derivatives is missing then it is assumed that the
     *          0-th order derivative is requested.
     *
     * @returns ChunkView A subset of the data.
     */
    template <class... QueryDDims>
    constexpr auto operator[](ddc::DiscreteElement<QueryDDims...> const& slice_spec) const
    {
        return get_internal_chunk(slice_spec).span_cview();
    }

    /**
     * @brief Get a ChunkSpan describing a subset of the data.
     *
     * @param slice_spec A discrete element describing the position at which these dimensions should be
     *          indexed. If information about the derivatives is missing then it is assumed that the
     *          0-th order derivative is requested.
     *
     * @returns ChunkSpan A subset of the data.
     */
    template <class... QueryDDims>
    constexpr auto operator[](ddc::DiscreteElement<QueryDDims...> const& slice_spec)
    {
        return get_internal_chunk(slice_spec);
    }

    /**
     * @brief Get a ChunkSpan describing a subset of the data.
     * This function allows a slice to be obtained however it is designed to return a ChunkSpan. It is
     * therefore not possible to request data from multiple chunks (e.g. derivatives from 0 to 3).
     *
     * @param odomain A discrete domain describing the position at which these dimensions should be
     *          indexed. If information about the derivatives is missing then it is assumed that the
     *          0-th order derivative is requested.
     *
     * @returns ChunkSpan A subset of the data.
     */
    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr auto operator[](ddc::DiscreteDomain<QueryDDims...> const& odomain)
    {
        return get_internal_chunk(odomain);
    }

    /**
     * @brief Get a ChunkView describing a subset of the data.
     * This function allows a slice to be obtained however it is designed to return a ChunkView. It is
     * therefore not possible to request data from multiple chunks (e.g. derivatives from 0 to 3).
     *
     * @param odomain A discrete domain describing the position at which these dimensions should be
     *          indexed. If information about the derivatives is missing then it is assumed that the
     *          0-th order derivative is requested.
     *
     * @returns ChunkView A subset of the data.
     */
    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr auto operator[](
            ddc::DiscreteDomain<QueryDDims...> const& odomain) const
    {
        return get_internal_chunk(odomain).span_cview();
    }

    /**
     * @brief Get one of the mdspans from the internal array internal_chunks.
     * This function takes domains on the derivative directions. Where derivatives are
     * missing it is assumed that the 0-th order derivative is requested. This dimension
     * is stripped from the resulting span. This is the recommended way to access the
     * internal chunks.
     *
     * @param provided_deriv_domain The derivative domain which should be retained.
     *
     * @returns ChunkSpan The chunk span on the physical domain and the requested domains.
     */
    template <class... ODims>
    auto get_mdspan(ddc::DiscreteDomain<ODims...> provided_deriv_domain)
    {
        static_assert(((ddc::in_tags_v<ODims, deriv_tags>)&&...));
        using provided_deriv_tags = ddc::detail::TypeSeq<ODims...>;
        using remaining_deriv_tags = ddc::type_seq_remove_t<deriv_tags, provided_deriv_tags>;
        using remaining_deriv_domain_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<remaining_deriv_tags>;
        using remaining_deriv_element_type =
                typename remaining_deriv_domain_type::discrete_element_type;

        discrete_deriv_element_type default_derivatives
                = detail::no_derivative_element<deriv_tags>();

        remaining_deriv_element_type deriv_elements(default_derivatives);
        discrete_deriv_element_type deriv_index(
                detail::select_default(provided_deriv_domain.front(), default_derivatives));

        int const array_idx = get_array_index(deriv_index);

        internal_mdspan_type internal_view = internal_chunks[array_idx];
        auto subview_all_dims = std::experimental::submdspan(
                internal_view,
                get_slicer_for<DDims>(provided_deriv_domain, array_idx)...);
        auto subview = std::experimental::
                submdspan(subview_all_dims, get_slicer_for<DDims>(deriv_elements, array_idx)...);
        return subview;
    }

    /**
     * @brief Get the mdspan holding the values of the function from the internal array internal_chunks.
     *
     * @returns ChunkSpan The chunk span on the physical domain and the requested domains.
     */
    auto get_mdspan()
    {
        ddc::DiscreteDomain<> no_specified_dims;
        return get_mdspan(no_specified_dims);
    }

    /**
     * @brief Get the ddc::ChunkSpan which holds the values of the function.
     *
     * This function is equivalent to calling operator[] with a Domain of rank 0.
     *
     * @returns ChunkSpan The chunk span on the physical domain.
     */
    auto get_values_span()
    {
        ddc::DiscreteDomain<> no_specified_dims;
        return get_internal_chunk(no_specified_dims);
    }

    /**
     * @brief Get the ddc::ChunkSpan which holds the values of the function.
     *
     * This function is equivalent to calling operator[] with a Domain of rank 0.
     *
     * @returns ChunkSpan The chunk span on the physical domain.
     */
    auto get_values_span() const
    {
        ddc::DiscreteDomain<> no_specified_dims;
        return get_internal_chunk(no_specified_dims).span_cview();
    }
};
