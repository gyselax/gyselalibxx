

# File derivative\_field\_common.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**derivative\_field\_common.hpp**](derivative__field__common_8hpp.md)

[Go to the documentation of this file](derivative__field__common_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <array>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp> // Needed for ddc::Deriv

#include "ddc_aliases.hpp"
#include "deriv_details.hpp"
#include "idx_range_slice.hpp"

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

template <class FieldType, class SupportType>
class DerivFieldCommon;

template <class FieldType, class... DDims>
class DerivFieldCommon<FieldType, IdxRange<DDims...>>
{
public:
    using chunk_type = FieldType;

    using deriv_tags = detail::deriv_sub_set_t<ddc::detail::TypeSeq<DDims...>>;

    using physical_deriv_grids = typename detail::strip_deriv_t<deriv_tags>;

    using physical_grids = ddc::type_seq_remove_t<ddc::detail::TypeSeq<DDims...>, deriv_tags>;

    using element_type = typename chunk_type::element_type;

    using discrete_domain_type = IdxRange<DDims...>;
    using index_range_type = discrete_domain_type;

    using discrete_element_type = Idx<DDims...>;
    using index_type = discrete_element_type;

protected:
    using internal_mdspan_type = Kokkos::mdspan<
            element_type,
            Kokkos::dextents<std::size_t, sizeof...(DDims)>,
            Kokkos::layout_stride>;

    using internal_mdview_type = Kokkos::mdspan<
            const element_type,
            Kokkos::dextents<std::size_t, sizeof...(DDims)>,
            Kokkos::layout_stride>;

    using chunk_span = typename FieldType::span_type;

    using chunk_view = typename FieldType::view_type;

    using physical_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<physical_grids>;

    using physical_index_type = typename physical_idx_range_type::discrete_element_type;

    using discrete_deriv_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<deriv_tags>;

    using discrete_deriv_index_type = typename discrete_deriv_idx_range_type::discrete_element_type;

    using discrete_deriv_vector_type = typename discrete_deriv_idx_range_type::discrete_vector_type;

    static constexpr int n_fields = 1 << ddc::type_seq_size_v<deriv_tags>;

    template <class, class, int, class>
    friend class DerivFieldMem;

    template <class, class, class, class>
    friend class DerivField;

protected:
    std::array<internal_mdspan_type, n_fields> internal_fields;

    physical_idx_range_type m_physical_idx_range;

    discrete_deriv_idx_range_type m_deriv_idx_range;

    to_subidx_range_collection<physical_deriv_grids> m_cross_derivative_idx_range;

protected:
    template <class DElem>
    KOKKOS_FUNCTION std::pair<int, index_type> get_index(DElem elem) const
    {
        discrete_deriv_index_type default_derivatives = detail::no_derivative_element<deriv_tags>();
        discrete_deriv_index_type deriv_index(detail::select_default(elem, default_derivatives));
        physical_index_type physical_index(elem);
        index_type index(physical_index, deriv_index);
        return std::pair<int, index_type>(get_array_index(deriv_index), index);
    }

    template <class... Tag>
    KOKKOS_FUNCTION int get_array_index(Idx<Tag...> idx) const
    {
        static_assert(std::is_same_v<Idx<Tag...>, discrete_deriv_index_type>);
        return (0 + ...
                + (int(ddc::select<Tag>(idx) != Idx<Tag>(0))
                   << ddc::type_seq_rank_v<Tag, deriv_tags>));
    }

    template <class QueryDDim, class... ODDims>
    KOKKOS_FUNCTION constexpr auto get_slicer_for(Idx<ODDims...> const& slice_idx, int array_idx)
            const
    {
        if constexpr (!ddc::in_tags_v<QueryDDim, ddc::detail::TypeSeq<ODDims...>>) {
            return Kokkos::full_extent;
        } else {
            if constexpr (ddc::in_tags_v<QueryDDim, physical_deriv_grids>) {
                // Physical dimension along which derivatives are known
                // If information is available about the physical index range
                if (array_idx & (1 << ddc::type_seq_rank_v<ddc::Deriv<QueryDDim>, deriv_tags>)) {
                    // If the derivative is being requested
                    return m_cross_derivative_idx_range.get_index(
                            ddc::select<QueryDDim>(slice_idx));
                }
            }
            if constexpr (ddc::in_tags_v<QueryDDim, physical_grids>) {
                // Physical dimension along which derivatives are not known
                return std::size_t(
                        (ddc::select<QueryDDim>(slice_idx)
                         - ddc::select<QueryDDim>(m_physical_idx_range).front()));
            }
            if constexpr (ddc::in_tags_v<QueryDDim, deriv_tags>) {
                // Derivative dimension
                if (array_idx & (1 << ddc::type_seq_rank_v<QueryDDim, deriv_tags>)) {
                    // If array contains derivatives
                    return std::size_t((ddc::select<QueryDDim>(slice_idx) - Idx<QueryDDim>(1)));
                } else {
                    // If array doesn't contain derivatives
                    assert(ddc::select<QueryDDim>(slice_idx) == Idx<QueryDDim>(0));
                    return std::size_t(0);
                }
            }
        }
    }

    template <class QueryDDim, class... ODDims>
    KOKKOS_FUNCTION constexpr auto get_slicer_for(
            IdxRange<ODDims...> const& slice_idx_range,
            int array_idx) const
    {
        if constexpr (!ddc::in_tags_v<QueryDDim, ddc::detail::TypeSeq<ODDims...>>) {
            return Kokkos::full_extent;
        } else {
            if constexpr (ddc::in_tags_v<QueryDDim, physical_deriv_grids>) {
                // Physical dimension along which derivatives are known
                // If information is available about the physical index range
                IdxRange<QueryDDim> idx_range_requested(slice_idx_range);
                if (array_idx & (1 << ddc::type_seq_rank_v<ddc::Deriv<QueryDDim>, deriv_tags>)) {
                    // If the derivative is being requested
                    assert(m_cross_derivative_idx_range.contains(idx_range_requested));
                    return std::pair<std::size_t, std::size_t>(
                            m_cross_derivative_idx_range.get_index(idx_range_requested.front()),
                            m_cross_derivative_idx_range.get_index(idx_range_requested.back()) + 1);
                }
            }
            if constexpr (ddc::in_tags_v<QueryDDim, physical_grids>) {
                // Physical dimension along which derivatives are not known
                return std::pair<std::size_t, std::size_t>(
                        ddc::select<QueryDDim>(slice_idx_range).front()
                                - ddc::select<QueryDDim>(m_physical_idx_range).front(),
                        ddc::select<QueryDDim>(slice_idx_range).back() + 1
                                - ddc::select<QueryDDim>(m_physical_idx_range).front());
            }
            if constexpr (ddc::in_tags_v<QueryDDim, deriv_tags>) {
                // Derivative dimension
                if (array_idx & (1 << ddc::type_seq_rank_v<QueryDDim, deriv_tags>)) {
                    // If array contains derivatives
                    return std::pair<std::size_t, std::size_t>(
                            ddc::select<QueryDDim>(slice_idx_range).front() - Idx<QueryDDim>(1),
                            ddc::select<QueryDDim>(slice_idx_range).back() + 1 - Idx<QueryDDim>(1));
                } else {
                    // If array doesn't contain derivatives
                    assert(ddc::select<QueryDDim>(slice_idx_range).front() == Idx<QueryDDim>(0));
                    assert(ddc::select<QueryDDim>(slice_idx_range).back()
                           == ddc::select<QueryDDim>(slice_idx_range).front());
                    return std::pair<std::size_t, std::size_t>(0, 1);
                }
            }
        }
    }

    template <class... ODims>
    KOKKOS_FUNCTION auto get_internal_field(IdxRange<ODims...> idx_range) const
    {
        // Get the types related to the provided information
        using provided_tags = ddc::detail::TypeSeq<ODims...>;
        using provided_deriv_tags = detail::deriv_sub_set_t<provided_tags>;

        // Get the types related to the implicit information
        using remaining_deriv_tags = ddc::type_seq_remove_t<deriv_tags, provided_deriv_tags>;
        using remaining_deriv_idx_range_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<remaining_deriv_tags>;

        // Find the index range of the derivatives (either provided or an index range containing only the 0-th derivative)
        remaining_deriv_idx_range_type no_deriv_idx_range = detail::get_idx_range_from_element(
                detail::no_derivative_element<remaining_deriv_tags>());
        discrete_deriv_idx_range_type deriv_idx_range(idx_range, no_deriv_idx_range);

        // Find the physical index range of the field
        physical_idx_range_type local_physical_idx_range(
                detail::select_default(idx_range, m_physical_idx_range));

        // Find the discrete index range of the field
        index_range_type full_idx_range(local_physical_idx_range, deriv_idx_range);

        // Find the index of the internal field
        int const array_idx = get_array_index(deriv_idx_range.front());

        // Get the relevant internal field
        internal_mdspan_type internal_view = internal_fields[array_idx];
        // Slice the relevant section of the internal field
        auto subview = Kokkos::
                submdspan(internal_view, get_slicer_for<DDims>(full_idx_range, array_idx)...);
        // Create a Field with the expected index range
        Field<element_type,
              index_range_type,
              typename chunk_type::memory_space,
              typename decltype(subview)::layout_type>
                local_field(subview, full_idx_range);

        // If necessary, slice off the derivative dimensions deduced implicitly
        if constexpr (ddc::type_seq_size_v<remaining_deriv_tags> == 0) {
            return local_field;
        } else {
            return local_field[no_deriv_idx_range.front()];
        }
    }

    template <class... ODims>
    KOKKOS_FUNCTION auto get_internal_field(Idx<ODims...> elem) const
    {
        // Get the types related to the provided information
        using provided_tags = ddc::detail::TypeSeq<ODims...>;
        using provided_deriv_tags = detail::deriv_sub_set_t<provided_tags>;
        using provided_physical_tags = ddc::type_seq_remove_t<provided_tags, provided_deriv_tags>;
        using provided_deriv_idx_range_type
                = ddc::detail::convert_type_seq_to_discrete_domain_t<provided_deriv_tags>;
        using provided_deriv_index_type =
                typename provided_deriv_idx_range_type::discrete_element_type;

        // Get the types related to the implicit information
        using remaining_deriv_tags = ddc::type_seq_remove_t<deriv_tags, provided_deriv_tags>;
        using remaining_deriv_idx_range_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<remaining_deriv_tags>;
        using remaining_deriv_index_type =
                typename remaining_deriv_idx_range_type::discrete_element_type;

        // Get the types related to the final field type
        using sliced_tags = ddc::type_seq_merge_t<provided_physical_tags, deriv_tags>;
        using sliced_idx_range_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<sliced_tags>;
        using sliced_index_type = typename sliced_idx_range_type::discrete_element_type;
        using final_tags = ddc::type_seq_remove_t<ddc::detail::TypeSeq<DDims...>, sliced_tags>;
        using final_idx_range_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<final_tags>;

        // Get the index of the relevant derivatives
        provided_deriv_index_type requested_derivs(elem);
        remaining_deriv_index_type no_deriv_elements(
                detail::no_derivative_element<remaining_deriv_tags>());
        discrete_deriv_index_type deriv_index(requested_derivs, no_deriv_elements);

        // Get the element which will slice the mdspan
        sliced_index_type slice_idx(elem, no_deriv_elements);

        // Find the index of the internal field
        int const array_idx = get_array_index(deriv_index);

        // Get the final index range
        final_idx_range_type final_idx_range(m_physical_idx_range);

        // Get the relevant internal field
        internal_mdspan_type internal_view = internal_fields[array_idx];
        // Slice the relevant section of the internal field
        auto subview
                = Kokkos::submdspan(internal_view, get_slicer_for<DDims>(slice_idx, array_idx)...);
        // Create a Field with the expected index range
        Field<element_type,
              final_idx_range_type,
              typename chunk_type::memory_space,
              typename decltype(subview)::layout_type>
                local_field(subview, final_idx_range);

        return local_field;
    }

    KOKKOS_FUNCTION DerivFieldCommon(
            physical_idx_range_type physical_idx_range,
            discrete_deriv_idx_range_type deriv_idx_range,
            to_subidx_range_collection<physical_deriv_grids> cross_derivative_idx_range)
        : m_physical_idx_range(physical_idx_range)
        , m_deriv_idx_range(deriv_idx_range)
        , m_cross_derivative_idx_range(cross_derivative_idx_range)
    {
    }

public:
    KOKKOS_DEFAULTED_FUNCTION ~DerivFieldCommon() = default;

    template <class... QueryDDims>
    constexpr auto operator[](Idx<QueryDDims...> const& slice_spec) const
    {
        return get_internal_field(slice_spec).span_cview();
    }

    template <class... QueryDDims>
    constexpr auto operator[](Idx<QueryDDims...> const& slice_spec)
    {
        return get_internal_field(slice_spec);
    }

    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr auto operator[](IdxRange<QueryDDims...> const& oidx_range)
    {
        return get_internal_field(oidx_range);
    }

    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr auto operator[](IdxRange<QueryDDims...> const& oidx_range) const
    {
        return get_internal_field(oidx_range).span_cview();
    }

    template <class... ODims>
    auto get_mdspan(IdxRange<ODims...> provided_deriv_idx_range)
    {
        static_assert(((ddc::in_tags_v<ODims, deriv_tags>)&&...));
        using provided_deriv_tags = ddc::detail::TypeSeq<ODims...>;
        using remaining_deriv_tags = ddc::type_seq_remove_t<deriv_tags, provided_deriv_tags>;
        using remaining_deriv_idx_range_type =
                typename ddc::detail::convert_type_seq_to_discrete_domain_t<remaining_deriv_tags>;
        using remaining_deriv_index_type =
                typename remaining_deriv_idx_range_type::discrete_element_type;

        discrete_deriv_index_type default_derivatives = detail::no_derivative_element<deriv_tags>();

        remaining_deriv_index_type deriv_elements(default_derivatives);
        discrete_deriv_index_type deriv_index(
                detail::select_default(provided_deriv_idx_range.front(), default_derivatives));

        int const array_idx = get_array_index(deriv_index);

        internal_mdspan_type internal_view = internal_fields[array_idx];
        auto subview_all_dims = Kokkos::submdspan(
                internal_view,
                get_slicer_for<DDims>(provided_deriv_idx_range, array_idx)...);
        auto subview = Kokkos::
                submdspan(subview_all_dims, get_slicer_for<DDims>(deriv_elements, array_idx)...);
        return subview;
    }

    auto get_mdspan()
    {
        IdxRange<> no_specified_dims;
        return get_mdspan(no_specified_dims);
    }

    auto get_values_field()
    {
        IdxRange<> no_specified_dims;
        return get_internal_field(no_specified_dims);
    }

    auto get_values_field() const
    {
        IdxRange<> no_specified_dims;
        return get_internal_field(no_specified_dims).span_cview();
    }
};
```


