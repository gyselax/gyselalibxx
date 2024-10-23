// SPDX-License-Identifier: MIT

#pragma once
#include <array>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp> // Needed for ddc::Deriv

#include <sll/math_tools.hpp>

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

/**
 * @brief An abstract class which holds a chunk of memory describing a field and its derivatives.
 * This is the superclass for DerivFieldMem and DerivField.
 *
 * @tparam FieldType The type of the object stored in the internal_fields array.
 * @tparam IdxRange<DDims...> The index range on which the internal fields are defined.
 *          This index range is the physical index range on which the values are defined combined with
 *          the index range of the derivatives of interest (e.g. IdxRange<Deriv<IDimX>, IDimX, IDimY>).
 */
template <class FieldType, class... DDims>
class DerivFieldCommon<FieldType, IdxRange<DDims...>>
{
public:
    /**
     * @brief The type of the field stored in the array
     *
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC FieldMem types are referred to as Chunk types and Field types are
     * referred to as ChunkSpan/ChunkView.
     */
    using chunk_type = FieldType;

    /// @brief A type sequence containing all derivatives present in this object.
    using deriv_tags = detail::deriv_sub_set_t<ddc::detail::TypeSeq<DDims...>>;

    /// @brief A type sequence containing all physical dimensions for which derivatives are present in this object.
    using physical_deriv_grids = typename detail::strip_deriv_t<deriv_tags>;

    /// @brief A type sequence containing all the physical dimensions on which the fields are defined.
    using physical_grids = ddc::type_seq_remove_t<ddc::detail::TypeSeq<DDims...>, deriv_tags>;

    /// @brief The type of the elements in the fields.
    using element_type = typename chunk_type::element_type;

    /**
     * @brief The IdxRange on which the fields in this object are defined.
     *
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC IdxRange types are referred to as DiscreteDomain types.
     */
    using discrete_domain_type = IdxRange<DDims...>;
    /// @brief The IdxRange on which the fields in this object are defined.
    using index_range_type = discrete_domain_type;

    /**
     * @brief The Idx which can be used to index this object.
     *
     * This is a DDC keyword used to make this class interchangeable with Field.
     * In DDC Idx types are referred to as DiscreteElement types.
     */
    using discrete_element_type = Idx<DDims...>;
    /// @brief The Idx which can be used to index this object.
    using index_type = discrete_element_type;

protected:
    /// @brief The type of the memory block stored in the array internal_fields
    using internal_mdspan_type = std::experimental::mdspan<
            element_type,
            std::experimental::dextents<std::size_t, sizeof...(DDims)>,
            std::experimental::layout_stride>;

    /// @brief The type of a constant view on the memory block stored in the array internal_fields
    using internal_mdview_type = std::experimental::mdspan<
            const element_type,
            std::experimental::dextents<std::size_t, sizeof...(DDims)>,
            std::experimental::layout_stride>;

    /// @brief The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field.
    using chunk_span = typename FieldType::span_type;

    /// @brief The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field.
    using chunk_view = typename FieldType::view_type;

    /// @brief The index range for the field excluding derivatives
    using physical_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<physical_grids>;

    /// @brief The Idx which describes the physical position where values are defined.
    using physical_index_type = typename physical_idx_range_type::discrete_element_type;

    /// @brief The IdxRange which describes the derivatives present on each field.
    using discrete_deriv_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<deriv_tags>;

    /** @brief The Idx which describes the order of the derivatives in each dimension.
     * (e.g. second-order derivative).
     */
    using discrete_deriv_index_type = typename discrete_deriv_idx_range_type::discrete_element_type;

    /** @brief The Idx which describes the order of the derivatives in each dimension.
     * (e.g. second-order derivative).
     */
    using discrete_deriv_vector_type = typename discrete_deriv_idx_range_type::mlength_type;

    /// @brief The number of fields which must be created to describe this object.
    static constexpr int n_fields = 1 << ddc::type_seq_size_v<deriv_tags>;

    template <class, class, int, class>
    friend class DerivFieldMem;

    template <class, class, class, class>
    friend class DerivField;

protected:
    /** @brief The internal fields describing the values and derivatives.
     *
     * The fields which contain the values have different index ranges to the fields containing derivatives
     * so a DDC object cannot be used directly.
     * E.g. for a 2D field (X,Y) with derivatives provided in both directions the elements of internal_fields
     * have the type : DFieldMem<IdxRange<Deriv<IDimX>, Deriv<IDimY>, IDimX, IDimY>
     * The derivative index ranges are then defined such that the  elements of internal_fields represent:
     * 0 : @f$f(x,y)@f$
     * 1 : @f$\partial_x^k f(x,y)  \quad \forall 1 \leq k \leq NDerivs@f$
     * 2 : @f$\partial_y^k f(x,y)  \quad \forall 1 \leq k \leq NDerivs@f$
     * 3 : @f$\partial_x^j \partial_y^k f(x,y)  \quad \forall 1 \leq j \leq NDerivs,  \forall 1 \leq k \leq NDerivs@f$
     */
    std::array<internal_mdspan_type, n_fields> internal_fields;

    /// @brief The physical index range on which the values are defined.
    physical_idx_range_type m_physical_idx_range;

    /// @brief The index range of available derivatives.
    discrete_deriv_idx_range_type m_deriv_idx_range;

    /// @brief The physical index ranges on which the derivatives are defined.
    to_subidx_range_collection<physical_deriv_grids> m_cross_derivative_idx_range;

protected:
    /**
     * @brief An internal function which provides the index of an element inside the internal_fields.
     * An Idx describes the element of interest. If information about the derivatives is
     * missing then it is assumed that the 0-th order derivative is requested.
     *
     * @param elem The element of interest.
     *
     * @returns int The index of the internal field inside the array internal_fields.
     * @returns index_type The index of the element of interest inside the field of interest.
     */
    template <class DElem>
    KOKKOS_FUNCTION std::pair<int, index_type> get_index(DElem elem) const
    {
        discrete_deriv_index_type default_derivatives = detail::no_derivative_element<deriv_tags>();
        discrete_deriv_index_type deriv_index(detail::select_default(elem, default_derivatives));
        physical_index_type physical_index(elem);
        index_type index(physical_index, deriv_index);
        return std::pair<int, index_type>(get_array_index(deriv_index), index);
    }

    /**
     * @brief An internal function which provides the index of a field inside the internal_fields array.
     * An Idx describes the derivatives of interest. n-th order derivatives are stored in the
     * same field for all n!=0 so it is sufficient to provide any valid element from the derivatives.
     *
     * The index is calculated using a bit mask.
     * Mathematically the equation which determines the index in `internal_fields` is:
     * @f$ \sum_{\xi} 2^{i_\xi} \delta_{\text{using derivs in }\xi} @f$
     * where @f$i_\xi@f$ is the index of the dimension @f$\xi@f$ in the tags.
     *
     * @param idx The derivatives of interest.
     *
     * @returns int The index of the internal field inside the array internal_fields.
     * @returns discrete_deriv_idx_range_type The index range of the derivatives at the field at the index.
     */
    template <class... Tag>
    KOKKOS_FUNCTION int get_array_index(Idx<Tag...> idx) const
    {
        static_assert(std::is_same_v<Idx<Tag...>, discrete_deriv_index_type>);
        return (0 + ...
                + (int(ddc::select<Tag>(idx) != Idx<Tag>(0))
                   << ddc::type_seq_rank_v<Tag, deriv_tags>));
    }

    /**
     * @brief Get an object which can be used to slice an mdspan.
     *
     * @param slice_idx The DDC element which should be used to slice the field.
     * @param array_idx The index of the mdspan in internal_fields that will be sliced.
     *
     * @tparam QueryDDim The dimension along which we want to slice.
     *
     * @returns An index or a slice which can be used to slice an mdspan.
     */
    template <class QueryDDim, class... ODDims>
    KOKKOS_FUNCTION constexpr auto get_slicer_for(Idx<ODDims...> const& slice_idx, int array_idx)
            const
    {
        if constexpr (!ddc::in_tags_v<QueryDDim, ddc::detail::TypeSeq<ODDims...>>) {
            return std::experimental::full_extent;
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

    /**
     * @brief Get an object which can be used to slice an mdspan.
     *
     * @param slice_idx_range The DDC index range which should be used to slice the field.
     * @param array_idx The index of the mdspan in internal_fields that will be sliced.
     *
     * @tparam QueryDDim The dimension along which the we want to slice.
     *
     * @returns A slice (often in the form of a (start, end) pair) which can be used to slice an mdspan.
     */
    template <class QueryDDim, class... ODDims>
    KOKKOS_FUNCTION constexpr auto get_slicer_for(
            IdxRange<ODDims...> const& slice_idx_range,
            int array_idx) const
    {
        if constexpr (!ddc::in_tags_v<QueryDDim, ddc::detail::TypeSeq<ODDims...>>) {
            return std::experimental::full_extent;
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

    /**
     * Get a Field from a subset of one of the mdspans in internal_fields. The provided
     * IdxRange is used to slice the mdspan in such a way that the resulting mdspan can be
     * saved in a Field. This means that all information about derivatives must be provided.
     * If information about a derivative is missing then it is assumed that the 0-th order derivative
     * is requested.
     *
     * @param idx_range The index range used to slice the mdspan.
     *
     * @returns Field The subset of the internal mdspan.
     */
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
        auto subview = std::experimental::
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

    /**
     * Get a Field from a subset of one of the mdspans in internal_fields. The provided
     * Idx is used to slice the mdspan in such a way that the resulting mdspan can be
     * saved in a Field. This means that all information about derivatives must be provided.
     * If information about a derivative is missing then it is assumed that the 0-th order derivative
     * is requested.
     *
     * @param elem The element used to slice the mdspan.
     *
     * @returns Field The subset of the internal mdspan.
     */
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
        auto subview = std::experimental::
                submdspan(internal_view, get_slicer_for<DDims>(slice_idx, array_idx)...);
        // Create a Field with the expected index range
        Field<element_type,
              final_idx_range_type,
              typename chunk_type::memory_space,
              typename decltype(subview)::layout_type>
                local_field(subview, final_idx_range);

        return local_field;
    }

    /**
     * @brief Protected constructor to be used by subclasses to initialise index ranges.
     *
     * @param physical_idx_range The index range on which the values of the function are defined.
     * @param deriv_idx_range The index range of the provided derivatives.
     * @param cross_derivative_idx_range The cross product of the index ranges on which the derivatives
     *              of the function are defined.
     */
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

    /**
     * @brief Get a ConstField describing a subset of the data.
     *
     * @param slice_spec A discrete element describing the position at which these dimensions should be
     *          indexed. If information about the derivatives is missing then it is assumed that the
     *          0-th order derivative is requested.
     *
     * @returns ConstField A subset of the data.
     */
    template <class... QueryDDims>
    constexpr auto operator[](Idx<QueryDDims...> const& slice_spec) const
    {
        return get_internal_field(slice_spec).span_cview();
    }

    /**
     * @brief Get a Field describing a subset of the data.
     *
     * @param slice_spec A discrete element describing the position at which these dimensions should be
     *          indexed. If information about the derivatives is missing then it is assumed that the
     *          0-th order derivative is requested.
     *
     * @returns Field A subset of the data.
     */
    template <class... QueryDDims>
    constexpr auto operator[](Idx<QueryDDims...> const& slice_spec)
    {
        return get_internal_field(slice_spec);
    }

    /**
     * @brief Get a Field describing a subset of the data.
     * This function allows a slice to be obtained however it is designed to return a Field. It is
     * therefore not possible to request data from multiple fields (e.g. derivatives from 0 to 3).
     *
     * @param oidx_range A discrete index range describing the position at which these dimensions should be
     *          indexed. If information about the derivatives is missing then it is assumed that the
     *          0-th order derivative is requested.
     *
     * @returns Field A subset of the data.
     */
    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr auto operator[](IdxRange<QueryDDims...> const& oidx_range)
    {
        return get_internal_field(oidx_range);
    }

    /**
     * @brief Get a ConstField describing a subset of the data.
     * This function allows a slice to be obtained however it is designed to return a ConstField. It is
     * therefore not possible to request data from multiple fields (e.g. derivatives from 0 to 3).
     *
     * @param oidx_range A discrete index range describing the position at which these dimensions should be
     *          indexed. If information about the derivatives is missing then it is assumed that the
     *          0-th order derivative is requested.
     *
     * @returns ConstField A subset of the data.
     */
    template <class... QueryDDims>
    KOKKOS_FUNCTION constexpr auto operator[](IdxRange<QueryDDims...> const& oidx_range) const
    {
        return get_internal_field(oidx_range).span_cview();
    }

    /**
     * @brief Get one of the mdspans from the internal array internal_fields.
     * This function takes index ranges on the derivative directions. Where derivatives are
     * missing it is assumed that the 0-th order derivative is requested. This dimension
     * is stripped from the resulting field. This is the recommended way to access the
     * internal fields.
     *
     * @param provided_deriv_idx_range The derivative index range which should be retained.
     *
     * @returns Field The field on the physical index range and the requested index ranges.
     */
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
        auto subview_all_dims = std::experimental::submdspan(
                internal_view,
                get_slicer_for<DDims>(provided_deriv_idx_range, array_idx)...);
        auto subview = std::experimental::
                submdspan(subview_all_dims, get_slicer_for<DDims>(deriv_elements, array_idx)...);
        return subview;
    }

    /**
     * @brief Get the mdspan holding the values of the function from the internal array internal_fields.
     *
     * @returns Field The field on the physical index range and the requested index ranges.
     */
    auto get_mdspan()
    {
        IdxRange<> no_specified_dims;
        return get_mdspan(no_specified_dims);
    }

    /**
     * @brief Get the Field which holds the values of the function.
     *
     * This function is equivalent to calling operator[] with a 0D IdxRange.
     *
     * @returns Field The field on the physical index range.
     */
    auto get_values_field()
    {
        IdxRange<> no_specified_dims;
        return get_internal_field(no_specified_dims);
    }

    /**
     * @brief Get the Field which holds the values of the function.
     *
     * This function is equivalent to calling operator[] with a 0D IdxRange.
     *
     * @returns Field The field on the physical index range.
     */
    auto get_values_field() const
    {
        IdxRange<> no_specified_dims;
        return get_internal_field(no_specified_dims).span_cview();
    }
};
