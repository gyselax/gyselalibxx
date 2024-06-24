// SPDX-License-Identifier: MIT

#pragma once

#include <ddc/ddc.hpp>

namespace detail {

/**
 * @brief A structure to determine if a tag is a derivative tag.
 *
 * @tparam Tag the class which may be a derivative.
 */
template <class Tag>
struct is_deriv_dim : std::false_type
{
};

template <class Tag>
struct is_deriv_dim<ddc::Deriv<Tag>> : std::true_type
{
};

/**
 * @brief A shortcut to get a boolean indicating if a tag is a derivative tag.
 *
 * @tparam Tag the class which may be a derivative.
 */
template <class T>
inline constexpr bool is_deriv_dim_v = is_deriv_dim<T>::value;

//-----------------------------------------------------------------------------

/**
 * @brief A class which extracts the derivative dimensions from a type sequence.
 */
template <class Tag>
struct DerivSubSet;

template <>
struct DerivSubSet<ddc::detail::TypeSeq<>>
{
    using type = ddc::detail::TypeSeq<>;
};

template <class Tag>
struct DerivSubSet<ddc::detail::TypeSeq<Tag>>
{
    using type = std::
            conditional_t<is_deriv_dim_v<Tag>, ddc::detail::TypeSeq<Tag>, ddc::detail::TypeSeq<>>;
};

template <class HeadTag, class... Tags>
struct DerivSubSet<ddc::detail::TypeSeq<HeadTag, Tags...>>
{
    using type = ddc::type_seq_merge_t<
            typename DerivSubSet<ddc::detail::TypeSeq<HeadTag>>::type,
            typename DerivSubSet<ddc::detail::TypeSeq<Tags...>>::type>;
};

/// A helper function to extract the derivative dimensions from a type sequence.
template <class Seq>
using deriv_sub_set_t = typename DerivSubSet<Seq>::type;

//-----------------------------------------------------------------------------

/**
 * @brief A class which gets the physical dimension associated with the provided derivative.
 */
template <class Tag>
struct StripDeriv;

template <class Tags>
struct StripDeriv<ddc::Deriv<Tags>>
{
    using type = Tags;
};

template <class... Tags>
struct StripDeriv<ddc::detail::TypeSeq<ddc::Deriv<Tags>...>>
{
    using type = ddc::detail::TypeSeq<Tags...>;
};

/// A helper function to get the physical dimension associated with the provided derivative.
template <class Seq>
using strip_deriv_t = typename StripDeriv<Seq>::type;

//-----------------------------------------------------------------------------

/**
 * @brief Select a tagged object by combining a limited number of known values
 * with a default value.
 *
 * For example this function can be used to build DDC objects.
 *
 * Example:
 * ddc::DiscreteElement<ddc::Deriv<X>> known_derivatives(2);
 * ddc::DiscreteElement<ddc::Deriv<X>, ddc::Deriv<Y>, ddc::Deriv<Z>> default_derivatives(0, 0, 0);
 * ddc::DiscreteElement<ddc::Deriv<X>, ddc::Deriv<Y>, ddc::Deriv<Z>> complete_derivative
 *         = select_default(known_derivatives, default_derivatives);
 * // Equivalent to:
 * // ddc::DiscreteElement<ddc::Deriv<X>, ddc::Deriv<Y>, ddc::Deriv<Z>> complete_derivative(2, 0, 0);
 *
 * @param known_values The values that are known and should appear in the result.
 * @param default_values The values that should be used if no value is provided for the related tag.
 *
 * @return A tagged object containing the known values at the tags where this
 *          information is provided and the default values everywhere else.
 */
template <
        template <typename...>
        class Container,
        class HeadQueryTag,
        class... QueryTags,
        class... Tags>
KOKKOS_FUNCTION Container<HeadQueryTag, QueryTags...> select_default(
        Container<Tags...> const& known_values,
        Container<HeadQueryTag, QueryTags...> const& default_values)
{
    if constexpr (sizeof...(QueryTags) == 0) {
        if constexpr (ddc::in_tags_v<HeadQueryTag, ddc::detail::TypeSeq<Tags...>>) {
            return ddc::select<HeadQueryTag>(known_values);
        } else {
            return ddc::select<HeadQueryTag>(default_values);
        }
    } else {
        return Container<HeadQueryTag, QueryTags...>(
                select_default(known_values, ddc::select<HeadQueryTag>(default_values)),
                select_default(known_values, ddc::select<QueryTags...>(default_values)));
    }
}

//-----------------------------------------------------------------------------

/**
 * @brief A class which provides a getter to collect the default discrete element
 * describing derivatives from a type sequence of derivative dimensions.
 */
template <class DerivSeq>
class DefaultDerivElem;

template <class... Tags>
struct DefaultDerivElem<ddc::detail::TypeSeq<Tags...>>
{
    static_assert((is_deriv_dim_v<Tags> && ...));
    KOKKOS_FUNCTION constexpr static ddc::DiscreteElement<Tags...> get_element()
    {
        return ddc::DiscreteElement<Tags...>(ddc::DiscreteElement<Tags>(0)...);
    }
};

/**
 * @brief Helper function to clean the call syntax of DefaultDerivElem.
 * A getter to collect the default discrete element describing derivatives
 * from a type sequence of derivative dimensions.
 *
 * @tparam A type sequence of derivative dimensions.
 *
 * @return A DiscreteElement of 0s.
 */
template <class DerivSeq>
KOKKOS_FUNCTION auto no_derivative_element()
{
    return DefaultDerivElem<DerivSeq>::get_element();
}

//-----------------------------------------------------------------------------

/**
 * @brief A helper function to get the domain which only contains the one specified element.
 *
 * @param idx The element that the domain should describe.
 *
 * @return The domain containing the point.
 */
template <class... Tag>
KOKKOS_FUNCTION ddc::DiscreteDomain<Tag...> get_domain_from_element(
        ddc::DiscreteElement<Tag...> idx)
{
    return ddc::DiscreteDomain<Tag...>(
            ddc::DiscreteDomain<Tag> {ddc::select<Tag>(idx), ddc::DiscreteVector<Tag>(1)}...);
}

} // namespace detail
