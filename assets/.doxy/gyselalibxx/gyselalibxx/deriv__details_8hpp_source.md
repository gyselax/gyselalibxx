

# File deriv\_details.hpp

[**File List**](files.md) **>** [**data\_types**](dir_eaa769653453aaefd8cc10e98e9bb3eb.md) **>** [**deriv\_details.hpp**](deriv__details_8hpp.md)

[Go to the documentation of this file](deriv__details_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

namespace detail {

template <class Tag>
struct is_deriv_dim : std::false_type
{
};

template <class Tag>
struct is_deriv_dim<ddc::Deriv<Tag>> : std::true_type
{
};

template <class T>
inline constexpr bool is_deriv_dim_v = is_deriv_dim<T>::value;

//-----------------------------------------------------------------------------

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

template <class Seq>
using deriv_sub_set_t = typename DerivSubSet<Seq>::type;

//-----------------------------------------------------------------------------

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

template <class Seq>
using strip_deriv_t = typename StripDeriv<Seq>::type;

//-----------------------------------------------------------------------------

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

template <class DerivSeq>
struct DefaultDerivElem;

template <class... Tags>
struct DefaultDerivElem<ddc::detail::TypeSeq<Tags...>>
{
    static_assert((is_deriv_dim_v<Tags> && ...));
    KOKKOS_FUNCTION constexpr static Idx<Tags...> get_element()
    {
        return Idx<Tags...>(Idx<Tags> {0}...);
    }
};

template <class DerivSeq>
KOKKOS_FUNCTION auto no_derivative_element()
{
    return DefaultDerivElem<DerivSeq>::get_element();
}

//-----------------------------------------------------------------------------

template <class... Tag>
KOKKOS_FUNCTION IdxRange<Tag...> get_idx_range_from_element(Idx<Tag...> idx)
{
    return IdxRange<Tag...>(IdxRange<Tag> {ddc::select<Tag>(idx), IdxStep<Tag> {1}}...);
}

//-----------------------------------------------------------------------------

template <class... Containers>
struct Combine;

template <class Container>
struct Combine<Container>
{
    using type = Container;
};

template <
        template <typename...>
        class Container,
        class... Tags,
        class... OTags,
        class... TailContainers>
struct Combine<Container<Tags...>, Container<OTags...>, TailContainers...>
{
    using type = Combine<Container<Tags..., OTags...>, TailContainers...>;
};

template <class... Containers>
using combine_t = typename detail::Combine<Containers...>::type;

} // namespace detail
```


