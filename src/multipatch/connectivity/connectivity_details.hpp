#pragma once

#include <ddc/ddc.hpp>

namespace connectivity_details {
/**
 * @brief A class which extracts the derivative dimensions from a type sequence.
 */
template <class Patch, class InterfaceTypes>
struct PatchConnection;

template <class Patch>
struct PatchConnection<Patch, ddc::detail::TypeSeq<>>
{
    using type = ddc::detail::TypeSeq<>;
};

template <class Patch, class InterfaceType>
struct PatchConnection<Patch, ddc::detail::TypeSeq<InterfaceType>>
{
    using type = std::conditional_t<
            InterfaceType::template connected_to_patch<Patch>(),
            ddc::detail::TypeSeq<InterfaceType>,
            ddc::detail::TypeSeq<>>;
};

template <class Patch, class InterfaceType1, class... RemainingInterfaceTypes>
struct PatchConnection<Patch, ddc::detail::TypeSeq<InterfaceType1, RemainingInterfaceTypes...>>
{
    using type = ddc::type_seq_merge_t<
            typename PatchConnection<Patch, ddc::detail::TypeSeq<InterfaceType1>>::type,
            typename PatchConnection<Patch, ddc::detail::TypeSeq<RemainingInterfaceTypes...>>::
                    type>;
};

template <class TypeSeq>
struct ToTuple;

template <class... I>
struct ToTuple<ddc::detail::TypeSeq<I...>>
{
    using type = std::tuple<I...>;
};

template <class TypeSeq>
using to_tuple_t = typename ToTuple<TypeSeq>::type;


template <typename T, T val1, T... vals>
struct AllSame
{
    static constexpr T val = ((val1 == vals) && ...);
};

} // end namespace connectivity_details
