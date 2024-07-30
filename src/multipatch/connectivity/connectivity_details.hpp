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
struct FilterEdges;

template <>
struct FilterEdges<ddc::detail::TypeSeq<>>
{
};

template <class EdgeType>
struct FilterEdges<ddc::detail::TypeSeq<EdgeType>>
{
    using type = std::conditional_t<
            std::is_same_v<EdgeType, OutsideEdge>,
            ddc::detail::TypeSeq<>,
            ddc::detail::TypeSeq<EdgeType>>;
};

template <class EdgeType1, class... RemainingEdgeTypes>
struct FilterEdges<ddc::detail::TypeSeq<EdgeType1, RemainingEdgeTypes...>>
{
    using type = ddc::type_seq_merge_t<
            typename FilterEdges<ddc::detail::TypeSeq<EdgeType1>>::type,
            typename FilterEdges<ddc::detail::TypeSeq<RemainingEdgeTypes...>>::type>;
};

template <class TypeSeq>
struct ToTuple;

template <class... I>
struct ToTuple<ddc::detail::TypeSeq<I...>>
{
    using type = std::tuple<I...>;
};

template <class TypeSeq>
struct ExtractPatches;

template <>
struct ExtractPatches<ddc::detail::TypeSeq<>>
{
    using type = ddc::detail::TypeSeq<>;
};

template <class EdgeType1, class... EdgeTypes>
struct ExtractPatches<ddc::detail::TypeSeq<EdgeType1, EdgeTypes...>>
{
    using type = ddc::type_seq_merge_t<
            ddc::detail::TypeSeq<typename EdgeType1::associated_patch>,
            typename ExtractPatches<ddc::detail::TypeSeq<EdgeTypes...>>::type>;
};

template <class Grid1D, class PatchTypeSeq>
struct FindPatch;

template <class Grid1D>
struct FindPatch<Grid1D, ddc::detail::TypeSeq<>>
{
    static_assert(!std::is_same_v<Grid1D, Grid1D>, "No patches found using dimension.");
};

template <class Grid1D, class Patch1, class... RemainingPatchTypes>
struct FindPatch<Grid1D, ddc::detail::TypeSeq<Patch1, RemainingPatchTypes...>>
{
    using type = std::conditional_t<
            std::is_same_v<Patch1::Grid1, Grid1D> || std::is_same_v<Patch1::Grid2, Grid1D>,
            Patch1,
            typename FindPatch<Grid1D, ddc::detail::TypeSeq<RemainingPatchTypes...>>::type>;
};


} // end namespace connectivity_details



template <class TypeSeq>
using to_tuple_t = typename connectivity_details::ToTuple<TypeSeq>::type;

template <class... T>
using filter_edges_t = typename connectivity_details::FilterEdges<ddc::detail::TypeSeq<T...>>::type;

template <class TypeSeq>
using extract_patches_t = typename connectivity_details::ExtractPatches<TypeSeq>::type;

template <class Patch, class InterfaceTypeSeq>
using direct_patch_connections_t =
        typename connectivity_details::PatchConnection<Patch, InterfaceTypeSeq>::type;

template <class Grid1D, class PatchTypeSeq>
using find_patch_t = typename FindPatch<Grid1D, PatchTypeSeq>::type;
