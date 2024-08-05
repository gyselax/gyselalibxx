#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "interface.hpp"
#include "patch.hpp"

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

template <class QueryGrid1D, class OGrid, class BSpl1, class BSpl2, class... RemainingPatchTypes>
struct FindPatch<
        QueryGrid1D,
        ddc::detail::TypeSeq<Patch<QueryGrid1D, OGrid, BSpl1, BSpl2>, RemainingPatchTypes...>>
{
    using type = Patch<QueryGrid1D, OGrid, BSpl1, BSpl2>;
};

template <class QueryGrid1D, class OGrid, class BSpl1, class BSpl2, class... RemainingPatchTypes>
struct FindPatch<
        QueryGrid1D,
        ddc::detail::TypeSeq<Patch<OGrid, QueryGrid1D, BSpl1, BSpl2>, RemainingPatchTypes...>>
{
    using type = Patch<OGrid, QueryGrid1D, BSpl1, BSpl2>;
};

template <class Grid1D, class Patch1, class... RemainingPatchTypes>
struct FindPatch<Grid1D, ddc::detail::TypeSeq<Patch1, RemainingPatchTypes...>>
{
    static_assert(!std::is_same_v<typename Patch1::Grid1, Grid1D>);
    static_assert(!std::is_same_v<typename Patch1::Grid2, Grid1D>);
    using type = typename FindPatch<Grid1D, ddc::detail::TypeSeq<RemainingPatchTypes...>>::type;
};

template <class Edge, class InterfaceTypeSeq>
struct FindInterface;

template <class Edge>
struct FindInterface<Edge, ddc::detail::TypeSeq<>>
{
    static_assert(!std::is_same_v<Edge, Edge>, "Edge not found in collection of interfaces");
};

template <class Edge, class Interface1, class... RemainingInterfaceTypes>
struct FindInterface<Edge, ddc::detail::TypeSeq<Interface1, RemainingInterfaceTypes...>>
{
    static_assert(!std::is_same_v<Edge, typename Interface1::Edge1>);
    static_assert(!std::is_same_v<Edge, typename Interface1::Edge2>);
    using type =
            typename FindInterface<Edge, ddc::detail::TypeSeq<RemainingInterfaceTypes...>>::type;
};

template <class Edge, class OEdge, bool Orientations, class... RemainingInterfaceTypes>
struct FindInterface<
        Edge,
        ddc::detail::TypeSeq<Interface<Edge, OEdge, Orientations>, RemainingInterfaceTypes...>>
{
    using type = Interface<Edge, OEdge, Orientations>;
};

template <class Edge, class OEdge, bool Orientations, class... RemainingInterfaceTypes>
struct FindInterface<
        Edge,
        ddc::detail::TypeSeq<Interface<OEdge, Edge, Orientations>, RemainingInterfaceTypes...>>
{
    using type = Interface<OEdge, Edge, Orientations>;
};

template <class EdgeType>
struct SwapEnd;

template <class Patch, class Grid1D>
struct SwapEnd<Edge<Patch, Grid1D, FRONT>>
{
    using type = Edge<Patch, Grid1D, BACK>;
};

template <class Patch, class Grid1D>
struct SwapEnd<Edge<Patch, Grid1D, BACK>>
{
    using type = Edge<Patch, Grid1D, FRONT>;
};

template <class Edge, class InterfaceTypeSeq>
struct GetNextEdge
{
    using current_interface = typename FindInterface<Edge, InterfaceTypeSeq>::type;
    using type = typename SwapEnd<typename current_interface::template OtherEdge<Edge>>::type;
};

template <class StartEdge, class InterfaceTypeSeq, class FoundGrids = ddc::detail::TypeSeq<>>
struct CollectGridsOnDimFrontToBack
{
    using current_front_to_back_interface =
            typename FindInterface<StartEdge, InterfaceTypeSeq>::type;
    using MatchingBackEdge =
            typename current_front_to_back_interface::template OtherEdge<StartEdge>;
    using type = std::conditional_t<
            std::is_same_v<
                    MatchingBackEdge,
                    OutsideEdge> || ddc::in_tags_v<StartEdge::associated_patch, FoundGrids>,
            FoundGrids,
            CollectGridsOnDimFrontToBack<
                    typename SwapEnd<MatchingBackEdge>::type,
                    InterfaceTypeSeq,
                    ddc::type_seq_merge_t<
                            FoundGrids,
                            ddc::detail::TypeSeq<typename StartEdge::grid>>>>;
};

template <class StartEdge, class InterfaceTypeSeq, class FoundGrids = ddc::detail::TypeSeq<>>
struct CollectGridsOnDimBackToFront
{
    using current_back_to_front_interface =
            typename FindInterface<StartEdge, InterfaceTypeSeq>::type;
    using MatchingFrontEdge =
            typename current_back_to_front_interface::template OtherEdge<StartEdge>;
    using type = std::conditional_t<
            std::is_same_v<
                    MatchingFrontEdge,
                    OutsideEdge> || ddc::in_tags_v<typename StartEdge::associated_patch, FoundGrids>,
            FoundGrids,
            CollectGridsOnDimBackToFront<
                    typename SwapEnd<MatchingFrontEdge>::type,
                    InterfaceTypeSeq,
                    ddc::type_seq_merge_t<
                            ddc::detail::TypeSeq<typename StartEdge::grid>,
                            FoundGrids>>>;
};

template <class StartPatch, class Grid1D, class InterfaceTypeSeq>
struct CollectGridsOnDim
{
    using type = ddc::type_seq_merge_t<
            typename CollectGridsOnDimBackToFront<
                    Edge<StartPatch, Grid1D, BACK>,
                    InterfaceTypeSeq>::type,
            typename CollectGridsOnDimFrontToBack<
                    Edge<StartPatch, Grid1D, FRONT>,
                    InterfaceTypeSeq>::type>;
};

template <class QueryGrid1D, class IdxRangeType>
struct SelectRelevantIdxRangeType;

template <class QueryGrid1D, class... IdxRangeGrids>
struct SelectRelevantIdxRangeType<QueryGrid1D, IdxRange<IdxRangeGrids...>>
{
    using type = std::conditional_t<
            ddc::in_tags_v<QueryGrid1D, ddc::detail::TypeSeq<IdxRangeGrids...>>,
            ddc::detail::TypeSeq<IdxRange<IdxRangeGrids...>>,
            ddc::detail::TypeSeq<>>;
};

template <class QueryGrid1D, class IdxRangeTuple>
struct FindRelevantIdxRangeType;

template <class QueryGrid1D>
struct FindRelevantIdxRangeType<QueryGrid1D, std::tuple<>>
{
    using type = ddc::detail::TypeSeq<>;
};

template <class QueryGrid1D, class HeadIdxRangeType, class... IdxRangeTypes>
struct FindRelevantIdxRangeType<QueryGrid1D, std::tuple<HeadIdxRangeType, IdxRangeTypes...>>
{
    using type = ddc::type_seq_merge_t<
            typename SelectRelevantIdxRangeType<QueryGrid1D, HeadIdxRangeType>::type,
            typename FindRelevantIdxRangeType<QueryGrid1D, std::tuple<IdxRangeTypes...>>::type>;
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
using find_patch_t = typename connectivity_details::FindPatch<Grid1D, PatchTypeSeq>::type;

template <class EdgeType, class InterfaceTypeSeq>
using find_associated_interface_t =
        typename connectivity_details::FindInterface<EdgeType, InterfaceTypeSeq>::type;

template <class StartPatch, class Grid1D, class InterfaceTypeSeq>
using collect_grids_on_dim_t = typename connectivity_details::
        CollectGridsOnDim<StartPatch, Grid1D, InterfaceTypeSeq>::type;

template <class QueryGrid1D, class IdxRangeTuple>
using find_relevant_idx_range_t = ddc::type_seq_element_t<
        0,
        typename connectivity_details::FindRelevantIdxRangeType<QueryGrid1D, IdxRangeTuple>::type>;
