

# File connectivity\_details.hpp

[**File List**](files.md) **>** [**connectivity**](dir_28b51abc9241105ab41b66c468e7d019.md) **>** [**connectivity\_details.hpp**](connectivity__details_8hpp.md)

[Go to the documentation of this file](connectivity__details_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include <array>

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "interface.hpp"
#include "patch.hpp"
#include "type_seq_tools.hpp"

namespace connectivity_details {

template <concepts::Edge Edge, concepts::Edge OEdge, bool Orientations>
using Interface = std::conditional_t<
        (get_edge_patch_id<Edge>() < get_edge_patch_id<OEdge>()),
        ::Interface<Edge, OEdge, Orientations>,
        ::Interface<OEdge, Edge, Orientations>>;

template <concepts::Patch Patch, ddc::concepts::type_seq InterfaceTypeSeq>
struct PatchConnection;

template <concepts::Patch Patch>
struct PatchConnection<Patch, ddc::TypeSeq<>>
{
    using type = ddc::TypeSeq<>;
};

template <concepts::Patch Patch, class InterfaceType>
struct PatchConnection<Patch, ddc::TypeSeq<InterfaceType>>
{
    using type = std::conditional_t<
            InterfaceType::template connected_to_patch<Patch>(),
            ddc::TypeSeq<InterfaceType>,
            ddc::TypeSeq<>>;
};

template <concepts::Patch Patch, class InterfaceType1, class... RemainingInterfaceTypes>
struct PatchConnection<Patch, ddc::TypeSeq<InterfaceType1, RemainingInterfaceTypes...>>
{
    using type = ddc::type_seq_merge_t<
            typename PatchConnection<Patch, ddc::TypeSeq<InterfaceType1>>::type,
            typename PatchConnection<Patch, ddc::TypeSeq<RemainingInterfaceTypes...>>::type>;
};

template <ddc::concepts::type_seq TypeSeq>
struct StripOutsideEdges;

template <concepts::Edge EdgeType>
struct StripOutsideEdges<ddc::TypeSeq<EdgeType>>
{
    using type = std::conditional_t<
            std::is_same_v<EdgeType, OutsideEdge>,
            ddc::TypeSeq<>,
            ddc::TypeSeq<EdgeType>>;
};

template <concepts::Edge EdgeType1, class... RemainingEdgeTypes>
struct StripOutsideEdges<ddc::TypeSeq<EdgeType1, RemainingEdgeTypes...>>
{
    using type = ddc::type_seq_merge_t<
            typename StripOutsideEdges<ddc::TypeSeq<EdgeType1>>::type,
            typename StripOutsideEdges<ddc::TypeSeq<RemainingEdgeTypes...>>::type>;
};

template <ddc::concepts::type_seq TypeSeq>
struct ToTuple;

template <class... I>
struct ToTuple<ddc::TypeSeq<I...>>
{
    using type = std::tuple<I...>;
};

template <ddc::concepts::type_seq TypeSeq>
struct ExtractPatches;

template <>
struct ExtractPatches<ddc::TypeSeq<>>
{
    using type = ddc::TypeSeq<>;
};

template <concepts::Edge EdgeType1, concepts::Edge... EdgeTypes>
struct ExtractPatches<ddc::TypeSeq<EdgeType1, EdgeTypes...>>
{
    using type = ddc::type_seq_merge_t<
            ddc::TypeSeq<typename EdgeType1::associated_patch>,
            typename ExtractPatches<ddc::TypeSeq<EdgeTypes...>>::type>;
};

template <ddc::concepts::discrete_dimension Grid1D, ddc::concepts::type_seq PatchTypeSeq>
struct FindPatch;

template <ddc::concepts::discrete_dimension Grid1D>
struct FindPatch<Grid1D, ddc::TypeSeq<>>
{
    // This class should not be reached in the recursion. Raise a readable error (based on the
    // template arguments) if it is instantiated.
    static_assert(!std::is_same_v<Grid1D, Grid1D>, "No patches found using dimension.");
};

template <
        class QueryGrid1D,
        std::size_t PatchId,
        class OGrid,
        class BSpl1,
        class BSpl2,
        class... RemainingPatchTypes>
struct FindPatch<
        QueryGrid1D,
        ddc::TypeSeq<Patch<PatchId, QueryGrid1D, OGrid, BSpl1, BSpl2>, RemainingPatchTypes...>>
{
    using type = Patch<PatchId, QueryGrid1D, OGrid, BSpl1, BSpl2>;
};

template <
        class QueryGrid1D,
        std::size_t PatchId,
        class OGrid,
        class BSpl1,
        class BSpl2,
        class... RemainingPatchTypes>
struct FindPatch<
        QueryGrid1D,
        ddc::TypeSeq<Patch<PatchId, OGrid, QueryGrid1D, BSpl1, BSpl2>, RemainingPatchTypes...>>
{
    using type = Patch<PatchId, OGrid, QueryGrid1D, BSpl1, BSpl2>;
};

template <
        ddc::concepts::discrete_dimension Grid1D,
        concepts::Patch Patch1,
        concepts::Patch... RemainingPatchTypes>
struct FindPatch<Grid1D, ddc::TypeSeq<Patch1, RemainingPatchTypes...>>
{
    // Should instantiate specialisations above instead of this specialisation
    static_assert(!std::is_same_v<typename Patch1::Grid1, Grid1D>);
    static_assert(!std::is_same_v<typename Patch1::Grid2, Grid1D>);
    using type = typename FindPatch<Grid1D, ddc::TypeSeq<RemainingPatchTypes...>>::type;
};

template <concepts::Edge Edge, ddc::concepts::type_seq InterfaceTypeSeq>
struct FindInterface;

template <concepts::Edge Edge>
struct FindInterface<Edge, ddc::TypeSeq<>>
{
    // This class should not be reached in the recursion. Raise a readable error (based on the
    // template arguments) if it is instantiated.
    static_assert(!std::is_same_v<Edge, Edge>, "Edge not found in collection of interfaces.");
};

template <concepts::Edge Edge, class Interface1, class... RemainingInterfaceTypes>
struct FindInterface<Edge, ddc::TypeSeq<Interface1, RemainingInterfaceTypes...>>
{
    // Should instantiate specialisations below instead of this specialisation
    static_assert(!std::is_same_v<Edge, typename Interface1::Edge1>);
    static_assert(!std::is_same_v<Edge, typename Interface1::Edge2>);
    using type = typename FindInterface<Edge, ddc::TypeSeq<RemainingInterfaceTypes...>>::type;
};

template <
        concepts::Edge Edge,
        concepts::Edge OEdge,
        bool Orientations,
        class... RemainingInterfaceTypes>
struct FindInterface<
        Edge,
        ddc::TypeSeq<::Interface<Edge, OEdge, Orientations>, RemainingInterfaceTypes...>>
{
    using type = Interface<Edge, OEdge, Orientations>;
};

template <
        concepts::Edge Edge,
        concepts::Edge OEdge,
        bool Orientations,
        class... RemainingInterfaceTypes>
struct FindInterface<
        Edge,
        ddc::TypeSeq<::Interface<OEdge, Edge, Orientations>, RemainingInterfaceTypes...>>
{
    using type = Interface<OEdge, Edge, Orientations>;
};

template <class InterfaceType, concepts::Edge FirstEdge>
struct EnforceFirstInterfaceEdge;

template <concepts::Edge FirstEdge, concepts::Edge Edge2, bool Orientations>
struct EnforceFirstInterfaceEdge<::Interface<FirstEdge, Edge2, Orientations>, FirstEdge>
{
    using type = Interface<FirstEdge, Edge2, Orientations>;
};

template <concepts::Edge FirstEdge, concepts::Edge Edge2, bool Orientations>
struct EnforceFirstInterfaceEdge<::Interface<Edge2, FirstEdge, Orientations>, FirstEdge>
{
    using type = Interface<FirstEdge, Edge2, Orientations>;
};

template <class InterfaceType, concepts::Edge FirstEdge>
using enforce_first_interface_edge_t =
        typename EnforceFirstInterfaceEdge<InterfaceType, FirstEdge>::type;

template <concepts::Edge EdgeType>
struct SwapExtremity;

template <concepts::Patch Patch, ddc::concepts::discrete_dimension Grid1D>
struct SwapExtremity<Edge<Patch, Grid1D, FRONT>>
{
    using type = Edge<Patch, Grid1D, BACK>;
};

template <concepts::Patch Patch, ddc::concepts::discrete_dimension Grid1D>
struct SwapExtremity<Edge<Patch, Grid1D, BACK>>
{
    using type = Edge<Patch, Grid1D, FRONT>;
};

} // namespace connectivity_details

template <concepts::Edge StartEdge, ddc::concepts::type_seq InterfaceTypeSeq>
using equivalent_edge_t = typename connectivity_details::
        FindInterface<StartEdge, InterfaceTypeSeq>::type::template OtherEdge<StartEdge>;

namespace connectivity_details {
enum InsertPosition { FrontInsert, BackInsert };

template <class ToInsert, ddc::concepts::type_seq TypeSeq, InsertPosition insert_pos>
struct AddToTypeSeq;

template <class ToInsert, ddc::concepts::type_seq TypeSeq>
struct AddToTypeSeq<ToInsert, TypeSeq, FrontInsert>
{
    using type = ddc::type_seq_merge_t<TypeSeq, ddc::TypeSeq<ToInsert>>;
};

template <class ToInsert, ddc::concepts::type_seq TypeSeq>
struct AddToTypeSeq<ToInsert, TypeSeq, BackInsert>
{
    using type = ddc::type_seq_merge_t<ddc::TypeSeq<ToInsert>, TypeSeq>;
};

template <
        concepts::Edge StartEdge,
        ddc::concepts::type_seq InterfaceTypeSeq,
        InsertPosition insert_pos,
        class FoundInterfaces = ddc::TypeSeq<>,
        concepts::Edge MatchingEdge = equivalent_edge_t<StartEdge, InterfaceTypeSeq>,
        bool interface_already_found // Periodic case
        = ddc::in_tags_v<
                enforce_first_interface_edge_t<
                        typename FindInterface<StartEdge, InterfaceTypeSeq>::type,
                        std::conditional_t<StartEdge::extremity == FRONT, MatchingEdge, StartEdge>>,
                FoundInterfaces>>
struct CollectInterfacesAlongDim;

template <
        concepts::Edge StartEdge,
        ddc::concepts::type_seq InterfaceTypeSeq,
        InsertPosition insert_pos,
        class FoundInterfaces,
        concepts::Edge MatchingEdge>
struct CollectInterfacesAlongDim<
        StartEdge,
        InterfaceTypeSeq,
        insert_pos,
        FoundInterfaces,
        MatchingEdge,
        false>
{
    using NewInterfaceList = typename AddToTypeSeq<
            enforce_first_interface_edge_t<
                    typename FindInterface<StartEdge, InterfaceTypeSeq>::type,
                    std::conditional_t<StartEdge::extremity == FRONT, MatchingEdge, StartEdge>>,
            FoundInterfaces,
            insert_pos>::type;
    using type = typename CollectInterfacesAlongDim<
            typename SwapExtremity<MatchingEdge>::type,
            InterfaceTypeSeq,
            insert_pos,
            NewInterfaceList>::type;
};

template <
        concepts::Edge StartEdge,
        ddc::concepts::type_seq InterfaceTypeSeq,
        InsertPosition insert_pos,
        class FoundInterfaces,
        concepts::Edge MatchingEdge>
struct CollectInterfacesAlongDim<
        StartEdge,
        InterfaceTypeSeq,
        insert_pos,
        FoundInterfaces,
        MatchingEdge,
        true>
{
    using type = FoundInterfaces;
};

template <
        concepts::Edge StartEdge,
        ddc::concepts::type_seq InterfaceTypeSeq,
        InsertPosition insert_pos,
        class FoundInterfaces>
struct CollectInterfacesAlongDim<
        StartEdge,
        InterfaceTypeSeq,
        insert_pos,
        FoundInterfaces,
        OutsideEdge,
        false>
{
    using type = typename AddToTypeSeq<
            enforce_first_interface_edge_t<
                    typename FindInterface<StartEdge, InterfaceTypeSeq>::type,
                    std::conditional_t<StartEdge::extremity == FRONT, OutsideEdge, StartEdge>>,
            FoundInterfaces,
            insert_pos>::type;
};

template <
        concepts::Edge StartEdge,
        ddc::concepts::type_seq InterfaceTypeSeq,
        InsertPosition insert_pos,
        ddc::concepts::type_seq FoundGrids = ddc::TypeSeq<>,
        concepts::Edge MatchingEdge = equivalent_edge_t<StartEdge, InterfaceTypeSeq>,
        bool grid_already_found // Periodic case
        = ddc::in_tags_v<typename StartEdge::perpendicular_grid, FoundGrids>>
struct CollectGridsAlongDim;

template <
        concepts::Edge StartEdge,
        ddc::concepts::type_seq InterfaceTypeSeq,
        InsertPosition insert_pos,
        ddc::concepts::type_seq FoundGrids,
        concepts::Edge MatchingEdge>
struct CollectGridsAlongDim<
        StartEdge,
        InterfaceTypeSeq,
        insert_pos,
        FoundGrids,
        MatchingEdge,
        false>
{
    using NewGridList =
            typename AddToTypeSeq<typename StartEdge::perpendicular_grid, FoundGrids, insert_pos>::
                    type;
    using type = typename CollectGridsAlongDim<
            typename SwapExtremity<MatchingEdge>::type,
            InterfaceTypeSeq,
            insert_pos,
            NewGridList>::type;
};

template <
        concepts::Edge StartEdge,
        ddc::concepts::type_seq InterfaceTypeSeq,
        InsertPosition insert_pos,
        ddc::concepts::type_seq FoundGrids,
        concepts::Edge MatchingEdge>
struct CollectGridsAlongDim<StartEdge, InterfaceTypeSeq, insert_pos, FoundGrids, MatchingEdge, true>
{
    using type = FoundGrids;
};

template <
        concepts::Edge StartEdge,
        ddc::concepts::type_seq InterfaceTypeSeq,
        InsertPosition insert_pos,
        ddc::concepts::type_seq FoundGrids>
struct CollectGridsAlongDim<StartEdge, InterfaceTypeSeq, insert_pos, FoundGrids, OutsideEdge, false>
{
    using type =
            typename AddToTypeSeq<typename StartEdge::perpendicular_grid, FoundGrids, insert_pos>::
                    type;
};

template <
        concepts::Patch StartPatch,
        ddc::concepts::discrete_dimension Grid1D,
        ddc::concepts::type_seq InterfaceTypeSeq>
struct CollectAllGridsOnDim
{
private:
    using BackwardTypeSeq = typename CollectGridsAlongDim<
            Edge<StartPatch, Grid1D, FRONT>,
            InterfaceTypeSeq,
            BackInsert>::type;
    // Remove the final element so ForwardTypeSeq doesn't find these elements again in the
    // periodic case and doesn't stop upon encountering Grid1D.
    using NonOverlappingBackwardTypeSeq
            = type_seq_range_t<BackwardTypeSeq, 0, ddc::type_seq_size_v<BackwardTypeSeq> - 1>;
    // Work forward from back (end) of grid inserting each new grid at the end of the sequence
    using ForwardTypeSeq = typename CollectGridsAlongDim<
            Edge<StartPatch, Grid1D, BACK>,
            InterfaceTypeSeq,
            FrontInsert,
            NonOverlappingBackwardTypeSeq>::type;

public:
    using type = ddc::type_seq_merge_t<BackwardTypeSeq, ForwardTypeSeq>;
};

template <
        concepts::Patch StartPatch,
        ddc::concepts::discrete_dimension Grid1D,
        ddc::concepts::type_seq InterfaceTypeSeq>
struct CollectAllInterfacesOnDim
{
    using BackwardTypeSeq = typename CollectInterfacesAlongDim<
            Edge<StartPatch, Grid1D, FRONT>,
            InterfaceTypeSeq,
            BackInsert>::type;
    using type = ddc::type_seq_merge_t<
            BackwardTypeSeq,
            // Work forward from back (end) of grid inserting each new grid at the end of the sequence
            typename CollectInterfacesAlongDim<
                    Edge<StartPatch, Grid1D, BACK>,
                    InterfaceTypeSeq,
                    FrontInsert,
                    BackwardTypeSeq>::type>;
};

template <
        ddc::concepts::discrete_dimension QueryGrid1D,
        ddc::concepts::discrete_domain IdxRangeType>
struct SelectRelevantIdxRangeType;

template <
        ddc::concepts::discrete_dimension QueryGrid1D,
        ddc::concepts::discrete_dimension... IdxRangeGrids>
struct SelectRelevantIdxRangeType<QueryGrid1D, IdxRange<IdxRangeGrids...>>
{
    using type = std::conditional_t<
            ddc::in_tags_v<QueryGrid1D, ddc::TypeSeq<IdxRangeGrids...>>,
            ddc::TypeSeq<IdxRange<IdxRangeGrids...>>,
            ddc::TypeSeq<>>;
};

template <ddc::concepts::discrete_dimension QueryGrid1D, class IdxRangeTuple>
struct FindRelevantIdxRangeType;

template <ddc::concepts::discrete_dimension QueryGrid1D>
struct FindRelevantIdxRangeType<QueryGrid1D, std::tuple<>>
{
    using type = ddc::TypeSeq<>;
};

template <
        ddc::concepts::discrete_dimension QueryGrid1D,
        ddc::concepts::discrete_domain IdxRangeHead,
        class... IdxRangeTypes>
struct FindRelevantIdxRangeType<QueryGrid1D, std::tuple<IdxRangeHead, IdxRangeTypes...>>
{
    using type = ddc::type_seq_merge_t<
            typename SelectRelevantIdxRangeType<QueryGrid1D, IdxRangeHead>::type,
            typename FindRelevantIdxRangeType<QueryGrid1D, std::tuple<IdxRangeTypes...>>::type>;
};

template <ddc::concepts::type_seq PatchTypeSeq>
struct PatchIdsAreUnique;

template <concepts::Patch... Patches>
struct PatchIdsAreUnique<ddc::TypeSeq<Patches...>>
{
    static constexpr bool value = []() {
        std::array<std::size_t, sizeof...(Patches)> const ids {Patches::id...};
        for (std::size_t i = 0; i < ids.size(); ++i) {
            for (std::size_t j = i + 1; j < ids.size(); ++j) {
                if (ids[i] == ids[j]) {
                    return false;
                }
            }
        }
        return true;
    }();
};

} // end namespace connectivity_details

template <ddc::concepts::type_seq PatchTypeSeq>
inline constexpr bool patch_ids_are_unique_v
        = connectivity_details::PatchIdsAreUnique<PatchTypeSeq>::value;

template <ddc::concepts::type_seq TypeSeq>
using to_tuple_t = typename connectivity_details::ToTuple<TypeSeq>::type;

template <concepts::Edge... EdgeType>
using strip_outside_edges_t =
        typename connectivity_details::StripOutsideEdges<ddc::TypeSeq<EdgeType...>>::type;

template <ddc::concepts::type_seq EdgeTypeSeq>
using extract_patches_t = typename connectivity_details::ExtractPatches<EdgeTypeSeq>::type;

template <concepts::Patch StartPatch, ddc::concepts::type_seq InterfaceTypeSeq>
using interfaces_of_patch_t =
        typename connectivity_details::PatchConnection<StartPatch, InterfaceTypeSeq>::type;

template <ddc::concepts::discrete_dimension Grid1D, ddc::concepts::type_seq PatchTypeSeq>
using find_patch_t = typename connectivity_details::FindPatch<Grid1D, PatchTypeSeq>::type;

template <concepts::Edge EdgeType, ddc::concepts::type_seq InterfaceTypeSeq>
using find_associated_interface_t =
        typename connectivity_details::FindInterface<EdgeType, InterfaceTypeSeq>::type;

template <
        concepts::Patch StartPatch,
        ddc::concepts::discrete_dimension Grid1D,
        ddc::concepts::type_seq InterfaceTypeSeq>
using collect_grids_on_dim_t = typename connectivity_details::
        CollectAllGridsOnDim<StartPatch, Grid1D, InterfaceTypeSeq>::type;

template <
        concepts::Patch StartPatch,
        ddc::concepts::discrete_dimension Grid1D,
        ddc::concepts::type_seq InterfaceTypeSeq>
using collect_interfaces_on_dim_t = typename connectivity_details::
        CollectAllInterfacesOnDim<StartPatch, Grid1D, InterfaceTypeSeq>::type;

template <ddc::concepts::discrete_dimension QueryGrid1D, class IdxRangeTuple>
using find_relevant_idx_range_t = ddc::type_seq_element_t<
        0,
        typename connectivity_details::FindRelevantIdxRangeType<QueryGrid1D, IdxRangeTuple>::type>;
```


