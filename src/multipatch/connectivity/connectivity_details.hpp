// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "interface.hpp"
#include "patch.hpp"

namespace connectivity_details {
/**
 * @brief A class which finds all interfaces connected to a given patch.
 *
 * @tparam Patch The patch whose connections are being searched for.
 * @tparam InterfaceTypeSeq A DDC type sequence containing all the possible Interfaces.
 */
template <class Patch, class InterfaceTypeSeq>
struct PatchConnection;

/// Specialisation of PatchConnection for an empty interface list.
template <class Patch>
struct PatchConnection<Patch, ddc::detail::TypeSeq<>>
{
    /// The type found by the class.
    using type = ddc::detail::TypeSeq<>;
};

/// Specialisation of PatchConnection for an interface list with one element.
template <class Patch, class InterfaceType>
struct PatchConnection<Patch, ddc::detail::TypeSeq<InterfaceType>>
{
    /// The type found by the class.
    using type = std::conditional_t<
            InterfaceType::template connected_to_patch<Patch>(),
            ddc::detail::TypeSeq<InterfaceType>,
            ddc::detail::TypeSeq<>>;
};

/// Specialisation of PatchConnection to iterate recursively over the interface type sequence.
template <class Patch, class InterfaceType1, class... RemainingInterfaceTypes>
struct PatchConnection<Patch, ddc::detail::TypeSeq<InterfaceType1, RemainingInterfaceTypes...>>
{
    /// The type found by the class.
    using type = ddc::type_seq_merge_t<
            typename PatchConnection<Patch, ddc::detail::TypeSeq<InterfaceType1>>::type,
            typename PatchConnection<Patch, ddc::detail::TypeSeq<RemainingInterfaceTypes...>>::
                    type>;
};

/**
 * @brief A class which finds all edges which are not OutsideEdge types.
 *
 * @tparam A DDC type sequence containing all the possible Edges.
 */
template <class TypeSeq>
struct StripOutsideEdges;

/// Specialisation of StripOutsideEdges for the case with one edge in the list.
template <class EdgeType>
struct StripOutsideEdges<ddc::detail::TypeSeq<EdgeType>>
{
    /// The type found by the class.
    using type = std::conditional_t<
            std::is_same_v<EdgeType, OutsideEdge>,
            ddc::detail::TypeSeq<>,
            ddc::detail::TypeSeq<EdgeType>>;
};

/// Specialisation of StripOutsideEdges to iterate recursively over the edge type sequence.
template <class EdgeType1, class... RemainingEdgeTypes>
struct StripOutsideEdges<ddc::detail::TypeSeq<EdgeType1, RemainingEdgeTypes...>>
{
    /// The type found by the class.
    using type = ddc::type_seq_merge_t<
            typename StripOutsideEdges<ddc::detail::TypeSeq<EdgeType1>>::type,
            typename StripOutsideEdges<ddc::detail::TypeSeq<RemainingEdgeTypes...>>::type>;
};

/**
 * @brief A class to convert a type sequence to a tuple type.
 * @tparam The type sequence.
 */
template <class TypeSeq>
struct ToTuple;

/// Specialisation of ToTuple for type sequences.
template <class... I>
struct ToTuple<ddc::detail::TypeSeq<I...>>
{
    /// The type found by the class.
    using type = std::tuple<I...>;
};

/**
 * @brief A class to find all the patches used by the various edges.
 * @tparam A DDC type sequence containing all the Edges.
 */
template <class TypeSeq>
struct ExtractPatches;

/// Specialisation of ExtractPatches for an empty patch list.
template <>
struct ExtractPatches<ddc::detail::TypeSeq<>>
{
    /// The type found by the class.
    using type = ddc::detail::TypeSeq<>;
};

/// Specialisation of ExtractPatches to iterate recursively over the edge type sequence.
template <class EdgeType1, class... EdgeTypes>
struct ExtractPatches<ddc::detail::TypeSeq<EdgeType1, EdgeTypes...>>
{
    /// The type found by the class.
    using type = ddc::type_seq_merge_t<
            ddc::detail::TypeSeq<typename EdgeType1::associated_patch>,
            typename ExtractPatches<ddc::detail::TypeSeq<EdgeTypes...>>::type>;
};

/**
 * @brief A class to locate a patch which contains the specified grid.
 * @tparam The identifying grid.
 * @tparam A DDC type sequence containing all possible patches.
 */
template <class Grid1D, class PatchTypeSeq>
struct FindPatch;

/// Specialisation of FindPatch for an empty patch list.
template <class Grid1D>
struct FindPatch<Grid1D, ddc::detail::TypeSeq<>>
{
    // This class should not be reached in the recursion. Raise a readable error (based on the
    // template arguments) if it is instantiated.
    static_assert(!std::is_same_v<Grid1D, Grid1D>, "No patches found using dimension.");
};

/**
 * Specialisation of FindPatch for the case where Grid1 from first patch in type sequence
 * matches QueryGrid1D
 */
template <class QueryGrid1D, class OGrid, class BSpl1, class BSpl2, class... RemainingPatchTypes>
struct FindPatch<
        QueryGrid1D,
        ddc::detail::TypeSeq<Patch<QueryGrid1D, OGrid, BSpl1, BSpl2>, RemainingPatchTypes...>>
{
    /// The type found by the class.
    using type = Patch<QueryGrid1D, OGrid, BSpl1, BSpl2>;
};

/**
 * Specialisation of FindPatch for the case where Grid2 from first patch in type sequence
 * matches QueryGrid1D
 */
template <class QueryGrid1D, class OGrid, class BSpl1, class BSpl2, class... RemainingPatchTypes>
struct FindPatch<
        QueryGrid1D,
        ddc::detail::TypeSeq<Patch<OGrid, QueryGrid1D, BSpl1, BSpl2>, RemainingPatchTypes...>>
{
    /// The type found by the class.
    using type = Patch<OGrid, QueryGrid1D, BSpl1, BSpl2>;
};

/// Specialisation of FindPatch to iterate recursively over the patch type sequence.
template <class Grid1D, class Patch1, class... RemainingPatchTypes>
struct FindPatch<Grid1D, ddc::detail::TypeSeq<Patch1, RemainingPatchTypes...>>
{
    // Should instantiate specialisations above instead of this specialisation
    static_assert(!std::is_same_v<typename Patch1::Grid1, Grid1D>);
    static_assert(!std::is_same_v<typename Patch1::Grid2, Grid1D>);
    /// The type found by the class.
    using type = typename FindPatch<Grid1D, ddc::detail::TypeSeq<RemainingPatchTypes...>>::type;
};

/**
 * @brief A class to locate an interface which contains the specified edge.
 * @tparam The edge being located.
 * @tparam A DDC type sequence of possible interfaces.
 */
template <class Edge, class InterfaceTypeSeq>
struct FindInterface;

/// Specialisation of FindInterface for an empty interface list.
template <class Edge>
struct FindInterface<Edge, ddc::detail::TypeSeq<>>
{
    // This class should not be reached in the recursion. Raise a readable error (based on the
    // template arguments) if it is instantiated.
    static_assert(!std::is_same_v<Edge, Edge>, "Edge not found in collection of interfaces.");
};

/// Specialisation of FindInterface to iterate recursively over the interface type sequence.
template <class Edge, class Interface1, class... RemainingInterfaceTypes>
struct FindInterface<Edge, ddc::detail::TypeSeq<Interface1, RemainingInterfaceTypes...>>
{
    // Should instantiate specialisations below instead of this specialisation
    static_assert(!std::is_same_v<Edge, typename Interface1::Edge1>);
    static_assert(!std::is_same_v<Edge, typename Interface1::Edge2>);
    /// The type found by the class.
    using type =
            typename FindInterface<Edge, ddc::detail::TypeSeq<RemainingInterfaceTypes...>>::type;
};

/// Specialisation of FindInterface for the case where Edge1 from the first interface matches Edge.
template <class Edge, class OEdge, bool Orientations, class... RemainingInterfaceTypes>
struct FindInterface<
        Edge,
        ddc::detail::TypeSeq<Interface<Edge, OEdge, Orientations>, RemainingInterfaceTypes...>>
{
    /// The type found by the class.
    using type = Interface<Edge, OEdge, Orientations>;
};

/// Specialisation of FindInterface for the case where Edge1 from the second interface matches Edge.
template <class Edge, class OEdge, bool Orientations, class... RemainingInterfaceTypes>
struct FindInterface<
        Edge,
        ddc::detail::TypeSeq<Interface<OEdge, Edge, Orientations>, RemainingInterfaceTypes...>>
{
    /// The type found by the class.
    using type = Interface<OEdge, Edge, Orientations>;
};

/**
 * @brief A class to get the opposite edge of a grid line from one of the edges.
 * @tparam EdgeType The start edge connected to a grid line.
 */
template <class EdgeType>
struct SwapExtremity;

/// Specialisation of SwapExtremity for an edge at the front end of a grid line.
template <class Patch, class Grid1D>
struct SwapExtremity<Edge<Patch, Grid1D, FRONT>>
{
    /// The type found by the class.
    using type = Edge<Patch, Grid1D, BACK>;
};

/// Specialisation of SwapExtremity for an edge at the back end of a grid line.
template <class Patch, class Grid1D>
struct SwapExtremity<Edge<Patch, Grid1D, BACK>>
{
    /// The type found by the class.
    using type = Edge<Patch, Grid1D, FRONT>;
};

/**
 * @brief A utility to get the other edge in an interface.
 * @tparam StartEdge The edge whose equivalent we are looking for.
 * @tparam InterfaceTypeSeq The DDC type sequence describing all possible interfaces which might
 *                  contain the start edge.
 */
template <class StartEdge, class InterfaceTypeSeq>
using equivalent_edge_t =
        typename FindInterface<StartEdge, InterfaceTypeSeq>::type::template OtherEdge<StartEdge>;

/**
 * @brief An enumerate to help when inserting elements into a type sequence.
 */
enum InsertPosition { FrontInsert, BackInsert };

/**
 * @brief A class which helps insert an element into a type sequence.
 * @tparam ToInsert The element to be inserted into the type sequence.
 * @tparam TypeSeq A type sequence into which elements will be inserted.
 * @tparam insert_pos The position where the element should be inserted (back/front).
 */
template <class ToInsert, class TypeSeq, InsertPosition insert_pos>
struct AddToTypeSeq;

/// Specialisation of AddToTypeSeq to add an element at the front of the type sequence.
template <class ToInsert, class TypeSeq>
struct AddToTypeSeq<ToInsert, TypeSeq, FrontInsert>
{
    /// The type found by the class.
    using type = ddc::type_seq_merge_t<TypeSeq, ddc::detail::TypeSeq<ToInsert>>;
};

/// Specialisation of AddToTypeSeq to add an element at the back of the type sequence.
template <class ToInsert, class TypeSeq>
struct AddToTypeSeq<ToInsert, TypeSeq, BackInsert>
{
    /// The type found by the class.
    using type = ddc::type_seq_merge_t<ddc::detail::TypeSeq<ToInsert>, TypeSeq>;
};

/**
 * @brief A class which collects grids along a given dimension on a specified direction
 *          from a starting edge.
 *
 * @tparam StartEdge The edge from which the collection should begin.
 * @tparam InterfaceTypeSeq A DDC type sequence containing all the possible Interfaces.
 * @tparam insert_pos The position where the element should be inserted (back/front).
 * @tparam FoundGrids The grids that have been discovered so far along the given
 *                      dimension and direction.
 */
template <
        class StartEdge,
        class InterfaceTypeSeq,
        InsertPosition insert_pos,
        class FoundGrids = ddc::detail::TypeSeq<>,
        class MatchingEdge = equivalent_edge_t<StartEdge, InterfaceTypeSeq>,
        bool grid_already_found // Periodic case
        = ddc::in_tags_v<typename StartEdge::perpendicular_grid, FoundGrids>>
struct CollectGridsAlongDim;

/// Specialisation of CollectGridsAlongDim to iterate recursively over the grids on the dimension.
template <
        class StartEdge,
        class InterfaceTypeSeq,
        InsertPosition insert_pos,
        class FoundGrids,
        class MatchingEdge>
struct CollectGridsAlongDim<
        StartEdge,
        InterfaceTypeSeq,
        insert_pos,
        FoundGrids,
        MatchingEdge,
        false>
{
    /// The new list of grids that have been found including the grid from the current patch.
    using NewGridList =
            typename AddToTypeSeq<typename StartEdge::perpendicular_grid, FoundGrids, insert_pos>::
                    type;
    /// The type found by the class.
    using type = typename CollectGridsAlongDim<
            typename SwapExtremity<MatchingEdge>::type,
            InterfaceTypeSeq,
            insert_pos,
            NewGridList>::type;
};

/// Specialisation of CollectGridsAlongDim to stop when the grid has already been identified (due to periodicity).
template <
        class StartEdge,
        class InterfaceTypeSeq,
        InsertPosition insert_pos,
        class FoundGrids,
        class MatchingEdge>
struct CollectGridsAlongDim<StartEdge, InterfaceTypeSeq, insert_pos, FoundGrids, MatchingEdge, true>
{
    /// The type found by the class.
    using type = FoundGrids;
};

/// Specialisation of CollectGridsAlongDim to stop when there are no more grids.
template <class StartEdge, class InterfaceTypeSeq, InsertPosition insert_pos, class FoundGrids>
struct CollectGridsAlongDim<StartEdge, InterfaceTypeSeq, insert_pos, FoundGrids, OutsideEdge, false>
{
    /// The type found by the class.
    using type =
            typename AddToTypeSeq<typename StartEdge::perpendicular_grid, FoundGrids, insert_pos>::
                    type;
};

/**
 * @brief A class which collects all grids along a given dimension in both directions.
 *
 * @tparam StartPatch The patch from which the collection should begin.
 * @tparam Grid1D The first grid to be included (this describes the dimension along which
 *                  grids are collected).
 * @tparam InterfaceTypeSeq A DDC type sequence containing all the possible Interfaces.
 */
template <class StartPatch, class Grid1D, class InterfaceTypeSeq>
struct CollectAllGridsOnDim
{
    /**
     * @brief The type sequence describing all grids found by iterating along this
     * dimension in the backwards direction.
     *
     * This is found by working backward from front (start) of grid inserting each
     * new grid at the start of the sequence.
     */
    using BackwardTypeSeq = typename CollectGridsAlongDim<
            Edge<StartPatch, Grid1D, FRONT>,
            InterfaceTypeSeq,
            BackInsert>::type;
    /// The type found by the class.
    using type = ddc::type_seq_merge_t<
            BackwardTypeSeq,
            // Work forward from back (end) of grid inserting each new grid at the end of the sequence
            typename CollectGridsAlongDim<
                    Edge<StartPatch, Grid1D, BACK>,
                    InterfaceTypeSeq,
                    FrontInsert,
                    BackwardTypeSeq>::type>;
};

/**
 * @brief A class to create a type sequence which contains the index range if it can be used to index the grid.
 * @tparam QueryGrid1D The grid which may or may not be present in the index range.
 * @tparam IdxRangeType The index range type being checked.
 */
template <class QueryGrid1D, class IdxRangeType>
struct SelectRelevantIdxRangeType;

/// Specialisation of SelectRelevantIdxRangeType to get access to the grids in the index range.
template <class QueryGrid1D, class... IdxRangeGrids>
struct SelectRelevantIdxRangeType<QueryGrid1D, IdxRange<IdxRangeGrids...>>
{
    /// The type found by the class.
    using type = std::conditional_t<
            ddc::in_tags_v<QueryGrid1D, ddc::detail::TypeSeq<IdxRangeGrids...>>,
            ddc::detail::TypeSeq<IdxRange<IdxRangeGrids...>>,
            ddc::detail::TypeSeq<>>;
};

/**
 * @brief A class to find any index range types which contain an index range defined on the
 *      provided grid.
 *      E.g. Grid1, std::tuple<IdxRange<Grid1, Grid2>, IdxRange<Grid3,Grid4>> will
 *      find: ddc::detail::TypeSeq<IdxRange<Grid1, Grid2>>
 *
 * @tparam QueryGrid1D The grid being searched for.
 * @tparam A tuple of index ranges which may contain the relevant grid.
 */
template <class QueryGrid1D, class IdxRangeTuple>
struct FindRelevantIdxRangeType;

/// Specialisation of FindRelevantIdxRangeType for an empty list of index range types.
template <class QueryGrid1D>
struct FindRelevantIdxRangeType<QueryGrid1D, std::tuple<>>
{
    /// The type found by the class.
    using type = ddc::detail::TypeSeq<>;
};

/// Specialisation of FindRelevantIdxRangeType to iterate recursively over the possible index range types.
template <class QueryGrid1D, class IdxRangeHead, class... IdxRangeTypes>
struct FindRelevantIdxRangeType<QueryGrid1D, std::tuple<IdxRangeHead, IdxRangeTypes...>>
{
    /// The type found by the class.
    using type = ddc::type_seq_merge_t<
            typename SelectRelevantIdxRangeType<QueryGrid1D, IdxRangeHead>::type,
            typename FindRelevantIdxRangeType<QueryGrid1D, std::tuple<IdxRangeTypes...>>::type>;
};

} // end namespace connectivity_details

/// A tool for converting a DDC type sequence to a tuple type.
template <class TypeSeq>
using to_tuple_t = typename connectivity_details::ToTuple<TypeSeq>::type;

/// A tool to find all edges which are not OutsideEdge types.
template <class... EdgeType>
using strip_outside_edges_t =
        typename connectivity_details::StripOutsideEdges<ddc::detail::TypeSeq<EdgeType...>>::type;

/// A tool to find all the patches used by the various edges.
template <class EdgeTypeSeq>
using extract_patches_t = typename connectivity_details::ExtractPatches<EdgeTypeSeq>::type;

/// A tool to find all interfaces directly connected to the start patch.
template <class StartPatch, class InterfaceTypeSeq>
using interfaces_of_patch_t =
        typename connectivity_details::PatchConnection<StartPatch, InterfaceTypeSeq>::type;

/// A tool to find a patch which contains the specified grid.
template <class Grid1D, class PatchTypeSeq>
using find_patch_t = typename connectivity_details::FindPatch<Grid1D, PatchTypeSeq>::type;

/// A tool to find the interface which contains a specified edge.
template <class EdgeType, class InterfaceTypeSeq>
using find_associated_interface_t =
        typename connectivity_details::FindInterface<EdgeType, InterfaceTypeSeq>::type;

/// A tool to collect all grids along a given line including the specified grid.
template <class StartPatch, class Grid1D, class InterfaceTypeSeq>
using collect_grids_on_dim_t = typename connectivity_details::
        CollectAllGridsOnDim<StartPatch, Grid1D, InterfaceTypeSeq>::type;

/// A tool to find the first multi-D index range which contains a specific grid.
template <class QueryGrid1D, class IdxRangeTuple>
using find_relevant_idx_range_t = ddc::type_seq_element_t<
        0,
        typename connectivity_details::FindRelevantIdxRangeType<QueryGrid1D, IdxRangeTuple>::type>;
