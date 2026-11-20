// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "geometry_descriptors.hpp"
#include "patch.hpp"


/**
 * @brief Define an edge of a given patch.
 * 
 * An edge is defined by a patch, a dimension and an extremity.
 * For example, in the patch defined on logical index range 
 * @f$ [a_x, b_x]\times[a_y, b_y] @f$, 
 * 
 * * the edge GridY, BACK refers to the set @f$ [a_x, b_x]\times\{b_y\} @f$, 
 * * and the edge GridY, FRONT refers to the set @f$ [a_x, b_x]\times\{a_y\} @f$. 
 * 
 * @tparam Patch Patch where the edge is defined. 
 * @tparam Grid1D Grid on the complementary dimension of the edge.     
 * @tparam extremity_val The BACK or FRONT value.
*/
template <concepts::Patch Patch, class Grid1D, Extremity extremity_val>
struct Edge
{
    /// @brief Patch where the edge is defined.
    using associated_patch = Patch;
    /// @brief Grid on the perpendicular dimension of the edge.
    using perpendicular_grid = Grid1D;
    /// @brief Grid parallel to the edge.
    using parallel_grid = std::conditional_t<
            std::is_same_v<Grid1D, typename Patch::Grid1>,
            typename Patch::Grid2,
            typename Patch::Grid1>;
    /// @brief Design if the edge is on the BACK or the FRONT of the other dimension.
    static constexpr Extremity extremity = extremity_val;
};

/**
 * @brief Define an edge for the outside index range.
 * OutsideEdge is a pseudo-edge outside the index range
 * used to define interfaces between patches and the 
 * outside index range.
*/
struct OutsideEdge
{
};


/**
 * @brief A type trait to check if a type is an edge (an Edge or an OutsideEdge).
 * @tparam EdgeType The type being checked.
 */
template <class EdgeType>
struct is_edge : std::false_type
{
};

/// @brief Specialisation of is_edge for an Edge of a patch.
template <concepts::Patch Patch, class Grid1D, Extremity extremity_val>
struct is_edge<Edge<Patch, Grid1D, extremity_val>> : std::true_type
{
};

/// @brief Specialisation of is_edge for the OutsideEdge.
template <>
struct is_edge<OutsideEdge> : std::true_type
{
};

/**
 * @brief A boolean indicating if a type is an edge (an Edge or an OutsideEdge).
 * @tparam EdgeType The type being checked.
 */
template <class EdgeType>
bool constexpr is_edge_v = is_edge<EdgeType>::value;

/**
 * @brief Get the id of the patch on which an edge is defined.
 *
 * This can be used to obtain a compile-time ordering of the edges.
 *
 * @tparam EdgeType The edge (an Edge or an OutsideEdge).
 * @return The id of the patch associated with the edge, or -1 for an OutsideEdge.
 */
template <class EdgeType>
constexpr int get_edge_patch_id()
{
    if constexpr (std::is_same_v<EdgeType, OutsideEdge>) {
        return -1;
    } else {
        return EdgeType::associated_patch::id;
    }
}

namespace concepts {
/**
 * @brief A concept describing an edge (an Edge or an OutsideEdge).
 * @tparam EdgeType The type being checked.
 */
template <class EdgeType>
concept Edge = is_edge_v<EdgeType>;
}; // namespace concepts
