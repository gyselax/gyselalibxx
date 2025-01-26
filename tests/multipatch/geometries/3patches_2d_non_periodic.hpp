// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 3 * 2D patches.
        - for all patches: 
            - non-periodic on X and Y;
            - uniform/non-uniform cubic splines;
            - uniform/non-uniform Grid1 and Grid2; 
*/
#pragma once

#define EXTEND_GEOMETRY
#include "2patches_2d_non_periodic.hpp"

namespace GEOM_NAMESPACE_NAME {
// PATCHES ---------------------------------------------------------------------------------------
/// @brief Third patch.
using Patch3 = Patch<GridX<3>, GridY<3>, BSplinesX<3>, BSplinesY<3>>;

/// @brief Sorted list of patches.
using PatchOrdering = ddc::detail::TypeSeq<Patch1, Patch2, Patch3>;


// EDGES -----------------------------------------------------------------------------------------
template <int PatchIdx>
using NorthEdge
        = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
               GridY<PatchIdx>,
               BACK>;
template <int PatchIdx>
using SouthEdge
        = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
               GridY<PatchIdx>,
               FRONT>;
template <int PatchIdx>
using EastEdge
        = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
               GridX<PatchIdx>,
               BACK>;
template <int PatchIdx>
using WestEdge
        = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
               GridX<PatchIdx>,
               FRONT>;


// INTERFACES ------------------------------------------------------------------------------------


// CONNECTIVITY ----------------------------------------------------------------------------------

} // namespace GEOM_NAMESPACE_NAME
#undef EXTEND_GEOMETRY
