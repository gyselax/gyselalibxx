// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 3 * 2D patches.
        - for all patches: 
            - non-periodic on X and Y;
            - uniform/non-uniform cubic splines;
            - uniform/non-uniform Grid1 and Grid2; 
*/
#ifdef GEOM_NAMESPACE_NAME
#define EXTEND_GEOMETRY
#include "2patches_2d_non_periodic.hpp"

namespace GEOM_NAMESPACE_NAME {
// PATCHES ---------------------------------------------------------------------------------------
/// @brief Third patch.
using Patch3 = PatchType<3>;

/// @brief Sorted list of patches.
using PatchOrdering = ddc::TypeSeq<Patch1, Patch2, Patch3>;


// INTERFACES ------------------------------------------------------------------------------------


// CONNECTIVITY ----------------------------------------------------------------------------------

} // namespace GEOM_NAMESPACE_NAME
#undef EXTEND_GEOMETRY
#endif
