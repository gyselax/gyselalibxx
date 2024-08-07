// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 5 * 2D patches.
        - for all patches: 
            - non-periodic on X and Y;
            - uniform cubic splines;
            - uniform Grid1 and Grid2; 

    The geometry describes a figure of eight (the shape of the number eight) such that:
        - The northern edge of Patch 1 is connected to the western edge of Patch 2
        - The eastern edge of Patch 4 is connected to the southern edge of Patch 5

    The outside describes a circle such that:
        - The X edges (eastern/western) of Patch 1 are connected to the northern edges of Patches 2 and 4
        - The X edges (eastern/western) of Patch 5 are connected to the southern edges of Patches 2 and 4

    → → → ↘
  ↗ 
 ↑     ↙ |  1  |↖ 
 ↑   -----------------
  ↖   2  |  3  |  4
    ----------------- ↖
       ↘ |  5  |↗      ↑
             ↘        ↗
                → → →
*/
#pragma once

#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

namespace figure_of_eight_5patches {

int constexpr BSplineDegree = 3;

// CONTINUOUS DIMENSIONS -------------------------------------------------------------------------
/**
 * @brief First continuous dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the pach. 
 */
template <int PatchIdx>
struct X
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/**
 * @brief Second continuous dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the pach. 
 */
template <int PatchIdx>
struct Y
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};



// DISCRETE DIMENSIONS ---------------------------------------------------------------------------
/**
 * @brief Points sequence on the second logical dimension of patch 1.
 * @tparam PatchIdx Index of the pach. 
 */
template <int PatchIdx>
struct GridX : UniformGridBase<X<PatchIdx>>
{
};

/**
 * @brief Points sequence on the second logical dimension of patch 1.
 * @tparam PatchIdx Index of the pach. 
 */
template <int PatchIdx>
struct GridY : UniformGridBase<Y<PatchIdx>>
{
};



// SPLINE DIMENSIONS -----------------------------------------------------------------------------
/**
 * @brief First spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the pach. 
 */
template <int PatchIdx>
struct BSplinesX : ddc::UniformBSplines<X<PatchIdx>, BSplineDegree>
{
};

/**
 * @brief Second spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the pach. 
 */
template <int PatchIdx>
struct BSplinesY : ddc::UniformBSplines<Y<PatchIdx>, BSplineDegree>
{
};


// PATCHES ---------------------------------------------------------------------------------------
/// @brief First patch.
using Patch1 = Patch<GridX<1>, GridY<1>, BSplinesX<1>, BSplinesY<1>>;

/// @brief Second patch.
using Patch2 = Patch<GridX<2>, GridY<2>, BSplinesX<2>, BSplinesY<2>>;

/// @brief Third patch.
using Patch3 = Patch<GridX<3>, GridY<3>, BSplinesX<3>, BSplinesY<3>>;

/// @brief Fourth patch.
using Patch4 = Patch<GridX<4>, GridY<4>, BSplinesX<4>, BSplinesY<4>>;

/// @brief Fifth patch.
using Patch5 = Patch<GridX<5>, GridY<5>, BSplinesX<5>, BSplinesY<5>>;


/// @brief Sorted list of patches.
using PatchOrdering = ddc::detail::TypeSeq<Patch1, Patch2, Patch3, Patch4, Patch5>;

using NorthEdge1 = Edge<Patch1, GridY<1>, BACK>;
using SouthEdge1 = Edge<Patch1, GridY<1>, FRONT>;
using EastEdge1 = Edge<Patch1, GridX<1>, BACK>;
using WestEdge1 = Edge<Patch1, GridX<1>, FRONT>;

using NorthEdge2 = Edge<Patch2, GridY<2>, BACK>;
using SouthEdge2 = Edge<Patch2, GridY<2>, FRONT>;
using EastEdge2 = Edge<Patch2, GridX<2>, BACK>;
using WestEdge2 = Edge<Patch2, GridX<2>, FRONT>;

using NorthEdge3 = Edge<Patch3, GridY<3>, BACK>;
using SouthEdge3 = Edge<Patch3, GridY<3>, FRONT>;
using EastEdge3 = Edge<Patch3, GridX<3>, BACK>;
using WestEdge3 = Edge<Patch3, GridX<3>, FRONT>;

using NorthEdge4 = Edge<Patch4, GridY<4>, BACK>;
using SouthEdge4 = Edge<Patch4, GridY<4>, FRONT>;
using EastEdge4 = Edge<Patch4, GridX<4>, BACK>;
using WestEdge4 = Edge<Patch4, GridX<4>, FRONT>;

using NorthEdge5 = Edge<Patch5, GridY<5>, BACK>;
using SouthEdge5 = Edge<Patch5, GridY<5>, FRONT>;
using EastEdge5 = Edge<Patch5, GridX<5>, BACK>;
using WestEdge5 = Edge<Patch5, GridX<5>, FRONT>;


using Patch3InterfaceNorth = Interface<NorthEdge3, SouthEdge1, true>;
using Patch3InterfaceSouth = Interface<SouthEdge3, NorthEdge5, false>;
using Patch3InterfaceEast = Interface<EastEdge3, WestEdge4, false>;
using Patch3InterfaceWest = Interface<WestEdge3, EastEdge2, true>;

using LoopInterface_2_1 = Interface<NorthEdge2, WestEdge1, false>;
using LoopInterface_1_4 = Interface<EastEdge1, NorthEdge4, true>;
using LoopInterface_4_5 = Interface<SouthEdge4, EastEdge5, false>;
using LoopInterface_5_2 = Interface<WestEdge5, SouthEdge2, true>;

using EightInterface_2_1 = Interface<WestEdge2, NorthEdge1, false>;
using EightInterface_5_4 = Interface<SouthEdge5, EastEdge4, false>;

using Connectivity = MultipatchConnectivity<
        Patch3InterfaceNorth,
        Patch3InterfaceSouth,
        Patch3InterfaceEast,
        Patch3InterfaceWest,
        LoopInterface_2_1,
        LoopInterface_1_4,
        LoopInterface_4_5,
        LoopInterface_5_2,
        EightInterface_2_1,
        EightInterface_5_4>;

}; // namespace figure_of_eight_5patches
