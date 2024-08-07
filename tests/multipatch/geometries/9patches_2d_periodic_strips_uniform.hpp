// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 9 * 2D patches.
        - for all patches: 
            - non-periodic on X,Y;
            - uniform cubic splines;
            - uniform Grid1 and Grid2; 
        - across patches:
            - periodic on X;
            - non-periodic on Y;


      1  |  2  |  3
    -----------------
      4  |  5  |  6
    -----------------
      7  |  8  |  9
*/
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "connectivity.hpp"
#include "edge.hpp"
#include "interface.hpp"
#include "patch.hpp"

namespace periodic_strips_uniform_2d_9patches {

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

/// @brief Sixth patch.
using Patch6 = Patch<GridX<6>, GridY<6>, BSplinesX<6>, BSplinesY<6>>;

/// @brief Seventh patch.
using Patch7 = Patch<GridX<7>, GridY<7>, BSplinesX<7>, BSplinesY<7>>;

/// @brief Eighth patch.
using Patch8 = Patch<GridX<8>, GridY<8>, BSplinesX<8>, BSplinesY<8>>;

/// @brief Ninth patch.
using Patch9 = Patch<GridX<9>, GridY<9>, BSplinesX<9>, BSplinesY<9>>;


/// @brief Sorted list of patches.
using PatchOrdering = ddc::detail::
        TypeSeq<Patch1, Patch2, Patch3, Patch4, Patch5, Patch6, Patch7, Patch8, Patch9>;

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

using NorthEdge6 = Edge<Patch6, GridY<6>, BACK>;
using SouthEdge6 = Edge<Patch6, GridY<6>, FRONT>;
using EastEdge6 = Edge<Patch6, GridX<6>, BACK>;
using WestEdge6 = Edge<Patch6, GridX<6>, FRONT>;

using NorthEdge7 = Edge<Patch7, GridY<7>, BACK>;
using SouthEdge7 = Edge<Patch7, GridY<7>, FRONT>;
using EastEdge7 = Edge<Patch7, GridX<7>, BACK>;
using WestEdge7 = Edge<Patch7, GridX<7>, FRONT>;

using NorthEdge8 = Edge<Patch8, GridY<8>, BACK>;
using SouthEdge8 = Edge<Patch8, GridY<8>, FRONT>;
using EastEdge8 = Edge<Patch8, GridX<8>, BACK>;
using WestEdge8 = Edge<Patch8, GridX<8>, FRONT>;

using NorthEdge9 = Edge<Patch9, GridY<9>, BACK>;
using SouthEdge9 = Edge<Patch9, GridY<9>, FRONT>;
using EastEdge9 = Edge<Patch9, GridX<9>, BACK>;
using WestEdge9 = Edge<Patch9, GridX<9>, FRONT>;


using NorthInterface1 = Interface<NorthEdge1, OutsideEdge, true>;
using NorthInterface2 = Interface<NorthEdge2, OutsideEdge, true>;
using NorthInterface3 = Interface<NorthEdge3, OutsideEdge, true>;

using Interface_1_4 = Interface<SouthEdge1, NorthEdge4, true>;
using Interface_2_5 = Interface<SouthEdge2, NorthEdge5, true>;
using Interface_3_6 = Interface<SouthEdge3, NorthEdge6, true>;

using Interface_4_7 = Interface<SouthEdge4, NorthEdge7, true>;
using Interface_5_8 = Interface<SouthEdge5, NorthEdge8, true>;
using Interface_6_9 = Interface<SouthEdge6, NorthEdge9, true>;

using SouthInterface7 = Interface<OutsideEdge, SouthEdge7, true>;
using SouthInterface8 = Interface<OutsideEdge, SouthEdge8, true>;
using SouthInterface9 = Interface<OutsideEdge, SouthEdge9, true>;


using Interface_1_2 = Interface<EastEdge1, WestEdge2, true>;
using Interface_4_5 = Interface<EastEdge4, WestEdge5, true>;
using Interface_7_8 = Interface<EastEdge7, WestEdge8, true>;

using Interface_2_3 = Interface<EastEdge2, WestEdge3, true>;
using Interface_5_6 = Interface<EastEdge5, WestEdge6, true>;
using Interface_8_9 = Interface<EastEdge8, WestEdge9, true>;

using Interface_3_1 = Interface<WestEdge1, EastEdge3, true>;
using Interface_6_4 = Interface<WestEdge4, EastEdge6, true>;
using Interface_9_7 = Interface<WestEdge7, EastEdge9, true>;

using Connectivity = MultipatchConnectivity<
        NorthInterface1,
        NorthInterface2,
        NorthInterface3,
        Interface_1_4,
        Interface_2_5,
        Interface_3_6,
        Interface_4_7,
        Interface_5_8,
        Interface_6_9,
        SouthInterface7,
        SouthInterface8,
        SouthInterface9,
        Interface_1_2,
        Interface_4_5,
        Interface_7_8,
        Interface_2_3,
        Interface_5_6,
        Interface_8_9,
        Interface_3_1,
        Interface_6_4,
        Interface_9_7>;

} // namespace periodic_strips_uniform_2d_9patches
