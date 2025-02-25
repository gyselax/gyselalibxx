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
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct X
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;
    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;
    /// A type-alias mapping to the covariant counterpart.
    using Dual = X;
};

/**
 * @brief Second continuous dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct Y
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
    /// A boolean indicating if dimension describes a covariant coordinate.
    static bool constexpr IS_COVARIANT = true;
    /// A boolean indicating if dimension describes a contravariant coordinate.
    static bool constexpr IS_CONTRAVARIANT = true;
    /// A type-alias mapping to the covariant counterpart.
    using Dual = Y;
};



// DISCRETE DIMENSIONS ---------------------------------------------------------------------------
/**
 * @brief Points sequence on the second logical dimension of patch 1.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct GridX : UniformGridBase<X<PatchIdx>>
{
};

/**
 * @brief Points sequence on the second logical dimension of patch 1.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct GridY : UniformGridBase<Y<PatchIdx>>
{
};



// SPLINE DIMENSIONS -----------------------------------------------------------------------------
/**
 * @brief First spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct BSplinesX : ddc::UniformBSplines<X<PatchIdx>, BSplineDegree>
{
};

/**
 * @brief Second spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
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
using NorthInterface1 = Interface<NorthEdge<1>, OutsideEdge, true>;
using NorthInterface2 = Interface<NorthEdge<2>, OutsideEdge, true>;
using NorthInterface3 = Interface<NorthEdge<3>, OutsideEdge, true>;

using Interface_1_4 = Interface<SouthEdge<1>, NorthEdge<4>, true>;
using Interface_2_5 = Interface<SouthEdge<2>, NorthEdge<5>, true>;
using Interface_3_6 = Interface<SouthEdge<3>, NorthEdge<6>, true>;

using Interface_4_7 = Interface<SouthEdge<4>, NorthEdge<7>, true>;
using Interface_5_8 = Interface<SouthEdge<5>, NorthEdge<8>, true>;
using Interface_6_9 = Interface<SouthEdge<6>, NorthEdge<9>, true>;

using SouthInterface7 = Interface<OutsideEdge, SouthEdge<7>, true>;
using SouthInterface8 = Interface<OutsideEdge, SouthEdge<8>, true>;
using SouthInterface9 = Interface<OutsideEdge, SouthEdge<9>, true>;


using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;
using Interface_4_5 = Interface<EastEdge<4>, WestEdge<5>, true>;
using Interface_7_8 = Interface<EastEdge<7>, WestEdge<8>, true>;

using Interface_2_3 = Interface<EastEdge<2>, WestEdge<3>, true>;
using Interface_5_6 = Interface<EastEdge<5>, WestEdge<6>, true>;
using Interface_8_9 = Interface<EastEdge<8>, WestEdge<9>, true>;

using Interface_3_1 = Interface<WestEdge<1>, EastEdge<3>, true>;
using Interface_6_4 = Interface<WestEdge<4>, EastEdge<6>, true>;
using Interface_9_7 = Interface<WestEdge<7>, EastEdge<9>, true>;


// CONNECTIVITY ----------------------------------------------------------------------------------
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
