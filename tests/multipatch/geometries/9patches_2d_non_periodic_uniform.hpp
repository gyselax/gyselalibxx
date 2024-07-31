// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 9 * 2D patches.
        - for all patches: 
            - non-periodic on X and Y;
            - uniform cubic splines;
            - uniform Grid1 and Grid2; 
*/

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "patch.hpp"


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
