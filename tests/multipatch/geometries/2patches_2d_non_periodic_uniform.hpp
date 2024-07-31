// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 2 * 2D patches.
        - patch 1 and patch 2: 
            - non-periodic on X and Y;
            - uniform cubic splines;
            - uniform Grid1 and Grid2; 
*/

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "patch.hpp"


int constexpr BSplineDegree = 3;

// CONTINUOUS DIMENSIONS -------------------------------------------------------------------------
/// @brief First continuous dimension of patch 1.
struct X1
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/// @brief Second continuous dimension of patch 1.
struct Y1
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/// @brief First continuous dimension of patch 2.
struct X2
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/// @brief Second continuous dimension of patch 2.
struct Y2
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};


// GRID: points sequences ------------------------------------------------------------------------
/// @brief Points sequence on the first logical dimension of patch 1.
struct GridX1 : ddc::UniformPointSampling<X1>
{
};
/// @brief Points sequence on the second logical dimension of patch 1.
struct GridY1 : ddc::UniformPointSampling<Y1>
{
};
/// @brief Points sequence on the first logical dimension of patch 2.
struct GridX2 : ddc::UniformPointSampling<X2>
{
};
/// @brief Points sequence on the second logical dimension of patch 2.
struct GridY2 : ddc::UniformPointSampling<Y2>
{
};


// SPLINE DIMENSIONS -----------------------------------------------------------------------------
/// @brief First spline dimension of patch 1.
struct BSplinesX1 : ddc::UniformBSplines<X1, BSplineDegree>
{
};
/// @brief Second spline dimension of patch 1.
struct BSplinesY1 : ddc::UniformBSplines<Y1, BSplineDegree>
{
};

/// @brief First spline dimension of patch 2.
struct BSplinesX2 : ddc::UniformBSplines<X2, BSplineDegree>
{
};
/// @brief Second spline dimension of patch 2.
struct BSplinesY2 : ddc::UniformBSplines<Y2, BSplineDegree>
{
};


// PATCHES ---------------------------------------------------------------------------------------
/// @brief First patch.
using Patch1 = Patch<GridX1, GridY1, BSplinesX1, BSplinesY1>;
/// @brief Second patch.
using Patch2 = Patch<GridX2, GridY2, BSplinesX2, BSplinesY2>;

/// @brief Sorted list of patches.
using PatchOrdering = ddc::detail::TypeSeq<Patch1, Patch2>;
