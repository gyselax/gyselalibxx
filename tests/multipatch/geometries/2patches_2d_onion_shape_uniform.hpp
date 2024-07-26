// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 2 * 2D patches.
        - patch 1 and patch 2: 
            - non-periodic on R and periodic P;
            - uniform cubic splines;
            - uniform Grid1 and Grid2; 
*/

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "patch.hpp"


int constexpr BSplineDegree = 3;

// CONTINUOUS DIMENSIONS -------------------------------------------------------------------------
/// @brief First continuous dimension of patch 1.
struct R1
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/// @brief Second continuous dimension of patch 1.
struct P1
{
    /// @brief Periodic dimension.
    static bool constexpr PERIODIC = true;
};

/// @brief First continuous dimension of patch 2.
struct R2
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/// @brief Second continuous dimension of patch 2.
struct P2
{
    /// @brief Periodic dimension.
    static bool constexpr PERIODIC = true;
};


// DISCRETE DIMENSIONS ---------------------------------------------------------------------------
/// @brief Points sequence on the first logical dimension of patch 1.
struct GridR1 : ddc::UniformPointSampling<R1>
{
};
/// @brief Points sequence on the second logical dimension of patch 1.
struct GridP1 : ddc::UniformPointSampling<P1>
{
};
/// @brief Points sequence on the first logical dimension of patch 2.
struct GridR2 : ddc::UniformPointSampling<R2>
{
};
/// @brief Points sequence on the second logical dimension of patch 2.
struct GridP2 : ddc::UniformPointSampling<P2>
{
};


// SPLINE DIMENSIONS -----------------------------------------------------------------------------
/// @brief First spline dimension of patch 1.
struct BSplinesR1 : ddc::UniformBSplines<R1, BSplineDegree>
{
};
/// @brief Second spline dimension of patch 1.
struct BSplinesP1 : ddc::UniformBSplines<P1, BSplineDegree>
{
};

/// @brief First spline dimension of patch 2.
struct BSplinesR2 : ddc::UniformBSplines<R2, BSplineDegree>
{
};
/// @brief Second spline dimension of patch 2.
struct BSplinesP2 : ddc::UniformBSplines<P2, BSplineDegree>
{
};


// PATCHES ---------------------------------------------------------------------------------------
/// @brief First patch.
using Patch1 = Patch<GridR1, GridP1, BSplinesR1, BSplinesP1>;
/// @brief Second patch.
using Patch2 = Patch<GridR2, GridP2, BSplinesR2, BSplinesP2>;


/// @brief Sorted list of patches.
using PatchOrdering = ddc::detail::TypeSeq<Patch1, Patch2>;
