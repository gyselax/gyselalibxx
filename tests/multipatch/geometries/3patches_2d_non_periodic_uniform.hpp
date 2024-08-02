// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 3 * 2D patches.
        - for all patches: 
            - non-periodic on X and Y;
            - uniform cubic splines;
            - uniform Grid1 and Grid2; 
*/

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "patch.hpp"


namespace non_periodic_uniform_2d_3patches {

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

/// @brief First continuous dimension of patch 3.
struct X3
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/// @brief Second continuous dimension of patch 3.
struct Y3
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};



// GRID: points sequences ------------------------------------------------------------------------
/// @brief Points sequence on the first logical dimension of patch 1.
struct GridX1 : UniformGridBase<X1>
{
};
/// @brief Points sequence on the second logical dimension of patch 1.
struct GridY1 : UniformGridBase<Y1>
{
};
/// @brief Points sequence on the first logical dimension of patch 2.
struct GridX2 : UniformGridBase<X2>
{
};
/// @brief Points sequence on the second logical dimension of patch 2.
struct GridY2 : UniformGridBase<Y2>
{
};
/// @brief Points sequence on the first logical dimension of patch 3.
struct GridX3 : UniformGridBase<X3>
{
};
/// @brief Points sequence on the second logical dimension of patch 3.
struct GridY3 : UniformGridBase<Y3>
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

/// @brief First spline dimension of patch 2.
struct BSplinesX3 : ddc::UniformBSplines<X3, BSplineDegree>
{
};
/// @brief Second spline dimension of patch 2.
struct BSplinesY3 : ddc::UniformBSplines<Y3, BSplineDegree>
{
};


// PATCHES ---------------------------------------------------------------------------------------
/// @brief First patch.
using Patch1 = Patch<GridX1, GridY1, BSplinesX1, BSplinesY1>;
/// @brief Second patch.
using Patch2 = Patch<GridX2, GridY2, BSplinesX2, BSplinesY2>;
/// @brief Third patch.
using Patch3 = Patch<GridX3, GridY3, BSplinesX3, BSplinesY3>;

/// @brief Sorted list of patches.
using PatchOrdering = ddc::detail::TypeSeq<Patch1, Patch2, Patch3>;


} // namespace non_periodic_uniform_2d_3patches