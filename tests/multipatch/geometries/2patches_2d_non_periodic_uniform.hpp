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

#include "connectivity.hpp"
#include "edge.hpp"
#include "interface.hpp"
#include "patch.hpp"

namespace non_periodic_uniform_2d_2patches {

namespace non_periodic_uniform_2d_2patches {

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

using NorthEdge1 = Edge<Patch1, GridY1, BACK>;
using SouthEdge1 = Edge<Patch1, GridY1, FRONT>;
using EastEdge1 = Edge<Patch1, GridX1, BACK>;
using WestEdge1 = Edge<Patch1, GridX1, FRONT>;

using NorthEdge2 = Edge<Patch2, GridY2, BACK>;
using SouthEdge2 = Edge<Patch2, GridY2, FRONT>;
using EastEdge2 = Edge<Patch2, GridX2, BACK>;
using WestEdge2 = Edge<Patch2, GridX2, FRONT>;

using NorthInterface1 = Interface<NorthEdge1, OutsideEdge, true>;
using SouthInterface1 = Interface<SouthEdge1, OutsideEdge, true>;
using WestInterface1 = Interface<WestEdge1, OutsideEdge, true>;

using NorthInterface2 = Interface<NorthEdge2, OutsideEdge, true>;
using SouthInterface2 = Interface<SouthEdge2, OutsideEdge, true>;
using EastInterface2 = Interface<EastEdge2, OutsideEdge, true>;

using Interface_1_2 = Interface<EastEdge1, WestEdge2, true>;

using Connectivity = MultipatchConnectivity<
        NorthInterface1,
        SouthInterface1,
        WestInterface1,
        NorthInterface2,
        SouthInterface2,
        EastInterface2,
        Interface_1_2>;

} // namespace non_periodic_uniform_2d_2patches
