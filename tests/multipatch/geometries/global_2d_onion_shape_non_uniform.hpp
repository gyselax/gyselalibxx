// SPDX-License-Identifier: MIT

/*
    Geometry defined here: global domain.
        - non-periodic on Rg and periodic Thetag;
        - non uniform cubic splines;
        - non uniform GridRg and GridThetag; 
*/
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"

namespace onion_shape_non_uniform_2d_global { // TODO: name it Global to be able to use it as the Patches?

int constexpr BSplineDegree = 3;

// CONTINUOUS DIMENSIONS -------------------------------------------------------------------------
/**
 * @brief First continuous dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
struct Rg
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/**
 * @brief Second continuous dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
struct Thetag
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = true;
};



// DISCRETE DIMENSIONS ---------------------------------------------------------------------------
/**
 * @brief Points sequence on the second logical dimension of patch 1.
 */
struct GridRg : NonUniformGridBase<Rg>
{
};

/**
 * @brief Points sequence on the second logical dimension of patch 1.
 */
struct GridThetag : NonUniformGridBase<Thetag>
{
};



// SPLINE DIMENSIONS -----------------------------------------------------------------------------
/**
 * @brief First spline dimension of patch PatchIdx.
 */
struct BSplinesRg : ddc::NonUniformBSplines<Rg, BSplineDegree>
{
};

/**
 * @brief Second spline dimension of patch PatchIdx.
 */
struct BSplinesThetag : ddc::NonUniformBSplines<Thetag, BSplineDegree>
{
};

} // namespace onion_shape_non_uniform_2d_global