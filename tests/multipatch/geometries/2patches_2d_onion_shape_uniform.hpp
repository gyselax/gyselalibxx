// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 2 * 2D patches.
        - patch 1 and patch 2: 
            - non-periodic on R and periodic Theta;
            - uniform cubic splines;
            - uniform GridR and GridTheta; 
*/
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "edge.hpp"
#include "patch.hpp"

namespace onion_shape_uniform_2d_2patches {


int constexpr BSplineDegree = 3;

// CONTINUOUS DIMENSIONS -------------------------------------------------------------------------
/**
 * @brief First continuous dimension of all the patches.
 */
struct R
{
    /// @brief Non periodic dimension.
    static bool constexpr PERIODIC = false;
};

/**
 * @brief Second continuous dimension of all the patches.
 */
struct Theta
{
    /// @brief Periodic dimension.
    static bool constexpr PERIODIC = true;
};



// DISCRETE DIMENSIONS ---------------------------------------------------------------------------
/**
 * @brief Points sequence on the second logical dimension of patch 1.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct GridR : UniformGridBase<R>
{
};

/**
 * @brief Points sequence on the second logical dimension of patch 1.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct GridTheta : UniformGridBase<Theta>
{
};



// SPLINE DIMENSIONS -----------------------------------------------------------------------------
/**
 * @brief First spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct BSplinesR : ddc::UniformBSplines<R, BSplineDegree>
{
};

/**
 * @brief Second spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct BSplinesTheta : ddc::UniformBSplines<Theta, BSplineDegree>
{
};


// PATCHES ---------------------------------------------------------------------------------------
/// @brief First patch.
using Patch1 = Patch<GridR<1>, GridTheta<1>, BSplinesR<1>, BSplinesTheta<1>>;

/// @brief Second patch.
using Patch2 = Patch<GridR<2>, GridTheta<2>, BSplinesR<2>, BSplinesTheta<2>>;

/// @brief Sorted list of patches.
using PatchOrdering = ddc::detail::TypeSeq<Patch1, Patch2>;


// EDGES -----------------------------------------------------------------------------------------
template <int PatchIdx>
using NorthEdge
        = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
               GridTheta<PatchIdx>,
               BACK>;
template <int PatchIdx>
using SouthEdge
        = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
               GridTheta<PatchIdx>,
               FRONT>;
template <int PatchIdx>
using EastEdge
        = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
               GridR<PatchIdx>,
               BACK>;
template <int PatchIdx>
using WestEdge
        = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
               GridR<PatchIdx>,
               FRONT>;


// INTERFACES ------------------------------------------------------------------------------------


// CONNECTIVITY ----------------------------------------------------------------------------------


} // namespace onion_shape_uniform_2d_2patches