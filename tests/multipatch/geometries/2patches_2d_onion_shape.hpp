// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 2 * 2D patches.
        - patch 1 and patch 2: 
            - non-periodic on R and periodic Theta;
            - uniform/non-uniform cubic splines;
            - uniform/non-uniform GridR and GridTheta; 
*/
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "edge.hpp"
#include "patch.hpp"

namespace GEOM_NAMESPACE_NAME {


int constexpr BSplineDegree = 3;

#ifdef UNIFORM
bool constexpr is_uniform = true;
#else
bool constexpr is_uniform = false;
#endif

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
struct GridR : std::conditional_t<is_uniform, UniformGridBase<R>, NonUniformGridBase<R>>
{
};

/**
 * @brief Points sequence on the second logical dimension of patch 1.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct GridTheta : std::conditional_t<is_uniform, UniformGridBase<Theta>, NonUniformGridBase<Theta>>
{
};



// SPLINE DIMENSIONS -----------------------------------------------------------------------------
/**
 * @brief First spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct BSplinesR
    : std::conditional_t<
              is_uniform,
              ddc::UniformBSplines<R, BSplineDegree>,
              ddc::NonUniformBSplines<R, BSplineDegree>>
{
};

/**
 * @brief Second spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct BSplinesTheta
    : std::conditional_t<
              is_uniform,
              ddc::UniformBSplines<Theta, BSplineDegree>,
              ddc::NonUniformBSplines<Theta, BSplineDegree>>
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

} // namespace GEOM_NAMESPACE_NAME
