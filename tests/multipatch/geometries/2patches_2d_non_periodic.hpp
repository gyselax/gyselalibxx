// SPDX-License-Identifier: MIT

/*
    Geometry defined here: 2 * 2D patches.
        - patch 1 and patch 2: 
            - non-periodic on X and Y;
            - uniform/non uniform cubic splines;
            - uniform/non uniform Grid1 and Grid2; 
*/
#pragma once
#include <type_traits>

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "connectivity.hpp"
#include "ddc_aliases.hpp"
#include "edge.hpp"
#include "interface.hpp"
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
struct GridX
    : std::conditional_t<is_uniform, UniformGridBase<X<PatchIdx>>, NonUniformGridBase<X<PatchIdx>>>
{
};

/**
 * @brief Points sequence on the second logical dimension of patch 1.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct GridY
    : std::conditional_t<is_uniform, UniformGridBase<Y<PatchIdx>>, NonUniformGridBase<Y<PatchIdx>>>
{
};



// SPLINE DIMENSIONS -----------------------------------------------------------------------------
/**
 * @brief First spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct BSplinesX
    : std::conditional_t<
              is_uniform,
              ddc::UniformBSplines<X<PatchIdx>, BSplineDegree>,
              ddc::NonUniformBSplines<X<PatchIdx>, BSplineDegree>>
{
};

/**
 * @brief Second spline dimension of patch PatchIdx.
 * @tparam PatchIdx Index of the patch. 
 */
template <int PatchIdx>
struct BSplinesY
    : std::conditional_t<
              is_uniform,
              ddc::UniformBSplines<Y<PatchIdx>, BSplineDegree>,
              ddc::NonUniformBSplines<Y<PatchIdx>, BSplineDegree>>
{
};


// PATCHES ---------------------------------------------------------------------------------------
/// @brief First patch.
using Patch1 = Patch<GridX<1>, GridY<1>, BSplinesX<1>, BSplinesY<1>>;

/// @brief Second patch.
using Patch2 = Patch<GridX<2>, GridY<2>, BSplinesX<2>, BSplinesY<2>>;

#ifndef EXTEND_GEOMETRY
/// @brief Sorted list of patches.
using PatchOrdering = ddc::detail::TypeSeq<Patch1, Patch2>;


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
using SouthInterface1 = Interface<SouthEdge<1>, OutsideEdge, true>;
using WestInterface1 = Interface<WestEdge<1>, OutsideEdge, true>;

using NorthInterface2 = Interface<NorthEdge<2>, OutsideEdge, true>;
using SouthInterface2 = Interface<SouthEdge<2>, OutsideEdge, true>;
using EastInterface2 = Interface<EastEdge<2>, OutsideEdge, true>;

using Interface_1_2 = Interface<EastEdge<1>, WestEdge<2>, true>;


// CONNECTIVITY ----------------------------------------------------------------------------------
using Connectivity = MultipatchConnectivity<
        NorthInterface1,
        SouthInterface1,
        WestInterface1,
        NorthInterface2,
        SouthInterface2,
        EastInterface2,
        Interface_1_2>;
#endif

} // namespace GEOM_NAMESPACE_NAME
