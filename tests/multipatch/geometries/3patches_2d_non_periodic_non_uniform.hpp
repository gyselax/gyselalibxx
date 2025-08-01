// SPDX-License-Identifier: MIT
#pragma once

#define GEOM_NAMESPACE_NAME non_periodic_non_uniform_2d_3patches
#include "3patches_2d_non_periodic.hpp"
#undef GEOM_NAMESPACE_NAME

// #pragma once

// #include <ddc/ddc.hpp>
// #include <ddc/kernels/splines.hpp>

// #include "connectivity.hpp"
// #include "edge.hpp"
// #include "interface.hpp"
// #include "patch.hpp"

// namespace non_periodic_non_uniform_2d_3patches {

// int constexpr BSplineDegree = 3;

// // CONTINUOUS DIMENSIONS -------------------------------------------------------------------------
// /**
//  * @brief First continuous dimension of patch PatchIdx.
//  */
// struct X
// {
//     /// @brief Non periodic dimension.
//     static bool constexpr PERIODIC = false;
//     /// A boolean indicating if dimension describes a covariant coordinate.
//     static bool constexpr IS_COVARIANT = true;
//     /// A boolean indicating if dimension describes a contravariant coordinate.
//     static bool constexpr IS_CONTRAVARIANT = true;
//     /// A type-alias mapping to the covariant counterpart.
//     using Dual = X;
// };

// /**
//  * @brief Second continuous dimension of patch PatchIdx.
//  */
// struct Y
// {
//     /// @brief Non periodic dimension.
//     static bool constexpr PERIODIC = false;
//     /// A boolean indicating if dimension describes a covariant coordinate.
//     static bool constexpr IS_COVARIANT = true;
//     /// A boolean indicating if dimension describes a contravariant coordinate.
//     static bool constexpr IS_CONTRAVARIANT = true;
//     /// A type-alias mapping to the covariant counterpart.
//     using Dual = Y;
// };



// // DISCRETE DIMENSIONS ---------------------------------------------------------------------------
// /**
//  * @brief Points sequence on the second logical dimension of patch 1.
//  * @tparam PatchIdx Index of the patch. 
//  */
// template <int PatchIdx>
// struct GridX : NonUniformGridBase<X>
// {
// };

// /**
//  * @brief Points sequence on the second logical dimension of patch 1.
//  * @tparam PatchIdx Index of the patch. 
//  */
// template <int PatchIdx>
// struct GridY : NonUniformGridBase<Y>
// {
// };



// // SPLINE DIMENSIONS -----------------------------------------------------------------------------
// /**
//  * @brief First spline dimension of patch PatchIdx.
//  * @tparam PatchIdx Index of the patch. 
//  */
// template <int PatchIdx>
// struct BSplinesX : ddc::NonUniformBSplines<X, BSplineDegree>
// {
// };

// /**
//  * @brief Second spline dimension of patch PatchIdx.
//  * @tparam PatchIdx Index of the patch. 
//  */
// template <int PatchIdx>
// struct BSplinesY : ddc::NonUniformBSplines<Y, BSplineDegree>
// {
// };


// // PATCHES ---------------------------------------------------------------------------------------
// /// @brief First patch.
// using Patch1 = Patch<GridX<1>, GridY<1>, BSplinesX<1>, BSplinesY<1>>;

// /// @brief Second patch.
// using Patch2 = Patch<GridX<2>, GridY<2>, BSplinesX<2>, BSplinesY<2>>;

// /// @brief Third patch.
// using Patch3 = Patch<GridX<3>, GridY<3>, BSplinesX<3>, BSplinesY<3>>;

// /// @brief Sorted list of patches.
// using PatchOrdering = ddc::detail::TypeSeq<Patch1, Patch2, Patch3>;


// // EDGES -----------------------------------------------------------------------------------------
// template <int PatchIdx>
// using NorthEdge
//         = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
//                GridY<PatchIdx>,
//                BACK>;
// template <int PatchIdx>
// using SouthEdge
//         = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
//                GridY<PatchIdx>,
//                FRONT>;
// template <int PatchIdx>
// using EastEdge
//         = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
//                GridX<PatchIdx>,
//                BACK>;
// template <int PatchIdx>
// using WestEdge
//         = Edge<typename ddc::type_seq_element_t<PatchIdx - 1, PatchOrdering>,
//                GridX<PatchIdx>,
//                FRONT>;



// } // namespace non_periodic_non_uniform_2d_3patches
