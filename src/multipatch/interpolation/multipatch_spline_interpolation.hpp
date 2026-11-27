// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_interpolation_2d.hpp"
#include "spline_interpolation.hpp"

namespace detail {

/**
 * @brief A helper to define the 2D spline interpolation used on each patch of a
 * MultipatchSplineInterpolator2D.
 *
 * @tparam ExecSpace The Kokkos execution space used for computations.
 * @tparam Solver The spline solver backend.
 * @tparam BoundaryClosures1 The spline boundary closures along the first dimension of each patch.
 * @tparam BoundaryClosures2 The spline boundary closures along the second dimension of each patch.
 */
template <
        class ExecSpace,
        ddc::SplineSolver Solver,
        class BoundaryClosures1,
        class BoundaryClosures2>
struct SplineInterpolatorOnPatch
{
    /// @brief The 2D spline interpolation on the given patch.
    template <class Patch>
    using type = ::SplineInterpolatorWSolver<
            ExecSpace,
            typename Patch::IdxRangeBS12,
            typename Patch::IdxRange12,
            Solver,
            PatchExtrapolationRules<typename Patch::Dim1>,
            PatchExtrapolationRules<typename Patch::Dim2>,
            BoundaryClosures1,
            BoundaryClosures2>;
};

} // namespace detail

/**
 * @brief A MultipatchInterpolation2D which uses a 2D spline interpolation on every patch.
 *
 * On each patch, the spline interpolation is defined on the B-splines (BSplines1, BSplines2)
 * and the grids (Grid1, Grid2) of the patch. The same boundary closures are used on every
 * patch. They are given as tags (e.g. SplineBoundaryClosure::Greville_Greville) which are
 * resolved separately on each patch.
 *
 * The extrapolation outside of the global domain is handled by MultipatchExtrapolationRule,
 * so the extrapolation rules on each patch are null (or periodic along periodic dimensions).
 *
 * The execution space is the one of the patch locator (PatchLocator::exec_space).
 *
 * @tparam PatchLocator The operator used to find the patch where a coordinate is physically located.
 * @tparam MultipatchExtrapolationRule The extrapolation rule used outside of the global domain
 *          (e.g. ConstantExtrapolationRuleOnion or NullExtrapolationRule).
 * @tparam BoundaryClosures1 A SplineBoundaryClosures describing the spline boundary closures
 *          along the first dimension of each patch.
 * @tparam BoundaryClosures2 A SplineBoundaryClosures describing the spline boundary closures
 *          along the second dimension of each patch.
 * @tparam Solver The spline solver backend (default: LAPACK).
 */
template <
        class PatchLocator,
        class MultipatchExtrapolationRule,
        class BoundaryClosures1,
        class BoundaryClosures2,
        ddc::SplineSolver Solver = ddc::SplineSolver::LAPACK>
using MultipatchSplineInterpolator2D = MultipatchInterpolation2D<
        detail::SplineInterpolatorOnPatch<
                typename PatchLocator::exec_space,
                Solver,
                BoundaryClosures1,
                BoundaryClosures2>::template type,
        PatchLocator,
        MultipatchExtrapolationRule>;
