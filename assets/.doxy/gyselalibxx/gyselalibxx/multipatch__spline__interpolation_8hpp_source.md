

# File multipatch\_spline\_interpolation.hpp

[**File List**](files.md) **>** [**interpolation**](dir_54175eb71bf3bdb1ca864f2a0d352472.md) **>** [**multipatch\_spline\_interpolation.hpp**](multipatch__spline__interpolation_8hpp.md)

[Go to the documentation of this file](multipatch__spline__interpolation_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once
#include <ddc/kernels/splines.hpp>

#include "ddc_aliases.hpp"
#include "multipatch_interpolation_2d.hpp"
#include "spline_interpolation.hpp"

namespace detail {

template <
        class ExecSpace,
        ddc::SplineSolver Solver,
        class BoundaryClosures1,
        class BoundaryClosures2>
struct SplineInterpolatorOnPatch
{
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
```


