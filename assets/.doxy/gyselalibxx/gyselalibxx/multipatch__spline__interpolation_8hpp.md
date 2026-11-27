

# File multipatch\_spline\_interpolation.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_54175eb71bf3bdb1ca864f2a0d352472.md) **>** [**multipatch\_spline\_interpolation.hpp**](multipatch__spline__interpolation_8hpp.md)

[Go to the source code of this file](multipatch__spline__interpolation_8hpp_source.md)



* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "multipatch_interpolation_2d.hpp"`
* `#include "spline_interpolation.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef MultipatchInterpolation2D&lt; detail::SplineInterpolatorOnPatch&lt; typename PatchLocator::exec\_space, Solver, BoundaryClosures1, BoundaryClosures2 &gt;::template type, PatchLocator, MultipatchExtrapolationRule &gt; | [**MultipatchSplineInterpolator2D**](#typedef-multipatchsplineinterpolator2d)  <br>_A MultipatchInterpolation2D which uses a 2D spline interpolation on every patch._  |
















































## Public Types Documentation




### typedef MultipatchSplineInterpolator2D 

_A MultipatchInterpolation2D which uses a 2D spline interpolation on every patch._ 
```C++
using MultipatchSplineInterpolator2D =  MultipatchInterpolation2D<
        detail::SplineInterpolatorOnPatch<
                typename PatchLocator::exec_space,
                Solver,
                BoundaryClosures1,
                BoundaryClosures2>::template type,
        PatchLocator,
        MultipatchExtrapolationRule>;
```



On each patch, the spline interpolation is defined on the B-splines (BSplines1, BSplines2) and the grids (Grid1, Grid2) of the patch. The same boundary closures are used on every patch. They are given as tags (e.g. [**SplineBoundaryClosure::Greville\_Greville**](namespaceSplineBoundaryClosure.md#typedef-greville_greville)) which are resolved separately on each patch.


The extrapolation outside of the global domain is handled by MultipatchExtrapolationRule, so the extrapolation rules on each patch are null (or periodic along periodic dimensions).


The execution space is the one of the patch locator (PatchLocator::exec\_space).




**Template parameters:**


* `PatchLocator` The operator used to find the patch where a coordinate is physically located. 
* `MultipatchExtrapolationRule` The extrapolation rule used outside of the global domain (e.g. [**ConstantExtrapolationRuleOnion**](structConstantExtrapolationRuleOnion.md) or [**NullExtrapolationRule**](structNullExtrapolationRule.md)). 
* `BoundaryClosures1` A [**SplineBoundaryClosures**](structSplineBoundaryClosures.md) describing the spline boundary closures along the first dimension of each patch. 
* `BoundaryClosures2` A [**SplineBoundaryClosures**](structSplineBoundaryClosures.md) describing the spline boundary closures along the second dimension of each patch. 
* `Solver` The spline solver backend (default: LAPACK). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interpolation/multipatch_spline_interpolation.hpp`

