

# File spline\_interpolation.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**spline\_interpolation.hpp**](spline__interpolation_8hpp.md)

[Go to the source code of this file](spline__interpolation_8hpp_source.md)



* `#include <utility>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "extrapolation_rule_choice.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**SplineBoundaryClosure**](namespaceSplineBoundaryClosure.md) <br>_Predefined_ [_**SplineBoundaryClosures**_](structSplineBoundaryClosures.md) _for the common case where the same closure applies at both boundaries._ |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**SplineBoundaryClosures**](structSplineBoundaryClosures.md) &lt;MinClosure, MaxClosure&gt;<br>_Groups the lower (min) and upper (max) ddc::SplineBuilderClosure of a spline builder into a single non-type template argument._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename detail::SplineInterpolatorResolver&lt; ExecSpace, IdxRangeBasis, IdxRangeInterpGrid, ddc::SplineSolver::LAPACK, TailTags... &gt;::type | [**SplineInterpolator**](#typedef-splineinterpolator)  <br>_Resolves to the 1D or 2D SplineInterpolator matching the given index ranges, using the default (LAPACK) spline solver backend._  |
| typedef typename detail::SplineInterpolatorResolver&lt; ExecSpace, IdxRangeBasis, IdxRangeInterpGrid, Solver, TailTags... &gt;::type | [**SplineInterpolatorWSolver**](#typedef-splineinterpolatorwsolver)  <br>_Resolves to the 1D or 2D SplineInterpolator matching the given index ranges, using an explicitly chosen spline solver backend._  |
















































## Public Types Documentation




### typedef SplineInterpolator 

_Resolves to the 1D or 2D SplineInterpolator matching the given index ranges, using the default (LAPACK) spline solver backend._ 
```C++
using SplineInterpolator =  typename detail::SplineInterpolatorResolver< ExecSpace, IdxRangeBasis, IdxRangeInterpGrid, ddc::SplineSolver::LAPACK, TailTags...>::type;
```



The resulting type is chosen by detail::SplineInterpolatorResolver based on the dimensionality of `IdxRangeBasis` and `IdxRangeInterpGrid:` detail::SplineInterpolator for the 1D case, detail::SplineInterpolator2D for the 2D case.




**Template parameters:**


* `ExecSpace` The Kokkos execution space used for computations. 
* `IdxRangeBasis` The index range of the B-spline basis (basis(es) for 2D). 
* `IdxRangeInterpGrid` The index range of the interpolation grid(s). 
* `TailTags` The extrapolation rule(s) and boundary closure(s) describing the boundary conditions, as expected by detail::SplineInterpolatorResolver. 




        

<hr>



### typedef SplineInterpolatorWSolver 

_Resolves to the 1D or 2D SplineInterpolator matching the given index ranges, using an explicitly chosen spline solver backend._ 
```C++
using SplineInterpolatorWSolver =  typename detail::SplineInterpolatorResolver< ExecSpace, IdxRangeBasis, IdxRangeInterpGrid, Solver, TailTags...>::type;
```



The resulting type is chosen by detail::SplineInterpolatorResolver based on the dimensionality of `IdxRangeBasis` and `IdxRangeInterpGrid:` detail::SplineInterpolator for the 1D case, detail::SplineInterpolator2D for the 2D case.




**Template parameters:**


* `ExecSpace` The Kokkos execution space used for computations. 
* `IdxRangeBasis` The index range of the B-spline basis (basis(es) for 2D). 
* `IdxRangeInterpGrid` The index range of the interpolation grid(s). 
* `Solver` The spline solver backend to use. 
* `TailTags` The extrapolation rule(s) and boundary closure(s) describing the boundary conditions, as expected by detail::SplineInterpolatorResolver. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/spline_interpolation.hpp`

