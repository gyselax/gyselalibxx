

# File spline\_interpolation.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**spline\_interpolation.hpp**](spline__interpolation_8hpp.md)

[Go to the source code of this file](spline__interpolation_8hpp_source.md)



* `#include <utility>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "extrapolation_rule_choice.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef detail::SplineInterpolator&lt; ExecSpace, Basis, InterpGrid, extrapolation\_rule\_t&lt; ExtrapRules, InterpGrid, double, Basis &gt;, MinBound, MaxBound, Solver &gt; | [**SplineInterpolator**](#typedef-splineinterpolator)  <br>_A helper alias to define an instance of detail::SplineInterpolator._  |
















































## Public Types Documentation




### typedef SplineInterpolator 

_A helper alias to define an instance of detail::SplineInterpolator._ 
```C++
using SplineInterpolator =  detail::SplineInterpolator< ExecSpace, Basis, InterpGrid, extrapolation_rule_t<ExtrapRules, InterpGrid, double, Basis>, MinBound, MaxBound, Solver>;
```



The helper allows ExtrapRules to be more general. It is a ddc::detail::TypeSeq&lt;MinExtrapolationRule, MaxExtrapolationRule&gt; pairing the extrapolation rules applied below/above the boundary. Each may be one of the tags in the [**ExtrapolationRule**](namespaceExtrapolationRule.md) namespace (e.g. [**ExtrapolationRule::Periodic**](namespaceExtrapolationRule.md#typedef-periodic)) or a custom, already-concrete extrapolation rule class. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/spline_interpolation.hpp`

