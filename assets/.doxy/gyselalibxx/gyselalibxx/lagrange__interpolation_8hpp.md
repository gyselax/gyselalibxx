

# File lagrange\_interpolation.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**lagrange\_interpolation.hpp**](lagrange__interpolation_8hpp.md)

[Go to the source code of this file](lagrange__interpolation_8hpp_source.md)



* `#include <utility>`
* `#include "extrapolation_rule_choice.hpp"`
* `#include "identity_interpolation_builder.hpp"`
* `#include "lagrange_basis_non_uniform.hpp"`
* `#include "lagrange_basis_uniform.hpp"`
* `#include "lagrange_evaluator.hpp"`
* `#include "nd_identity_interpolation_builder.hpp"`
* `#include "nd_lagrange_evaluator.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename detail::LagrangeInterpolatorResolver&lt; ExecSpace, DataType, IdxRangeBasis, IdxRangeInterpGrid, ExtrapRules... &gt;::type | [**LagrangeInterpolator**](#typedef-lagrangeinterpolator)  <br>_A helper alias to define an instance of detail::LagrangeInterpolator._  |
















































## Public Types Documentation




### typedef LagrangeInterpolator 

_A helper alias to define an instance of detail::LagrangeInterpolator._ 
```C++
using LagrangeInterpolator =  typename detail::LagrangeInterpolatorResolver< ExecSpace, DataType, IdxRangeBasis, IdxRangeInterpGrid, ExtrapRules...>::type;
```



The helper allows ExtrapRules to be more general. It is a ddc::detail::TypeSeq&lt;MinExtrapolationRule, MaxExtrapolationRule&gt; pairing the extrapolation rules applied below/above the boundary. Each may be one of the tags in the [**ExtrapolationRule**](namespaceExtrapolationRule.md) namespace (e.g. [**ExtrapolationRule::Periodic**](namespaceExtrapolationRule.md#typedef-periodic)) or a custom, already-concrete extrapolation rule class. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/lagrange_interpolation.hpp`

