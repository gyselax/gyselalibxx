

# File lagrange\_interpolation.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**lagrange\_interpolation.hpp**](lagrange__interpolation_8hpp.md)

[Go to the source code of this file](lagrange__interpolation_8hpp_source.md)



* `#include <utility>`
* `#include "extrapolation_rule_choice.hpp"`
* `#include "identity_interpolation_builder.hpp"`
* `#include "lagrange_basis_non_uniform.hpp"`
* `#include "lagrange_basis_uniform.hpp"`
* `#include "lagrange_evaluator.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef detail::LagrangeInterpolator&lt; ExecSpace, Basis, InterpGrid, extrapolation\_rule\_t&lt; ExtrapRules, typename [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md)&lt; ExecSpace, typename ExecSpace::memory\_space, DataType, InterpGrid, Basis &gt;::basis\_domain\_type, DataType, Basis &gt;, DataType &gt; | [**LagrangeInterpolator**](#typedef-lagrangeinterpolator)  <br>_A helper alias to define an instance of detail::LagrangeInterpolator._  |
















































## Public Types Documentation




### typedef LagrangeInterpolator 

_A helper alias to define an instance of detail::LagrangeInterpolator._ 
```C++
using LagrangeInterpolator =  detail::LagrangeInterpolator< ExecSpace, Basis, InterpGrid, extrapolation_rule_t< ExtrapRules, typename IdentityInterpolationBuilder< ExecSpace, typename ExecSpace::memory_space, DataType, InterpGrid, Basis>::basis_domain_type, DataType, Basis>, DataType>;
```



The helper allows ExtrapRules to be more general. It is a ddc::detail::TypeSeq&lt;MinExtrapolationRule, MaxExtrapolationRule&gt; pairing the extrapolation rules applied below/above the boundary. Each may be one of the tags in the [**ExtrapolationRule**](namespaceExtrapolationRule.md) namespace (e.g. [**ExtrapolationRule::Periodic**](namespaceExtrapolationRule.md#typedef-periodic)) or a custom, already-concrete extrapolation rule class. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/lagrange_interpolation.hpp`

