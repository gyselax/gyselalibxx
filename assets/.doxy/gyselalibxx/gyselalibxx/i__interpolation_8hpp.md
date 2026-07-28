

# File i\_interpolation.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**i\_interpolation.hpp**](i__interpolation_8hpp.md)

[Go to the source code of this file](i__interpolation_8hpp_source.md)



* `#include <ddc/kernels/splines.hpp>`
* `#include "constant_identity_interpolation_extrapolation_rule.hpp"`
* `#include "geometry_descriptors.hpp"`
* `#include "i_interpolation_builder.hpp"`
* `#include "i_interpolation_evaluator.hpp"`
* `#include "lagrange_basis_non_uniform.hpp"`
* `#include "lagrange_basis_uniform.hpp"`
* `#include "type_seq_tools.hpp"`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**ExtrapolationRule**](namespaceExtrapolationRule.md) <br>_A namespace containing tag types describing how a function is extrapolated outside the interpolation domain._  |
| namespace | [**concepts**](namespaceconcepts.md) <br> |


## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Constant**](structExtrapolationRule_1_1Constant.md) <br>_Tag selecting constant extrapolation._  |
| struct | [**NullValue**](structExtrapolationRule_1_1NullValue.md) <br>_Tag selecting null extrapolation: the function evaluates to zero outside the domain._  |
| struct | [**Periodic**](structExtrapolationRule_1_1Periodic.md) <br>_Tag selecting periodic extrapolation._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**is\_lagrange\_basis\_v**](#variable-is_lagrange_basis_v)   = `is\_uniform\_lagrange\_basis\_v&lt;Basis&gt; \|\| is\_non\_uniform\_lagrange\_basis\_v&lt;Basis&gt;`<br>_A type trait that is true when Basis is a Lagrange basis type._  |
|  constexpr bool | [**is\_spline\_basis\_v**](#variable-is_spline_basis_v)   = `ddc::is\_uniform\_bsplines\_v&lt;Basis&gt; \|\| ddc::is\_non\_uniform\_bsplines\_v&lt;Basis&gt;`<br>_A type trait that is true when Basis is a B-spline basis type._  |












































## Public Attributes Documentation




### variable is\_lagrange\_basis\_v 

_A type trait that is true when Basis is a Lagrange basis type._ 
```C++
constexpr bool is_lagrange_basis_v;
```



Evaluates to true for both [**UniformLagrangeBasis**](classUniformLagrangeBasis.md) and [**NonUniformLagrangeBasis**](classNonUniformLagrangeBasis.md) instantiations, false for all other types.




**Template parameters:**


* `Basis` The basis type to test. 




        

<hr>



### variable is\_spline\_basis\_v 

_A type trait that is true when Basis is a B-spline basis type._ 
```C++
constexpr bool is_spline_basis_v;
```



Evaluates to true for both ddc::UniformBSplines and ddc::NonUniformBSplines instantiations, false for all other types.




**Template parameters:**


* `Basis` The basis type to test. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation.hpp`

