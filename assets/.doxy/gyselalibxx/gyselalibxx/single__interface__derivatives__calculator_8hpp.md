

# File single\_interface\_derivatives\_calculator.hpp



[**FileList**](files.md) **>** [**interface\_derivatives**](dir_d1bd52a3e76a422151eefdcc4e15c189.md) **>** [**single\_interface\_derivatives\_calculator.hpp**](single__interface__derivatives__calculator_8hpp.md)

[Go to the source code of this file](single__interface__derivatives__calculator_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "edge.hpp"`
* `#include "edge_transformation.hpp"`
* `#include "geometry_descriptors.hpp"`
* `#include "types.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**SingleInterfaceDerivativesCalculator**](classSingleInterfaceDerivativesCalculator.md) &lt;class InterfaceType&gt;<br>_Compute the derivative of an equivalent global spline at the interface between two patches._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_single\_derivative\_calculator**](#variable-enable_single_derivative_calculator)   = `false`<br> |
|  constexpr bool | [**enable\_single\_derivative\_calculator&lt; SingleInterfaceDerivativesCalculator&lt; InterfaceType &gt; &gt;**](#variable-enable_single_derivative_calculator-singleinterfacederivativescalculator-interfacetype)   = `true`<br> |
|  constexpr bool | [**is\_single\_derivative\_calculator\_v**](#variable-is_single_derivative_calculator_v)   = `enable\_single\_derivative\_calculator&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;T&gt;&gt;&gt;`<br> |












































## Public Attributes Documentation




### variable enable\_single\_derivative\_calculator 

```C++
constexpr bool enable_single_derivative_calculator;
```




<hr>



### variable enable\_single\_derivative\_calculator&lt; SingleInterfaceDerivativesCalculator&lt; InterfaceType &gt; &gt; 

```C++
constexpr bool enable_single_derivative_calculator< SingleInterfaceDerivativesCalculator< InterfaceType > >;
```




<hr>



### variable is\_single\_derivative\_calculator\_v 

```C++
constexpr bool is_single_derivative_calculator_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/single_interface_derivatives_calculator.hpp`

