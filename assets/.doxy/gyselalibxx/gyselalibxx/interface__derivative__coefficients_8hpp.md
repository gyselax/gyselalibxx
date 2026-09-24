

# File interface\_derivative\_coefficients.hpp



[**FileList**](files.md) **>** [**interface\_derivatives**](dir_d1bd52a3e76a422151eefdcc4e15c189.md) **>** [**interface\_derivative\_coefficients.hpp**](interface__derivative__coefficients_8hpp.md)

[Go to the source code of this file](interface__derivative__coefficients_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "edge.hpp"`
* `#include "edge_transformation.hpp"`
* `#include "geometry_descriptors.hpp"`
* `#include "types.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md) &lt;class InterfaceType&gt;<br>_Compute the coefficients a, b and c of the interface derivative reconstruction method._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_interface\_derivative\_coefficients**](#variable-enable_interface_derivative_coefficients)   = `false`<br> |
|  constexpr bool | [**enable\_interface\_derivative\_coefficients&lt; InterfaceDerivCoeffs&lt; InterfaceType &gt; &gt;**](#variable-enable_interface_derivative_coefficients-interfacederivcoeffs-interfacetype)   = `true`<br> |
|  constexpr bool | [**is\_single\_derivative\_calculator\_v**](#variable-is_single_derivative_calculator_v)   = `enable\_interface\_derivative\_coefficients&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;T&gt;&gt;&gt;`<br> |












































## Public Attributes Documentation




### variable enable\_interface\_derivative\_coefficients 

```C++
constexpr bool enable_interface_derivative_coefficients;
```




<hr>



### variable enable\_interface\_derivative\_coefficients&lt; InterfaceDerivCoeffs&lt; InterfaceType &gt; &gt; 

```C++
constexpr bool enable_interface_derivative_coefficients< InterfaceDerivCoeffs< InterfaceType > >;
```




<hr>



### variable is\_single\_derivative\_calculator\_v 

```C++
constexpr bool is_single_derivative_calculator_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/interface_derivative_coefficients.hpp`

