

# File interface\_derivative\_coefficients\_collection.hpp



[**FileList**](files.md) **>** [**interface\_derivatives**](dir_d1bd52a3e76a422151eefdcc4e15c189.md) **>** [**interface\_derivative\_coefficients\_collection.hpp**](interface__derivative__coefficients__collection_8hpp.md)

[Go to the source code of this file](interface__derivative__coefficients__collection_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "interface_derivative_coefficients.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**InterfaceDerivCoeffsCollection**](classInterfaceDerivCoeffsCollection.md) &lt;Interfaces&gt;<br>_A class to store a collection of interface derivative calculators templated on the interfaces._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_interface\_derivative\_coefficients\_collection**](#variable-enable_interface_derivative_coefficients_collection)   = `false`<br> |
|  constexpr bool | [**enable\_interface\_derivative\_coefficients\_collection&lt; InterfaceDerivCoeffsCollection&lt; DerivCalculatorType... &gt; &gt;**](#variable-enable_interface_derivative_coefficients_collection-interfacederivcoeffscollection-derivcalculatortype)   = `true`<br> |
|  constexpr bool | [**is\_single\_derivative\_calculator\_collection\_v**](#variable-is_single_derivative_calculator_collection_v)   = `/* multi line expression */`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**InterfaceDerivCoeffsCollection**](#function-interfacederivcoeffscollection) (DerivCalculatorType const &... derivative\_calculators) <br> |




























## Public Attributes Documentation




### variable enable\_interface\_derivative\_coefficients\_collection 

```C++
constexpr bool enable_interface_derivative_coefficients_collection;
```




<hr>



### variable enable\_interface\_derivative\_coefficients\_collection&lt; InterfaceDerivCoeffsCollection&lt; DerivCalculatorType... &gt; &gt; 

```C++
constexpr bool enable_interface_derivative_coefficients_collection< InterfaceDerivCoeffsCollection< DerivCalculatorType... > >;
```




<hr>



### variable is\_single\_derivative\_calculator\_collection\_v 

```C++
constexpr bool is_single_derivative_calculator_collection_v;
```




<hr>
## Public Functions Documentation




### function InterfaceDerivCoeffsCollection 

```C++
template<class... DerivCalculatorType>
InterfaceDerivCoeffsCollection (
    DerivCalculatorType const &... derivative_calculators
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/interface_derivative_coefficients_collection.hpp`

