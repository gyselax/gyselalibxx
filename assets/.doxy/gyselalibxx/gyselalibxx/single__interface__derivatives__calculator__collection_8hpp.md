

# File single\_interface\_derivatives\_calculator\_collection.hpp



[**FileList**](files.md) **>** [**interface\_derivatives**](dir_d1bd52a3e76a422151eefdcc4e15c189.md) **>** [**single\_interface\_derivatives\_calculator\_collection.hpp**](single__interface__derivatives__calculator__collection_8hpp.md)

[Go to the source code of this file](single__interface__derivatives__calculator__collection_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "single_interface_derivatives_calculator.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**SingleInterfaceDerivativesCalculatorCollection**](classSingleInterfaceDerivativesCalculatorCollection.md) &lt;Interfaces&gt;<br>_A class to store a collection of interface derivative calculators templated on the interfaces._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**enable\_single\_derivative\_calculator\_collection**](#variable-enable_single_derivative_calculator_collection)   = `false`<br> |
|  constexpr bool | [**enable\_single\_derivative\_calculator\_collection&lt; SingleInterfaceDerivativesCalculatorCollection&lt; DerivCalculatorType... &gt; &gt;**](#variable-enable_single_derivative_calculator_collection-singleinterfacederivativescalculatorcollection-derivcalculatortype)   = `true`<br> |
|  constexpr bool | [**is\_single\_derivative\_calculator\_collection\_v**](#variable-is_single_derivative_calculator_collection_v)   = `/* multi line expression */`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SingleInterfaceDerivativesCalculatorCollection**](#function-singleinterfacederivativescalculatorcollection) (DerivCalculatorType const &... derivative\_calculators) <br> |




























## Public Attributes Documentation




### variable enable\_single\_derivative\_calculator\_collection 

```C++
constexpr bool enable_single_derivative_calculator_collection;
```




<hr>



### variable enable\_single\_derivative\_calculator\_collection&lt; SingleInterfaceDerivativesCalculatorCollection&lt; DerivCalculatorType... &gt; &gt; 

```C++
constexpr bool enable_single_derivative_calculator_collection< SingleInterfaceDerivativesCalculatorCollection< DerivCalculatorType... > >;
```




<hr>



### variable is\_single\_derivative\_calculator\_collection\_v 

```C++
constexpr bool is_single_derivative_calculator_collection_v;
```




<hr>
## Public Functions Documentation




### function SingleInterfaceDerivativesCalculatorCollection 

```C++
template<class... DerivCalculatorType>
SingleInterfaceDerivativesCalculatorCollection (
    DerivCalculatorType const &... derivative_calculators
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/single_interface_derivatives_calculator_collection.hpp`

