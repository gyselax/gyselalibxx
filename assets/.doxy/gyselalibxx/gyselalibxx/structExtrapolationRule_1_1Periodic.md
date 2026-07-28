

# Struct ExtrapolationRule::Periodic



[**ClassList**](annotated.md) **>** [**ExtrapolationRule**](namespaceExtrapolationRule.md) **>** [**Periodic**](structExtrapolationRule_1_1Periodic.md)



_Tag selecting periodic extrapolation._ [More...](#detailed-description)

* `#include <i_interpolation.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::PeriodicExtrapolationRule&lt; typename CoeffGrid::continuous\_dimension\_type &gt; | [**type**](#typedef-type)  <br>_The concrete extrapolation rule class for a given CoeffGrid/DataType._  |
















































## Detailed Description


The value at a point outside the domain is taken as the value at the equivalent point inside the domain. 


    
## Public Types Documentation




### typedef type 

_The concrete extrapolation rule class for a given CoeffGrid/DataType._ 
```C++
using ExtrapolationRule::Periodic::type =  ddc::PeriodicExtrapolationRule<typename CoeffGrid::continuous_dimension_type>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation.hpp`

