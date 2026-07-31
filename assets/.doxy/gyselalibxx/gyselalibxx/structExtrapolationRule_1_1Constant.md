

# Struct ExtrapolationRule::Constant



[**ClassList**](annotated.md) **>** [**ExtrapolationRule**](namespaceExtrapolationRule.md) **>** [**Constant**](structExtrapolationRule_1_1Constant.md)



_Tag selecting constant extrapolation._ [More...](#detailed-description)

* `#include <i_interpolation.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::conditional\_t&lt; is\_spline\_basis\_v&lt; Basis &gt;, ddc::ConstantExtrapolationRule&lt; typename CoeffGrid::continuous\_dimension\_type &gt;, [**ConstantIdentityInterpolationExtrapolationRule**](classConstantIdentityInterpolationExtrapolationRule.md)&lt; CoeffGrid, DataType &gt; &gt; | [**type**](#typedef-type)  <br>_The concrete extrapolation rule class for a given CoeffGrid/DataType._  |
















































## Detailed Description


The function is clamped to the value at the nearest boundary point. 


    
## Public Types Documentation




### typedef type 

_The concrete extrapolation rule class for a given CoeffGrid/DataType._ 
```C++
using ExtrapolationRule::Constant::type =  std::conditional_t< is_spline_basis_v<Basis>, ddc::ConstantExtrapolationRule<typename CoeffGrid::continuous_dimension_type>, ConstantIdentityInterpolationExtrapolationRule<CoeffGrid, DataType> >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation.hpp`

