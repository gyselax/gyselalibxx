

# Struct ExtrapolationRule::Constant



[**ClassList**](annotated.md) **>** [**ExtrapolationRule**](namespaceExtrapolationRule.md) **>** [**Constant**](structExtrapolationRule_1_1Constant.md)



_Tag selecting constant extrapolation._ [More...](#detailed-description)

* `#include <i_interpolation.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::conditional\_t&lt; is\_spline\_basis\_v&lt; Basis &gt;, typename detail::DDCConstantExtrapolationRuleBuilder&lt; ddc::remove\_dims\_of\_t&lt; IdxRangeCoeff, find\_grid\_t&lt; typename Basis::continuous\_dimension\_type, ddc::to\_type\_seq\_t&lt; IdxRangeCoeff &gt; &gt; &gt;, typename Basis::continuous\_dimension\_type &gt;::type, [**ConstantIdentityInterpolationExtrapolationRule**](classConstantIdentityInterpolationExtrapolationRule.md)&lt; find\_grid\_t&lt; typename Basis::continuous\_dimension\_type, ddc::to\_type\_seq\_t&lt; IdxRangeCoeff &gt; &gt;, DataType &gt; &gt; | [**type**](#typedef-type)  <br>_The concrete extrapolation rule class for a given CoeffGrid/DataType._  |
















































## Detailed Description


The function is clamped to the value at the nearest boundary point. 


    
## Public Types Documentation




### typedef type 

_The concrete extrapolation rule class for a given CoeffGrid/DataType._ 
```C++
using ExtrapolationRule::Constant::type =  std::conditional_t< is_spline_basis_v<Basis>, typename detail::DDCConstantExtrapolationRuleBuilder< ddc::remove_dims_of_t< IdxRangeCoeff, find_grid_t< typename Basis::continuous_dimension_type, ddc::to_type_seq_t<IdxRangeCoeff> >>, typename Basis::continuous_dimension_type>::type, ConstantIdentityInterpolationExtrapolationRule< find_grid_t< typename Basis::continuous_dimension_type, ddc::to_type_seq_t<IdxRangeCoeff> >, DataType> >;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation.hpp`

