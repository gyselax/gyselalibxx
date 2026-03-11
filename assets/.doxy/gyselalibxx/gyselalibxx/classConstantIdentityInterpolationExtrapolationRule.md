

# Class ConstantIdentityInterpolationExtrapolationRule

**template &lt;class CoeffGrid, class DataType&gt;**



[**ClassList**](annotated.md) **>** [**ConstantIdentityInterpolationExtrapolationRule**](classConstantIdentityInterpolationExtrapolationRule.md)



_A constant extrapolation rule for identity-based (Lagrange) interpolation._ [More...](#detailed-description)

* `#include <constant_identity_interpolation_extrapolation_rule.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ConstantIdentityInterpolationExtrapolationRule**](#function-constantidentityinterpolationextrapolationrule) (Idx&lt; CoeffGrid &gt; coeff\_idx) <br>_Construct the rule, clamping to the coefficient at_ `coeff_idx` _._ |
|  KOKKOS\_FUNCTION double | [**operator()**](#function-operator) (CoordType pos, ConstField&lt; DataType, IdxRange&lt; CoeffGrid &gt;, Layout, MemorySpace &gt; const interp\_coef) const<br>_Return the boundary interpolation coefficient for a coordinate outside the domain._  |




























## Detailed Description


When an evaluation coordinate falls outside the interpolation domain, this rule returns the interpolation coefficient at a fixed boundary index (`m_coeff_idx`) rather than extrapolating the polynomial. This is the Lagrange analog of ddc::ConstantExtrapolationRule used with spline evaluators.


It is intended to be constructed via get\_extrapolation&lt;CONSTANT, CoeffGrid, Basis&gt;(), which supplies the appropriate boundary index automatically.




**Template parameters:**


* `CoeffGrid` The discrete grid type of the interpolation coefficients. 
* `DataType` The floating-point type of the function values (default: double). 




    
## Public Functions Documentation




### function ConstantIdentityInterpolationExtrapolationRule 

_Construct the rule, clamping to the coefficient at_ `coeff_idx` _._
```C++
inline explicit ConstantIdentityInterpolationExtrapolationRule::ConstantIdentityInterpolationExtrapolationRule (
    Idx< CoeffGrid > coeff_idx
) 
```





**Parameters:**


* `coeff_idx` The index of the boundary coefficient to return for all out-of-domain evaluation coordinates. 




        

<hr>



### function operator() 

_Return the boundary interpolation coefficient for a coordinate outside the domain._ 
```C++
template<class CoordType, class Layout, class MemorySpace>
inline KOKKOS_FUNCTION double ConstantIdentityInterpolationExtrapolationRule::operator() (
    CoordType pos,
    ConstField< DataType, IdxRange< CoeffGrid >, Layout, MemorySpace > const interp_coef
) const
```





**Parameters:**


* `pos` The coordinate where the function is to be evaluated (unused). 
* `interp_coef` The full field of interpolation coefficients.



**Returns:**

The value of `interp_coef` at the fixed boundary index. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/constant_identity_interpolation_extrapolation_rule.hpp`

