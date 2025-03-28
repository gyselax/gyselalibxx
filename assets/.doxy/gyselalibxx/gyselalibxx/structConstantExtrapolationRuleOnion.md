

# Struct ConstantExtrapolationRuleOnion

**template &lt;class PatchLocator&gt;**



[**ClassList**](annotated.md) **>** [**ConstantExtrapolationRuleOnion**](structConstantExtrapolationRuleOnion.md)



_Define constant extrapolation rule for onion shape geometries. Struct useful for the MultipatchSplineEvaluator types._  __[More...](#detailed-description)

* `#include <constant_extrapolation_rules_onion.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ConstantExtrapolationRuleOnion**](#function-constantextrapolationruleonion) (Coord&lt; R\_min &gt; const & r\_min, Coord&lt; R\_max &gt; const & r\_max) <br>_Instantiate a_ [_**ConstantExtrapolationRuleOnion**_](structConstantExtrapolationRuleOnion.md) _. The R1 and R2 templates are needed for GPU._ |
|  KOKKOS\_FUNCTION double | [**operator()**](#function-operator) (Coord&lt; Dim... &gt; const & coord\_extrap, [**MultipatchField**](classMultipatchField.md)&lt; SplinesOnPatch, Patches... &gt; const & patches\_splines, int const out\_of\_bounds\_idx) const<br>_Evaluate at a given outside coordinate._  |




























## Detailed Description




**Template parameters:**


* `PatchLocator` A patch locator specialised for onion shape geometries. 




    
## Public Functions Documentation




### function ConstantExtrapolationRuleOnion 

_Instantiate a_ [_**ConstantExtrapolationRuleOnion**_](structConstantExtrapolationRuleOnion.md) _. The R1 and R2 templates are needed for GPU._
```C++
inline explicit ConstantExtrapolationRuleOnion::ConstantExtrapolationRuleOnion (
    Coord< R_min > const & r_min,
    Coord< R_max > const & r_max
) 
```





**Template parameters:**


* `R1` Continuous dimension of the minimum radial coordinate of the domain. 
* `R2` Continuous dimension of the maximum radial coordinate of the domain. 



**Parameters:**


* `r_min` Minimum radial coordinate of the domain. 
* `r_max` Maximum radial coordinate of the domain. 




        

<hr>



### function operator() 

_Evaluate at a given outside coordinate._ 
```C++
template<class... Dim, template< typename P > typename SplinesOnPatch, class... Patches>
inline KOKKOS_FUNCTION double ConstantExtrapolationRuleOnion::operator() (
    Coord< Dim... > const & coord_extrap,
    MultipatchField < SplinesOnPatch, Patches... > const & patches_splines,
    int const out_of_bounds_idx
) const
```





**Template parameters:**


* `Dim` Continuous dimensions where the given coordinate is defined. 
* `SplinesOnPatch` Field of spline coefficients template on the [**Patch**](structPatch.md). 
* `Patches` [**Patch**](structPatch.md) types. 



**Parameters:**


* `coord_extrap` Coordinate where we want to evaluate. 
* `patches_splines` Splines stored in a [**MultipatchType**](classMultipatchType.md). 
* `out_of_bounds_idx` Index of the localisation of the coordinate. It is supposed to be negative to be considered as outside of the domain. 



**Returns:**

A double with the evaluation outside of the domain. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/constant_extrapolation_rules_onion.hpp`

