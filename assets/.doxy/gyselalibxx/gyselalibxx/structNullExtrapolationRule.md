

# Struct NullExtrapolationRule



[**ClassList**](annotated.md) **>** [**NullExtrapolationRule**](structNullExtrapolationRule.md)



_Define null extrapolation rule common to all geometries._ [More...](#detailed-description)

* `#include <null_extrapolation_rules.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION double | [**operator()**](#function-operator) (Coord&lt; Dim... &gt; const & coord\_extrap, [**MultipatchType**](classMultipatchType.md)&lt; SplinesOnPatch, Patches... &gt; const & patches\_splines, int const out\_of\_bounds\_idx) const<br>_Evaluate to zero the splines at a given coordinate outside of the splines domain._  |




























## Detailed Description




**See also:** [**MultipatchSplineEvaluator2D**](classMultipatchSplineEvaluator2D.md). 



    
## Public Functions Documentation




### function operator() 

_Evaluate to zero the splines at a given coordinate outside of the splines domain._ 
```C++
template<class... Dim, template< typename P > typename SplinesOnPatch, class... Patches>
inline KOKKOS_FUNCTION double NullExtrapolationRule::operator() (
    Coord< Dim... > const & coord_extrap,
    MultipatchType < SplinesOnPatch, Patches... > const & patches_splines,
    int const out_of_bounds_idx
) const
```





**Parameters:**


* `coord_extrap` Coordinate where we want to evaluate. 
* `patches_splines` Splines stored in a [**MultipatchType**](classMultipatchType.md). 
* `out_of_bounds_idx` Index of the localisation of the coordinate. It is supposed to be negative to be considered as outside of the domain. 



**Returns:**

A null double. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/null_extrapolation_rules.hpp`

