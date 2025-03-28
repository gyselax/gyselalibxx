

# Class MultipatchSplineBuilder

**template &lt;class ExecSpace, class MemorySpace, template&lt; typename P &gt; typename BSplineOnPatch, template&lt; typename P &gt; typename GridOnPatch, ddc::BoundCond BcLower, ddc::BoundCond BcUpper, ddc::BoundCond BcTransition, class Connectivity, ddc::SplineSolver Solver, template&lt; typename P &gt; typename ValuesOnPatch, class... Patches&gt;**



[**ClassList**](annotated.md) **>** [**MultipatchSplineBuilder**](classMultipatchSplineBuilder.md)



_A class to call all the builders of all the patches once._ [More...](#detailed-description)

* `#include <multipatch_spline_builder.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MultipatchSplineBuilder**](#function-multipatchsplinebuilder) (BuilderOnPatch&lt; Patches &gt; const &... builders) <br>_Instantiate the_ [_**MultipatchSplineBuilder**_](classMultipatchSplineBuilder.md) _from a std::tuple of all the builder on each patch._ |
|  void | [**operator()**](#function-operator) ([**MultipatchSplineCoeffs**](classMultipatchField.md) splines, [**MultipatchValues**](classMultipatchField.md) const & values, std::optional&lt; [**MultipatchDerivs**](classMultipatchField.md) &gt; derivs\_xmin=std::nullopt, std::optional&lt; [**MultipatchDerivs**](classMultipatchField.md) &gt; derivs\_xmax=std::nullopt) const<br>_Build the spline representation of each given function._  |
|   | [**~MultipatchSplineBuilder**](#function-multipatchsplinebuilder) () = default<br> |




























## Detailed Description


We need to instantiate all the builders for all the patches in the main code. We process the same way for the Fields containing the spline coefficients and the values of the function on each patch. The Fields are stored in [**MultipatchField**](classMultipatchField.md) objects. The builders are stored in this class. This class is instantiated with all the builders. The operator() allows to call all the builders stored in the member of this class in one single line.


This function is useful to avoid calling all the builders individually, especially in a multipatch geometry with several patches.




**Template parameters:**


* `ExecSpace` The space (CPU/GPU) where the calculations are carried out. 
* `MemorySpace` The space (CPU/GPU) where the coefficients and values are stored. 
* `BSplineOnPatch` A type alias which provides the BSpline type along which the splines are built. 
* `GridOnPatch` A type alias which provides the Grid type along which the interpolation points of the splines are found. 
* `BcLower` The lower boundary condition. 
* `BcUpper` The upper boundary condition. 
* `BcTransition` The boundary condition used at the interface between 2 patches. 
* `Connectivity` A [**MultipatchConnectivity**](classMultipatchConnectivity.md) object describing the interfaces between patches. 
* `Solver` The SplineSolver giving the backend used to perform the spline approximation. See DDC for more details. 
* `ValuesOnPatch` A type alias which provides the field type which will be used to pass the values of the function at the interpolation points. The index range of this field type should contain any batch dimensions. 




    
## Public Functions Documentation




### function MultipatchSplineBuilder 

_Instantiate the_ [_**MultipatchSplineBuilder**_](classMultipatchSplineBuilder.md) _from a std::tuple of all the builder on each patch._
```C++
inline explicit MultipatchSplineBuilder::MultipatchSplineBuilder (
    BuilderOnPatch< Patches > const &... builders
) 
```





**Warning:**

The builders have to be sorted in the same order as the patches in the tuple.




**Parameters:**


* `builders` Spline builders for each patch. 




        

<hr>



### function operator() 

_Build the spline representation of each given function._ 
```C++
inline void MultipatchSplineBuilder::operator() (
    MultipatchSplineCoeffs splines,
    MultipatchValues const & values,
    std::optional< MultipatchDerivs > derivs_xmin=std::nullopt,
    std::optional< MultipatchDerivs > derivs_xmax=std::nullopt
) const
```





**Parameters:**


* `splines` [**MultipatchField**](classMultipatchField.md) of all the Fields pointing to the spline representations. 
* `values` [**MultipatchField**](classMultipatchField.md) of all the Fields pointing to the function values. 
* `derivs_xmin` [**MultipatchField**](classMultipatchField.md) of all the ConstFields describing the function derivatives at the lower bound. 
* `derivs_xmax` [**MultipatchField**](classMultipatchField.md) of all the ConstFields describing the function derivatives at the upper bound. 




        

<hr>



### function ~MultipatchSplineBuilder 

```C++
MultipatchSplineBuilder::~MultipatchSplineBuilder () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/multipatch_spline_builder.hpp`

