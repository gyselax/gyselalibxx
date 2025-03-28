

# Class MultipatchSplineBuilder2D

**template &lt;class ExecSpace, class MemorySpace, template&lt; typename P &gt; typename BSpline1OnPatch, template&lt; typename P &gt; typename BSpline2OnPatch, template&lt; typename P &gt; typename Grid1OnPatch, template&lt; typename P &gt; typename Grid2OnPatch, ddc::BoundCond BcLower1, ddc::BoundCond BcUpper1, ddc::BoundCond BcLower2, ddc::BoundCond BcUpper2, ddc::BoundCond BcTransition, class Connectivity, ddc::SplineSolver Solver, template&lt; typename P &gt; typename ValuesOnPatch, class... Patches&gt;**



[**ClassList**](annotated.md) **>** [**MultipatchSplineBuilder2D**](classMultipatchSplineBuilder2D.md)



_A class to call all the builders of all the patches once._ [More...](#detailed-description)

* `#include <multipatch_spline_builder_2d.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MultipatchSplineBuilder2D**](#function-multipatchsplinebuilder2d) (BuilderOnPatch&lt; Patches &gt; const &... builders) <br>_Instantiate the_ [_**MultipatchSplineBuilder**_](classMultipatchSplineBuilder.md) _from a std::tuple of all the builder on each patch._ |
|  void | [**operator()**](#function-operator) ([**MultipatchSplineCoeffs**](classMultipatchField.md) splines, [**MultipatchValues**](classMultipatchField.md) const & values, std::optional&lt; [**MultipatchDerivs1**](classMultipatchField.md) &gt; derivs\_min1=std::nullopt, std::optional&lt; [**MultipatchDerivs1**](classMultipatchField.md) &gt; derivs\_max1=std::nullopt, std::optional&lt; [**MultipatchDerivs2**](classMultipatchField.md) &gt; derivs\_min2=std::nullopt, std::optional&lt; [**MultipatchDerivs2**](classMultipatchField.md) &gt; derivs\_max2=std::nullopt, std::optional&lt; [**MultipatchDerivs12**](classMultipatchField.md) &gt; mixed\_derivs\_min1\_min2=std::nullopt, std::optional&lt; [**MultipatchDerivs12**](classMultipatchField.md) &gt; mixed\_derivs\_max1\_min2=std::nullopt, std::optional&lt; [**MultipatchDerivs12**](classMultipatchField.md) &gt; mixed\_derivs\_min1\_max2=std::nullopt, std::optional&lt; [**MultipatchDerivs12**](classMultipatchField.md) &gt; mixed\_derivs\_max1\_max2=std::nullopt) const<br>_Build the spline representation of each given function._  |
|   | [**~MultipatchSplineBuilder2D**](#function-multipatchsplinebuilder2d) () = default<br> |




























## Detailed Description


We need to instantiate all the builders for all the patches in the main code. We process the same way for the Field containing the spline coefficients and the values of the function on each patch. The Fields are stored in [**MultipatchField**](classMultipatchField.md) objects. The builders are stored in this class. This class is instantiated with all the builders. The operator() allows to call all the builders stored in the member of this class in one single line.


This function is useful to avoid calling all the builders individually, especially in a multipatch geometry with several patches.




**Template parameters:**


* `ExecSpace` The space (CPU/GPU) where the calculations are carried out. 
* `MemorySpace` The space (CPU/GPU) where the coefficients and values are stored. 
* `BSpline1OnPatch` A type alias which provides the first BSpline type along which the splines are built. 
* `BSpline2OnPatch` A type alias which provides the second BSpline type along which the splines are built. 
* `Grid1OnPatch` A type alias which provides the first Grid type along which the interpolation points of the splines are found. 
* `Grid2OnPatch` A type alias which provides the second Grid type along which the interpolation points of the splines are found. 
* `BcLower1` The lower boundary condition on the first dimension. 
* `BcUpper1` The upper boundary condition on the first dimension. 
* `BcLower2` The lower boundary condition on the second dimension. 
* `BcUpper2` The upper boundary condition on the second dimension. 
* `BcTransition` The boundary condition used at the interface between 2 patches. 
* `Connectivity` A [**MultipatchConnectivity**](classMultipatchConnectivity.md) object describing the interfaces between patches. 
* `Solver` The SplineSolver giving the backend used to perform the spline approximation. See DDC for more details. 
* `ValuesOnPatch` A type alias which provides the field type which will be used to pass the values of the function at the interpolation points. The index range of this field type should contain any batch dimensions. 




    
## Public Functions Documentation




### function MultipatchSplineBuilder2D 

_Instantiate the_ [_**MultipatchSplineBuilder**_](classMultipatchSplineBuilder.md) _from a std::tuple of all the builder on each patch._
```C++
inline explicit MultipatchSplineBuilder2D::MultipatchSplineBuilder2D (
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
inline void MultipatchSplineBuilder2D::operator() (
    MultipatchSplineCoeffs splines,
    MultipatchValues const & values,
    std::optional< MultipatchDerivs1 > derivs_min1=std::nullopt,
    std::optional< MultipatchDerivs1 > derivs_max1=std::nullopt,
    std::optional< MultipatchDerivs2 > derivs_min2=std::nullopt,
    std::optional< MultipatchDerivs2 > derivs_max2=std::nullopt,
    std::optional< MultipatchDerivs12 > mixed_derivs_min1_min2=std::nullopt,
    std::optional< MultipatchDerivs12 > mixed_derivs_max1_min2=std::nullopt,
    std::optional< MultipatchDerivs12 > mixed_derivs_min1_max2=std::nullopt,
    std::optional< MultipatchDerivs12 > mixed_derivs_max1_max2=std::nullopt
) const
```





**Parameters:**


* `splines` [**MultipatchField**](classMultipatchField.md) of all the Fields pointing to the spline representations. 
* `values` [**MultipatchField**](classMultipatchField.md) of all the Fields pointing to the function values. 
* `derivs_min1` [**MultipatchField**](classMultipatchField.md) of all the ConstFields describing the function derivatives in the first dimension at the lower bound of the second dimension. 
* `derivs_max1` [**MultipatchField**](classMultipatchField.md) of all the ConstFields describing the function derivatives in the first dimension at the upper bound of the second dimension. 
* `derivs_min2` [**MultipatchField**](classMultipatchField.md) of all the ConstFields describing the function derivatives in the second dimension at the lower bound of the first dimension. 
* `derivs_max2` [**MultipatchField**](classMultipatchField.md) of all the ConstFields describing the function derivatives in the second dimension at the upper bound of the first dimension. 
* `mixed_derivs_min1_min2` The values of the the cross-derivatives at the lower boundary in the first dimension and the lower boundary in the second dimension. 
* `mixed_derivs_max1_min2` The values of the the cross-derivatives at the upper boundary in the first dimension and the lower boundary in the second dimension. 
* `mixed_derivs_min1_max2` The values of the the cross-derivatives at the lower boundary in the first dimension and the upper boundary in the second dimension. 
* `mixed_derivs_max1_max2` The values of the the cross-derivatives at the upper boundary in the first dimension and the upper boundary in the second dimension. 




        

<hr>



### function ~MultipatchSplineBuilder2D 

```C++
MultipatchSplineBuilder2D::~MultipatchSplineBuilder2D () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/multipatch_spline_builder_2d.hpp`

