

# Class MultipatchSplineBuilder2D

**template &lt;class ExecSpace, class MemorySpace, template&lt; typename P &gt; typename BSpline1OnPatch, template&lt; typename P &gt; typename BSpline2OnPatch, template&lt; typename P &gt; typename Grid1OnPatch, template&lt; typename P &gt; typename Grid2OnPatch, ddc::SplineBuilderClosure SBCLower1, ddc::SplineBuilderClosure SBCUpper1, ddc::SplineBuilderClosure SBCLower2, ddc::SplineBuilderClosure SBCUpper2, ddc::SplineBuilderClosure BcTransition, class Connectivity, ddc::SplineSolver Solver, template&lt; typename P &gt; typename ValuesOnPatch, class... Patches&gt;**



[**ClassList**](annotated.md) **>** [**MultipatchSplineBuilder2D**](classMultipatchSplineBuilder2D.md)



_A class to call all the builders of all the patches once._ [More...](#detailed-description)

* `#include <multipatch_spline_builder_2d.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef MultipatchIdxRangeBatchedBS&lt; MultipatchIdxRangeBatched &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_Multi-patch batched domain with InterpolationDDim replaced by BSplines._  |
| typedef MultipatchIdxRangeBatchedDeriv&lt; MultipatchIdxRangeBatched &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_Multi-patch batched domain with InterpolationDDim replaced by deriv\_type._  |
| typedef MultipatchIdxRangeBS | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The index range for the interpolation coefficients._  |
| typedef double | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_Execution space._  |
| typedef MultipatchIdxRange | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The ND index range for the interpolation mesh._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_Memory space._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MultipatchSplineBuilder2D**](#function-multipatchsplinebuilder2d) (BuilderOnPatch&lt; Patches &gt; const &... builders) <br>_Instantiate the_ [_**MultipatchSplineBuilder**_](classMultipatchSplineBuilder.md) _from a std::tuple of all the builder on each patch._ |
|  MultipatchIdxRangeBatchedBS&lt; MultipatchIdxRangeBatched &gt; const | [**batched\_basis\_idx\_ranges**](#function-batched_basis_idx_ranges) (MultipatchIdxRangeBatched idx\_ranges) const<br>_Get a MultipatchType collecting the index ranges typed on the whole B-splines dimensions. These index ranges are adapted to allocate splines on the multi-patch domain._  |
|  void | [**operator()**](#function-operator) (MultipatchSplineCoeffs splines, MultipatchValues const & values, std::optional&lt; MultipatchDerivs1 &gt; derivs\_min1=std::nullopt, std::optional&lt; MultipatchDerivs1 &gt; derivs\_max1=std::nullopt, std::optional&lt; MultipatchDerivs2 &gt; derivs\_min2=std::nullopt, std::optional&lt; MultipatchDerivs2 &gt; derivs\_max2=std::nullopt, std::optional&lt; MultipatchDerivs12 &gt; cross\_derivs\_min1\_min2=std::nullopt, std::optional&lt; MultipatchDerivs12 &gt; cross\_derivs\_max1\_min2=std::nullopt, std::optional&lt; MultipatchDerivs12 &gt; cross\_derivs\_min1\_max2=std::nullopt, std::optional&lt; MultipatchDerivs12 &gt; cross\_derivs\_max1\_max2=std::nullopt) const<br>_Build the spline representation of each given function._  |
|  void | [**operator()**](#function-operator_1) (MultipatchSplineCoeffs splines, MultipatchDerivField functions\_and\_derivs) const<br>_Build the spline representation of each given function._  |
|  MultipatchIdxRangeBS const | [**spline\_idx\_ranges**](#function-spline_idx_ranges) () const<br>_Get a MultipatchType collecting the index ranges typed on the B-splines dimensions of the patches. These index ranges are adapted to allocate splines on the multi-patch domain._  |
|   | [**~MultipatchSplineBuilder2D**](#function-multipatchsplinebuilder2d) () = default<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions._  |


























## Detailed Description


We need to instantiate all the builders for all the patches in the main code. We process the same way for the Field containing the spline coefficients and the values of the function on each patch. The Fields are stored in MultipatchField objects. The builders are stored in this class. This class is instantiated with all the builders. The operator() allows to call all the builders stored in the member of this class in one single line.


This function is useful to avoid calling all the builders individually, especially in a multipatch geometry with several patches.




**Warning:**

[**MultipatchSplineBuilder2D**](classMultipatchSplineBuilder2D.md) is implemented for tensor-product multi-patch decompositions.




**Template parameters:**


* `ExecSpace` The space (CPU/GPU) where the calculations are carried out. 
* `MemorySpace` The space (CPU/GPU) where the coefficients and values are stored. 
* `BSpline1OnPatch` A type alias which provides the first BSpline type along which the splines are built. 
* `BSpline2OnPatch` A type alias which provides the second BSpline type along which the splines are built. 
* `Grid1OnPatch` A type alias which provides the first Grid type along which the interpolation points of the splines are found. 
* `Grid2OnPatch` A type alias which provides the second Grid type along which the interpolation points of the splines are found. 
* `SBCLower1` The lower spline closure on the first dimension. 
* `SBCUpper1` The upper spline closure on the first dimension. 
* `SBCLower2` The lower spline closure on the second dimension. 
* `SBCUpper2` The upper spline closure on the second dimension. 
* `BcTransition` The boundary condition used at the interface between 2 patches. 
* `Connectivity` A [**MultipatchConnectivity**](classMultipatchConnectivity.md) object describing the interfaces between patches. 
* `Solver` The SplineSolver giving the backend used to perform the spline approximation. See DDC for more details. 
* `ValuesOnPatch` A type alias which provides the field type which will be used to pass the values of the function at the interpolation points. The index range of this field type should contain any batch dimensions. 




    
## Public Types Documentation




### typedef batched\_basis\_idx\_range\_type 

_Multi-patch batched domain with InterpolationDDim replaced by BSplines._ 
```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::batched_basis_idx_range_type =  MultipatchIdxRangeBatchedBS<MultipatchIdxRangeBatched>;
```




<hr>



### typedef batched\_derivs\_idx\_range\_type 

_Multi-patch batched domain with InterpolationDDim replaced by deriv\_type._ 
```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::batched_derivs_idx_range_type =  MultipatchIdxRangeBatchedDeriv<MultipatchIdxRangeBatched>;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The index range for the interpolation coefficients._ 
```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::coeff_idx_range_type =  MultipatchIdxRangeBS;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::data_type =  double;
```




<hr>



### typedef exec\_space 

_Execution space._ 
```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::exec_space =  ExecSpace;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The ND index range for the interpolation mesh._ 
```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::interpolation_idx_range_type =  MultipatchIdxRange;
```




<hr>



### typedef memory\_space 

_Memory space._ 
```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::memory_space =  MemorySpace;
```




<hr>
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



### function batched\_basis\_idx\_ranges 

_Get a MultipatchType collecting the index ranges typed on the whole B-splines dimensions. These index ranges are adapted to allocate splines on the multi-patch domain._ 
```C++
template<class MultipatchIdxRangeBatched>
inline MultipatchIdxRangeBatchedBS< MultipatchIdxRangeBatched > const MultipatchSplineBuilder2D::batched_basis_idx_ranges (
    MultipatchIdxRangeBatched idx_ranges
) const
```





**Parameters:**


* `idx_ranges` MultipatchType collecting the index ranges on the grids of the patches. 



**Returns:**

a MultipatchType collecting the index ranges typed on the whole B-splines dimensions. 





        

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
    std::optional< MultipatchDerivs12 > cross_derivs_min1_min2=std::nullopt,
    std::optional< MultipatchDerivs12 > cross_derivs_max1_min2=std::nullopt,
    std::optional< MultipatchDerivs12 > cross_derivs_min1_max2=std::nullopt,
    std::optional< MultipatchDerivs12 > cross_derivs_max1_max2=std::nullopt
) const
```





**Parameters:**


* `splines` MultipatchField of all the Fields pointing to the spline representations. 
* `values` MultipatchField of all the Fields pointing to the function values. 
* `derivs_min1` MultipatchField of all the ConstFields describing the function derivatives in the first dimension at the lower bound of the second dimension. 
* `derivs_max1` MultipatchField of all the ConstFields describing the function derivatives in the first dimension at the upper bound of the second dimension. 
* `derivs_min2` MultipatchField of all the ConstFields describing the function derivatives in the second dimension at the lower bound of the first dimension. 
* `derivs_max2` MultipatchField of all the ConstFields describing the function derivatives in the second dimension at the upper bound of the first dimension. 
* `cross_derivs_min1_min2` The values of the the cross-derivatives at the lower boundary in the first dimension and the lower boundary in the second dimension. 
* `cross_derivs_max1_min2` The values of the the cross-derivatives at the upper boundary in the first dimension and the lower boundary in the second dimension. 
* `cross_derivs_min1_max2` The values of the the cross-derivatives at the lower boundary in the first dimension and the upper boundary in the second dimension. 
* `cross_derivs_max1_max2` The values of the the cross-derivatives at the upper boundary in the first dimension and the upper boundary in the second dimension. 




        

<hr>



### function operator() 

_Build the spline representation of each given function._ 
```C++
inline void MultipatchSplineBuilder2D::operator() (
    MultipatchSplineCoeffs splines,
    MultipatchDerivField functions_and_derivs
) const
```





**Parameters:**


* `splines` MultipatchField of all the Fields pointing to the spline representations. 
* `functions_and_derivs` MultipatchField of all the DerivFields pointing to the function values and their derivatives. 




        

<hr>



### function spline\_idx\_ranges 

_Get a MultipatchType collecting the index ranges typed on the B-splines dimensions of the patches. These index ranges are adapted to allocate splines on the multi-patch domain._ 
```C++
inline MultipatchIdxRangeBS const MultipatchSplineBuilder2D::spline_idx_ranges () const
```





**Returns:**

a MultipatchType collecting the index ranges typed on the B-splines dimensions of the patches. 





        

<hr>



### function ~MultipatchSplineBuilder2D 

```C++
MultipatchSplineBuilder2D::~MultipatchSplineBuilder2D () = default
```




<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions._ 
```C++
static inline std::size_t MultipatchSplineBuilder2D::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/multipatch_spline_builder_2d.hpp`

