

# Class InterpolationBuilder2DCollection&lt; ddc::TypeSeq&lt; Patches... &gt;, ddc::TypeSeq&lt; Builders... &gt; &gt;

**template &lt;concepts::Patch... Patches, concepts::InterpolationBuilder... Builders&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationBuilder2DCollection&lt; ddc::TypeSeq&lt; Patches... &gt;, ddc::TypeSeq&lt; Builders... &gt; &gt;**](classInterpolationBuilder2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Type8e3f4e34d093a5390cc74d679af413db.md)



_A class to call all the builders of all the patches once._ [More...](#detailed-description)

* `#include <interpolation_builder_2d_collection.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef MultipatchIdxRangeBatchedBasis&lt; MultipatchIdxRangeBatched &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_Multi-patch batched domain with the interpolation grids replaced by the basis grids._  |
| typedef MultipatchIdxRangeBatchedDeriv&lt; MultipatchIdxRangeBatched &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_Multi-patch batched domain with the interpolation grids replaced by deriv types._  |
| typedef MultipatchType&lt; IdxRangeCoeffOnPatch, Patches... &gt; | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The index range for the interpolation coefficients._  |
| typedef typename [**InterpolationBuilderTraits**](structInterpolationBuilderTraits.md)&lt; FirstBuilder &gt;::data\_type | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef typename FirstBuilder::exec\_space | [**exec\_space**](#typedef-exec_space)  <br>_Execution space._  |
| typedef MultipatchType&lt; IdxRangeInterpolationOnPatch, Patches... &gt; | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The ND index range for the interpolation mesh._  |
| typedef typename FirstBuilder::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_Memory space._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**InterpolationBuilder2DCollection**](#function-interpolationbuilder2dcollection) (Builders const &... builders) <br>_Instantiate the_ [_**InterpolationBuilder2DCollection**_](classInterpolationBuilder2DCollection.md) _from the builders on each patch._ |
|  MultipatchIdxRangeBatchedBasis&lt; MultipatchIdxRangeBatched &gt; const | [**batched\_basis\_idx\_ranges**](#function-batched_basis_idx_ranges) (MultipatchIdxRangeBatched idx\_ranges) const<br>_Get a MultipatchType collecting the index ranges typed on the whole basis dimensions. These index ranges are adapted to allocate interpolation coefficients on the multi-patch domain._  |
|  void | [**operator()**](#function-operator) (MultipatchCoeffs&lt; MultipatchValues &gt; coeffs, MultipatchValues const & values, std::optional&lt; MultipatchDerivs1&lt; MultipatchValues &gt; &gt; derivs\_min1=std::nullopt, std::optional&lt; MultipatchDerivs1&lt; MultipatchValues &gt; &gt; derivs\_max1=std::nullopt, std::optional&lt; MultipatchDerivs2&lt; MultipatchValues &gt; &gt; derivs\_min2=std::nullopt, std::optional&lt; MultipatchDerivs2&lt; MultipatchValues &gt; &gt; derivs\_max2=std::nullopt, std::optional&lt; MultipatchDerivs12&lt; MultipatchValues &gt; &gt; cross\_derivs\_min1\_min2=std::nullopt, std::optional&lt; MultipatchDerivs12&lt; MultipatchValues &gt; &gt; cross\_derivs\_max1\_min2=std::nullopt, std::optional&lt; MultipatchDerivs12&lt; MultipatchValues &gt; &gt; cross\_derivs\_min1\_max2=std::nullopt, std::optional&lt; MultipatchDerivs12&lt; MultipatchValues &gt; &gt; cross\_derivs\_max1\_max2=std::nullopt) const<br>_Build the interpolation coefficients of each given function._  |
|  void | [**operator()**](#function-operator_1) (MultipatchCoeffs&lt; MultipatchDerivField &gt; coeffs, MultipatchDerivField functions\_and\_derivs) const<br>_Build the interpolation coefficients of each given function._  |
|   | [**~InterpolationBuilder2DCollection**](#function-interpolationbuilder2dcollection) () = default<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions._  |


























## Detailed Description


See the documentation of the primary template.




**Template parameters:**


* `Patches` The patches. 
* `Builders` The 2D interpolation builders on each patch, in the same order as the patches. 




    
## Public Types Documentation




### typedef batched\_basis\_idx\_range\_type 

_Multi-patch batched domain with the interpolation grids replaced by the basis grids._ 
```C++
using InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::batched_basis_idx_range_type =  MultipatchIdxRangeBatchedBasis<MultipatchIdxRangeBatched>;
```




<hr>



### typedef batched\_derivs\_idx\_range\_type 

_Multi-patch batched domain with the interpolation grids replaced by deriv types._ 
```C++
using InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::batched_derivs_idx_range_type =  MultipatchIdxRangeBatchedDeriv<MultipatchIdxRangeBatched>;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The index range for the interpolation coefficients._ 
```C++
using InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::coeff_idx_range_type =  MultipatchType<IdxRangeCoeffOnPatch, Patches...>;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::data_type =  typename InterpolationBuilderTraits<FirstBuilder>::data_type;
```




<hr>



### typedef exec\_space 

_Execution space._ 
```C++
using InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::exec_space =  typename FirstBuilder::exec_space;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The ND index range for the interpolation mesh._ 
```C++
using InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::interpolation_idx_range_type =  MultipatchType<IdxRangeInterpolationOnPatch, Patches...>;
```




<hr>



### typedef memory\_space 

_Memory space._ 
```C++
using InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::memory_space =  typename FirstBuilder::memory_space;
```




<hr>
## Public Functions Documentation




### function InterpolationBuilder2DCollection 

_Instantiate the_ [_**InterpolationBuilder2DCollection**_](classInterpolationBuilder2DCollection.md) _from the builders on each patch._
```C++
inline explicit InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::InterpolationBuilder2DCollection (
    Builders const &... builders
) 
```





**Warning:**

The builders have to be sorted in the same order as the patches.




**Parameters:**


* `builders` Interpolation builders for each patch. 




        

<hr>



### function batched\_basis\_idx\_ranges 

_Get a MultipatchType collecting the index ranges typed on the whole basis dimensions. These index ranges are adapted to allocate interpolation coefficients on the multi-patch domain._ 
```C++
template<class MultipatchIdxRangeBatched>
inline MultipatchIdxRangeBatchedBasis< MultipatchIdxRangeBatched > const InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::batched_basis_idx_ranges (
    MultipatchIdxRangeBatched idx_ranges
) const
```





**Parameters:**


* `idx_ranges` MultipatchType collecting the index ranges on the grids of the patches. 



**Returns:**

a MultipatchType collecting the index ranges typed on the whole basis dimensions. 





        

<hr>



### function operator() 

_Build the interpolation coefficients of each given function._ 
```C++
template<class MultipatchValues>
inline void InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::operator() (
    MultipatchCoeffs< MultipatchValues > coeffs,
    MultipatchValues const & values,
    std::optional< MultipatchDerivs1< MultipatchValues > > derivs_min1=std::nullopt,
    std::optional< MultipatchDerivs1< MultipatchValues > > derivs_max1=std::nullopt,
    std::optional< MultipatchDerivs2< MultipatchValues > > derivs_min2=std::nullopt,
    std::optional< MultipatchDerivs2< MultipatchValues > > derivs_max2=std::nullopt,
    std::optional< MultipatchDerivs12< MultipatchValues > > cross_derivs_min1_min2=std::nullopt,
    std::optional< MultipatchDerivs12< MultipatchValues > > cross_derivs_max1_min2=std::nullopt,
    std::optional< MultipatchDerivs12< MultipatchValues > > cross_derivs_min1_max2=std::nullopt,
    std::optional< MultipatchDerivs12< MultipatchValues > > cross_derivs_max1_max2=std::nullopt
) const
```





**Parameters:**


* `coeffs` MultipatchField of all the Fields pointing to the interpolation coefficients. 
* `values` MultipatchField of all the Fields pointing to the function values. 
* `derivs_min1` MultipatchField of all the ConstFields describing the function derivatives in the first dimension at the lower bound of the second dimension. 
* `derivs_max1` MultipatchField of all the ConstFields describing the function derivatives in the first dimension at the upper bound of the second dimension. 
* `derivs_min2` MultipatchField of all the ConstFields describing the function derivatives in the second dimension at the lower bound of the first dimension. 
* `derivs_max2` MultipatchField of all the ConstFields describing the function derivatives in the second dimension at the upper bound of the first dimension. 
* `cross_derivs_min1_min2` The values of the the cross-derivatives at the lower boundary in the first dimension and the lower boundary in the second dimension. 
* `cross_derivs_max1_min2` The values of the the cross-derivatives at the upper boundary in the first dimension and the lower boundary in the second dimension. 
* `cross_derivs_min1_max2` The values of the the cross-derivatives at the lower boundary in the first dimension and the upper boundary in the second dimension. 
* `cross_derivs_max1_max2` The values of the the cross-derivatives at the upper boundary in the first dimension and the upper boundary in the second dimension.

The derivatives are only used by builders which are ddc::SplineBuilder2D. They are ignored by any other builder.




**Template parameters:**


* `MultipatchValues` The type of the MultipatchField containing the function values. The index range of each field may contain batch dimensions. 




        

<hr>



### function operator() 

_Build the interpolation coefficients of each given function._ 
```C++
inline void InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::operator() (
    MultipatchCoeffs< MultipatchDerivField > coeffs,
    MultipatchDerivField functions_and_derivs
) const
```



The derivatives are only used by builders which are ddc::SplineBuilder2D. Any other builder only uses the function values.




**Parameters:**


* `coeffs` MultipatchField of all the Fields pointing to the interpolation coefficients. 
* `functions_and_derivs` MultipatchField of all the DerivFields pointing to the function values and their derivatives. 




        

<hr>



### function ~InterpolationBuilder2DCollection 

```C++
InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::~InterpolationBuilder2DCollection () = default
```




<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions._ 
```C++
static inline std::size_t InterpolationBuilder2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interpolation/interpolation_builder_2d_collection.hpp`

