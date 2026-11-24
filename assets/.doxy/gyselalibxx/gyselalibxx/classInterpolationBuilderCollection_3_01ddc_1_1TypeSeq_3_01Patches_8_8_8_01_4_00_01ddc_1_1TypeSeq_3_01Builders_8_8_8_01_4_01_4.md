

# Class InterpolationBuilderCollection&lt; ddc::TypeSeq&lt; Patches... &gt;, ddc::TypeSeq&lt; Builders... &gt; &gt;

**template &lt;concepts::Patch... Patches, concepts::InterpolationBuilder1D... Builders&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationBuilderCollection&lt; ddc::TypeSeq&lt; Patches... &gt;, ddc::TypeSeq&lt; Builders... &gt; &gt;**](classInterpolationBuilderCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1TypeSeq_3_01Builders_8_8_8_01_4_01_4.md)



_A class to call all the builders of all the patches once._ [More...](#detailed-description)

* `#include <interpolation_builder_collection.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef MultipatchIdxRangeBatchedBasis&lt; MultipatchIdxRangeBatched &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_Multi-patch batched domain with the interpolation grids replaced by the basis grids._  |
| typedef MultipatchIdxRangeBatchedDeriv&lt; MultipatchIdxRangeBatched &gt; | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_Multi-patch batched domain with the interpolation grids replaced by deriv types._  |
| typedef MultipatchType&lt; IdxRangeCoeffOnPatch, Patches... &gt; | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The index range for the interpolation coefficients._  |
| typedef typename [**InterpolationBuilderTraits**](structInterpolationBuilderTraits.md)&lt; FirstBuilder &gt;::data\_type | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef typename FirstBuilder::exec\_space | [**exec\_space**](#typedef-exec_space)  <br>_Execution space._  |
| typedef MultipatchType&lt; IdxRangeInterpolationOnPatch, Patches... &gt; | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The 1D index range for the interpolation mesh._  |
| typedef typename FirstBuilder::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_Memory space._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**InterpolationBuilderCollection**](#function-interpolationbuildercollection) (Builders const &... builders) <br>_Instantiate the_ [_**InterpolationBuilderCollection**_](classInterpolationBuilderCollection.md) _from the builders on each patch._ |
|  MultipatchIdxRangeBatchedBasis&lt; MultipatchIdxRangeBatched &gt; const | [**batched\_basis\_idx\_ranges**](#function-batched_basis_idx_ranges) (MultipatchIdxRangeBatched idx\_ranges) const<br>_Get a MultipatchType collecting the index ranges typed on the whole basis dimensions. These index ranges are adapted to allocate interpolation coefficients on the multi-patch domain._  |
|  void | [**operator()**](#function-operator) (MultipatchCoeffs&lt; MultipatchValues &gt; coeffs, MultipatchValues const & values, std::optional&lt; MultipatchDerivs&lt; MultipatchValues &gt; &gt; derivs\_xmin=std::nullopt, std::optional&lt; MultipatchDerivs&lt; MultipatchValues &gt; &gt; derivs\_xmax=std::nullopt) const<br>_Build the interpolation coefficients of each given function._  |
|   | [**~InterpolationBuilderCollection**](#function-interpolationbuildercollection) () = default<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions._  |


























## Detailed Description


See the documentation of the primary template.




**Template parameters:**


* `Patches` The patches. 
* `Builders` The 1D interpolation builders on each patch, in the same order as the patches. 




    
## Public Types Documentation




### typedef batched\_basis\_idx\_range\_type 

_Multi-patch batched domain with the interpolation grids replaced by the basis grids._ 
```C++
using InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::batched_basis_idx_range_type =  MultipatchIdxRangeBatchedBasis<MultipatchIdxRangeBatched>;
```




<hr>



### typedef batched\_derivs\_idx\_range\_type 

_Multi-patch batched domain with the interpolation grids replaced by deriv types._ 
```C++
using InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::batched_derivs_idx_range_type =  MultipatchIdxRangeBatchedDeriv<MultipatchIdxRangeBatched>;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The index range for the interpolation coefficients._ 
```C++
using InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::coeff_idx_range_type =  MultipatchType<IdxRangeCoeffOnPatch, Patches...>;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::data_type =  typename InterpolationBuilderTraits<FirstBuilder>::data_type;
```




<hr>



### typedef exec\_space 

_Execution space._ 
```C++
using InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::exec_space =  typename FirstBuilder::exec_space;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The 1D index range for the interpolation mesh._ 
```C++
using InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::interpolation_idx_range_type =  MultipatchType<IdxRangeInterpolationOnPatch, Patches...>;
```




<hr>



### typedef memory\_space 

_Memory space._ 
```C++
using InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::memory_space =  typename FirstBuilder::memory_space;
```




<hr>
## Public Functions Documentation




### function InterpolationBuilderCollection 

_Instantiate the_ [_**InterpolationBuilderCollection**_](classInterpolationBuilderCollection.md) _from the builders on each patch._
```C++
inline explicit InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::InterpolationBuilderCollection (
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
inline MultipatchIdxRangeBatchedBasis< MultipatchIdxRangeBatched > const InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::batched_basis_idx_ranges (
    MultipatchIdxRangeBatched idx_ranges
) const
```





**Parameters:**


* `idx_ranges` MultipatchType collecting the (batched) index ranges on the grids of the patches. 



**Returns:**

a MultipatchType collecting the index ranges typed on the whole basis dimensions. 





        

<hr>



### function operator() 

_Build the interpolation coefficients of each given function._ 
```C++
template<class MultipatchValues>
inline void InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::operator() (
    MultipatchCoeffs< MultipatchValues > coeffs,
    MultipatchValues const & values,
    std::optional< MultipatchDerivs< MultipatchValues > > derivs_xmin=std::nullopt,
    std::optional< MultipatchDerivs< MultipatchValues > > derivs_xmax=std::nullopt
) const
```





**Parameters:**


* `coeffs` MultipatchField of all the Fields pointing to the interpolation coefficients. 
* `values` MultipatchField of all the Fields pointing to the function values. 
* `derivs_xmin` MultipatchField of all the ConstFields describing the function derivatives at the lower bound. 
* `derivs_xmax` MultipatchField of all the ConstFields describing the function derivatives at the upper bound.



**Template parameters:**


* `MultipatchValues` The type of the MultipatchField containing the function values. The index range of each field may contain batch dimensions. 




        

<hr>



### function ~InterpolationBuilderCollection 

```C++
InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::~InterpolationBuilderCollection () = default
```




<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions._ 
```C++
static inline std::size_t InterpolationBuilderCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Builders... > >::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interpolation/interpolation_builder_collection.hpp`

