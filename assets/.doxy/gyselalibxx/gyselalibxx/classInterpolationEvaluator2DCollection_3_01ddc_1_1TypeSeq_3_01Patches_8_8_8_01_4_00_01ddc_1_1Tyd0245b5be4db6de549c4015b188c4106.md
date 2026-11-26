

# Class InterpolationEvaluator2DCollection&lt; ddc::TypeSeq&lt; Patches... &gt;, ddc::TypeSeq&lt; Evaluators... &gt;, PatchLocator, ExtrapolationRule &gt;

**template &lt;concepts::Patch... Patches, concepts::InterpolationEvaluator... Evaluators, class PatchLocator, class [**ExtrapolationRule**](namespaceExtrapolationRule.md)&gt;**



[**ClassList**](annotated.md) **>** [**InterpolationEvaluator2DCollection&lt; ddc::TypeSeq&lt; Patches... &gt;, ddc::TypeSeq&lt; Evaluators... &gt;, PatchLocator, ExtrapolationRule &gt;**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md)



_A class to evaluate all the interpolations of all the patches at once._ [More...](#detailed-description)

* `#include <interpolation_evaluator_2d_collection.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ConstField&lt; [**data\_type**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-data_type), [**coeff\_idx\_range\_type**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-coeff_idx_range_type)&lt; [**Patch**](structPatch.md) &gt;, [**memory\_space**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-memory_space) &gt; | [**CoeffOnPatch**](#typedef-coeffonpatch)  <br>_Type for MultipatchType: A field of 2D interpolation coefficients for a non-batched interpolation defined on both of the_ [_**Patch**_](structPatch.md) _'s logical dimensions. Needed public for functions on GPU._ |
| typedef ddc:: type\_seq\_element\_t&lt; 0, ddc::to\_type\_seq\_t&lt; typename Traits&lt; [**Patch**](structPatch.md) &gt;::basis\_idx\_range\_type &gt; &gt; | [**basis\_type1**](#typedef-basis_type1)  <br>_The discrete dimension representing the basis along first dimension._  |
| typedef ddc:: type\_seq\_element\_t&lt; 1, ddc::to\_type\_seq\_t&lt; typename Traits&lt; [**Patch**](structPatch.md) &gt;::basis\_idx\_range\_type &gt; &gt; | [**basis\_type2**](#typedef-basis_type2)  <br>_The discrete dimension representing the basis along second dimension._  |
| typedef typename ddc::detail::convert\_type\_seq\_to\_discrete\_domain\_t&lt; [**physical\_grids**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-physical_grids)&lt; ValuesField &gt; &gt; | [**batched\_evaluation\_idx\_range\_type**](#typedef-batched_evaluation_idx_range_type)  <br>_The type of the whole domain representing evaluation points._  |
| typedef typename Traits&lt; [**Patch**](structPatch.md) &gt;::coeff\_idx\_range\_type | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The type of the 2D index range on which the interpolation coefficients are defined._  |
| typedef typename [**evaluation\_discrete\_dimension\_type1**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-evaluation_discrete_dimension_type1)&lt; [**Patch**](structPatch.md) &gt;::continuous\_dimension\_type | [**continuous\_dimension\_type1**](#typedef-continuous_dimension_type1)  <br>_The type of the first evaluation continuous dimension used by this class._  |
| typedef typename [**evaluation\_discrete\_dimension\_type2**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-evaluation_discrete_dimension_type2)&lt; [**Patch**](structPatch.md) &gt;::continuous\_dimension\_type | [**continuous\_dimension\_type2**](#typedef-continuous_dimension_type2)  <br>_The type of the second evaluation continuous dimension used by this class._  |
| typedef typename [**InterpolationEvaluatorTraits**](structInterpolationEvaluatorTraits.md)&lt; FirstEvaluator &gt;::coord\_type | [**coord\_type**](#typedef-coord_type)  <br>_The ND coordinate type corresponding to the evaluation mesh._  |
| typedef typename [**InterpolationEvaluatorTraits**](structInterpolationEvaluatorTraits.md)&lt; FirstEvaluator &gt;::data\_type | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef detail::deriv\_sub\_set\_t&lt; ddc::to\_type\_seq\_t&lt; [**idx\_range\_type**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-idx_range_type)&lt; ValuesField &gt; &gt; &gt; | [**deriv\_tags**](#typedef-deriv_tags)  <br>_TypeSeq of the derivative dimensions (contained in_ [_**DerivField**_](classDerivField.md) _)._ |
| typedef ddc::type\_seq\_element\_t&lt; 0, ddc::to\_type\_seq\_t&lt; [**evaluation\_idx\_range\_type**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-evaluation_idx_range_type)&lt; [**Patch**](structPatch.md) &gt; &gt; &gt; | [**evaluation\_discrete\_dimension\_type1**](#typedef-evaluation_discrete_dimension_type1)  <br>_The type of the first discrete dimension of interest used by this class._  |
| typedef ddc::type\_seq\_element\_t&lt; 1, ddc::to\_type\_seq\_t&lt; [**evaluation\_idx\_range\_type**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-evaluation_idx_range_type)&lt; [**Patch**](structPatch.md) &gt; &gt; &gt; | [**evaluation\_discrete\_dimension\_type2**](#typedef-evaluation_discrete_dimension_type2)  <br>_The type of the second discrete dimension of interest used by this class._  |
| typedef typename Traits&lt; [**Patch**](structPatch.md) &gt;::evaluation\_idx\_range\_type | [**evaluation\_idx\_range\_type**](#typedef-evaluation_idx_range_type)  <br>_The type of the 2D evaluation mesh used by this class._  |
| typedef IdxRange&lt; [**evaluation\_discrete\_dimension\_type1**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-evaluation_discrete_dimension_type1)&lt; [**Patch**](structPatch.md) &gt; &gt; | [**evaluation\_idx\_range\_type1**](#typedef-evaluation_idx_range_type1)  <br>_The type of the domain for the 1D evaluation mesh along first dimension used by this class._  |
| typedef IdxRange&lt; [**evaluation\_discrete\_dimension\_type2**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-evaluation_discrete_dimension_type2)&lt; [**Patch**](structPatch.md) &gt; &gt; | [**evaluation\_idx\_range\_type2**](#typedef-evaluation_idx_range_type2)  <br>_The type of the domain for the 1D evaluation mesh along second dimension used by this class._  |
| typedef typename FirstEvaluator::exec\_space | [**exec\_space**](#typedef-exec_space)  <br>_The type of the Kokkos execution space used by this class._  |
| typedef typename ValuesField::discrete\_domain\_type | [**idx\_range\_type**](#typedef-idx_range_type)  <br>_The type of the whole discrete domain (dimensions for grids and derivatives)._  |
| typedef typename FirstEvaluator::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The type of the Kokkos memory space used by this class._  |
| typedef ddc::type\_seq\_remove\_t&lt; ddc::to\_type\_seq\_t&lt; [**idx\_range\_type**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-idx_range_type)&lt; ValuesField &gt; &gt;, [**deriv\_tags**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-deriv_tags)&lt; ValuesField &gt; &gt; | [**physical\_grids**](#typedef-physical_grids)  <br>_TypeSeq of the grid dimensions (without derivative dimensions)._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**n\_patches**](#variable-n_patches)   = `ddc::type\_seq\_size\_v&lt;PatchOrdering&gt;`<br>_The number of patches._  |














## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**InterpolationEvaluator2DCollection**](#function-interpolationevaluator2dcollection-12) (Evaluators const &... evaluators, PatchLocator const & patch\_locator, [**ExtrapolationRule**](namespaceExtrapolationRule.md) const & extrapolation\_rule) <br>_Instantiate a_ [_**InterpolationEvaluator2DCollection**_](classInterpolationEvaluator2DCollection.md) _._ |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**InterpolationEvaluator2DCollection**](#function-interpolationevaluator2dcollection-22) ([**InterpolationEvaluator2DCollection**](classInterpolationEvaluator2DCollection.md) const &) = default<br>_Default copy operator callable from GPU._  |
|  void | [**apply\_evaluator\_deriv\_field**](#function-apply_evaluator_deriv_field) (Idx&lt; DerivDims... &gt; const & deriv\_order, ValuesField const & patch\_values, CoordConstFieldOnPatch&lt; StoringPatch &gt; const & patch\_coords, MultipatchCoeffs const patches\_coeffs) const<br>_Compute the derivatives of a given patch at the coordinates defined on the given patch. Specialisation for patch\_values on_ [_**DerivField**_](classDerivField.md) _. Needed public for functions on GPU._ |
|  void | [**apply\_evaluator\_field**](#function-apply_evaluator_field) (Idx&lt; DerivDims... &gt; const & deriv\_order, ValuesField const & patch\_values, CoordConstFieldOnPatch&lt; StoringPatch &gt; const & patch\_coords, MultipatchCoeffs const patches\_coeffs) const<br>_Compute the values or the derivatives of a given patch at the coordinates defined on the given patch. Specialisation for patch\_values on DField. Needed public for functions on GPU._  |
|  void | [**apply\_integrate**](#function-apply_integrate) (double & integral, [**CoeffOnPatch**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-coeffonpatch)&lt; [**Patch**](structPatch.md) &gt; const & coeffs) const<br>_Integrate the spline defined on the given patch._  |
|  KOKKOS\_FUNCTION [**data\_type**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-data_type) | [**deriv**](#function-deriv-12) (Idx&lt; DerivDims... &gt; const & deriv\_order, Coord&lt; CDim1, CDim2 &gt; const & coord\_eval, MultipatchCoeffs const & patches\_coeffs) const<br>_Differentiate 2D interpolations (described by their interpolation coefficients) at a given coordinate along a specified dimension of interest._  |
|  void | [**deriv**](#function-deriv-22) (Idx&lt; DerivDims... &gt; const & deriv\_order, MultipatchValues const & patches\_deriv, MultipatchCoordField const & patches\_coords, MultipatchCoeffs const & patches\_coeffs) const<br>_Differentiate 2D interpolations (described by their interpolation coefficients) on meshes along first dimension of interest._  |
|  void | [**integrate**](#function-integrate) (DKokkosView\_h&lt; [**n\_patches**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#variable-n_patches) &gt; const & integrals, MultipatchCoeffs const & patches\_coeffs) const<br>_Integration of splines (described by their spline coefficients)._  |
|  KOKKOS\_FUNCTION [**data\_type**](classInterpolationEvaluator2DCollection_3_01ddc_1_1TypeSeq_3_01Patches_8_8_8_01_4_00_01ddc_1_1Tyd0245b5be4db6de549c4015b188c4106.md#typedef-data_type) | [**operator()**](#function-operator) (Coord const coord\_eval, MultipatchCoeffs const & patches\_coeffs) const<br>_Evaluate the 2D interpolations (described by their interpolation coefficients) at a given coordinate._  |
|  void | [**operator()**](#function-operator_1) (MultipatchValues const & patches\_values, MultipatchCoordField const & patches\_coords, MultipatchCoeffs const & patches\_coeffs) const<br>_Evaluate 2D interpolations (described by their interpolation coefficients) on meshes._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~InterpolationEvaluator2DCollection**](#function-interpolationevaluator2dcollection) () = default<br>_Default destructor callable from GPU._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**rank**](#function-rank) () <br>_The number of interpolation dimensions._  |


























## Detailed Description


See the documentation of the primary template.




**Template parameters:**


* `Patches` The patches. 
* `Evaluators` The 2D interpolation evaluators on each patch, in the same order as the patches. 
* `PatchLocator` A operator that finds the patch where a given coordinate is physically located. 
* [**ExtrapolationRule**](namespaceExtrapolationRule.md) The extrapolation rule type for outside of the global domain. 




    
## Public Types Documentation




### typedef CoeffOnPatch 

_Type for MultipatchType: A field of 2D interpolation coefficients for a non-batched interpolation defined on both of the_ [_**Patch**_](structPatch.md) _'s logical dimensions. Needed public for functions on GPU._
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::CoeffOnPatch =  ConstField<data_type, coeff_idx_range_type<Patch>, memory_space>;
```




<hr>



### typedef basis\_type1 

_The discrete dimension representing the basis along first dimension._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::basis_type1 =  ddc::
            type_seq_element_t<0, ddc::to_type_seq_t<typename Traits<Patch>::basis_idx_range_type>>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef basis\_type2 

_The discrete dimension representing the basis along second dimension._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::basis_type2 =  ddc::
            type_seq_element_t<1, ddc::to_type_seq_t<typename Traits<Patch>::basis_idx_range_type>>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef batched\_evaluation\_idx\_range\_type 

_The type of the whole domain representing evaluation points._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::batched_evaluation_idx_range_type = 
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<
                    physical_grids<ValuesField>>;
```





**Template parameters:**


* `ValuesField` A Field or [**DerivField**](classDerivField.md) type storing the evaluated values. 




        

<hr>



### typedef coeff\_idx\_range\_type 

_The type of the 2D index range on which the interpolation coefficients are defined._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::coeff_idx_range_type =  typename Traits<Patch>::coeff_idx_range_type;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef continuous\_dimension\_type1 

_The type of the first evaluation continuous dimension used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::continuous_dimension_type1 = 
            typename evaluation_discrete_dimension_type1<Patch>::continuous_dimension_type;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef continuous\_dimension\_type2 

_The type of the second evaluation continuous dimension used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::continuous_dimension_type2 = 
            typename evaluation_discrete_dimension_type2<Patch>::continuous_dimension_type;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef coord\_type 

_The ND coordinate type corresponding to the evaluation mesh._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::coord_type =  typename InterpolationEvaluatorTraits<FirstEvaluator>::coord_type;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::data_type =  typename InterpolationEvaluatorTraits<FirstEvaluator>::data_type;
```




<hr>



### typedef deriv\_tags 

_TypeSeq of the derivative dimensions (contained in_ [_**DerivField**_](classDerivField.md) _)._
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::deriv_tags =  detail::deriv_sub_set_t<ddc::to_type_seq_t<idx_range_type<ValuesField>>>;
```





**Template parameters:**


* `ValuesField` A Field or [**DerivField**](classDerivField.md) type storing the evaluated values. 




        

<hr>



### typedef evaluation\_discrete\_dimension\_type1 

_The type of the first discrete dimension of interest used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::evaluation_discrete_dimension_type1 =  ddc::type_seq_element_t<0, ddc::to_type_seq_t<evaluation_idx_range_type<Patch>>>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_discrete\_dimension\_type2 

_The type of the second discrete dimension of interest used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::evaluation_discrete_dimension_type2 =  ddc::type_seq_element_t<1, ddc::to_type_seq_t<evaluation_idx_range_type<Patch>>>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_idx\_range\_type 

_The type of the 2D evaluation mesh used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::evaluation_idx_range_type =  typename Traits<Patch>::evaluation_idx_range_type;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_idx\_range\_type1 

_The type of the domain for the 1D evaluation mesh along first dimension used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::evaluation_idx_range_type1 =  IdxRange<evaluation_discrete_dimension_type1<Patch>>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef evaluation\_idx\_range\_type2 

_The type of the domain for the 1D evaluation mesh along second dimension used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::evaluation_idx_range_type2 =  IdxRange<evaluation_discrete_dimension_type2<Patch>>;
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type. 




        

<hr>



### typedef exec\_space 

_The type of the Kokkos execution space used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::exec_space =  typename FirstEvaluator::exec_space;
```




<hr>



### typedef idx\_range\_type 

_The type of the whole discrete domain (dimensions for grids and derivatives)._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::idx_range_type =  typename ValuesField::discrete_domain_type;
```





**Template parameters:**


* `ValuesField` A Field or [**DerivField**](classDerivField.md) type storing the evaluated values. 




        

<hr>



### typedef memory\_space 

_The type of the Kokkos memory space used by this class._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::memory_space =  typename FirstEvaluator::memory_space;
```




<hr>



### typedef physical\_grids 

_TypeSeq of the grid dimensions (without derivative dimensions)._ 
```C++
using InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::physical_grids =  ddc::type_seq_remove_t<
            ddc::to_type_seq_t<idx_range_type<ValuesField>>,
            deriv_tags<ValuesField>>;
```





**Template parameters:**


* `ValuesField` A Field or [**DerivField**](classDerivField.md) type storing the evaluated values. 




        

<hr>
## Public Static Attributes Documentation




### variable n\_patches 

_The number of patches._ 
```C++
std::size_t InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::n_patches;
```




<hr>
## Public Functions Documentation




### function InterpolationEvaluator2DCollection [1/2]

_Instantiate a_ [_**InterpolationEvaluator2DCollection**_](classInterpolationEvaluator2DCollection.md) _._
```C++
inline KOKKOS_FUNCTION InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::InterpolationEvaluator2DCollection (
    Evaluators const &... evaluators,
    PatchLocator const & patch_locator,
    ExtrapolationRule const & extrapolation_rule
) 
```





**Warning:**

The evaluators have to be sorted in the same order as the patches.




**Parameters:**


* `evaluators` The interpolation evaluators for each patch. 
* `patch_locator` An operator to locate a coordinate. The mapping stored in this class has to be invertible. 
* `extrapolation_rule` The extrapolation rules. 




        

<hr>



### function InterpolationEvaluator2DCollection [2/2]

_Default copy operator callable from GPU._ 
```C++
KOKKOS_DEFAULTED_FUNCTION InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::InterpolationEvaluator2DCollection (
    InterpolationEvaluator2DCollection const &
) = default
```




<hr>



### function apply\_evaluator\_deriv\_field 

_Compute the derivatives of a given patch at the coordinates defined on the given patch. Specialisation for patch\_values on_ [_**DerivField**_](classDerivField.md) _. Needed public for functions on GPU._
```C++
template<class StoringPatch, class ValuesField, class... DerivDims>
inline void InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::apply_evaluator_deriv_field (
    Idx< DerivDims... > const & deriv_order,
    ValuesField const & patch_values,
    CoordConstFieldOnPatch< StoringPatch > const & patch_coords,
    MultipatchCoeffs const patches_coeffs
) const
```





**Template parameters:**


* `StoringPatch` [**Patch**](structPatch.md) type where the given coordinates are stored. They are not especially physically located on this patch. 
* `ValuesField` The type of the [**DerivField**](classDerivField.md) of values. 
* `DerivDims` Derivative dimension types, it can be on the first dimension, the second or both. 



**Parameters:**


* `deriv_order` Order (k,l) of the derivatives (e.g. d\_x^k d\_y^l f). 
* `patch_values` [**DerivField**](classDerivField.md) of values and derivatives of the function. 
* `patch_coords` ConstField of coordinates defined on the StoringPatch and where we want to evaluate the function or derivative. 
* `patches_coeffs` MultipatchType of interpolation coefficients on every patches. 




        

<hr>



### function apply\_evaluator\_field 

_Compute the values or the derivatives of a given patch at the coordinates defined on the given patch. Specialisation for patch\_values on DField. Needed public for functions on GPU._ 
```C++
template<class StoringPatch, class ValuesField, class... DerivDims>
inline void InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::apply_evaluator_field (
    Idx< DerivDims... > const & deriv_order,
    ValuesField const & patch_values,
    CoordConstFieldOnPatch< StoringPatch > const & patch_coords,
    MultipatchCoeffs const patches_coeffs
) const
```





**Template parameters:**


* `StoringPatch` [**Patch**](structPatch.md) type where the given coordinates are stored. They are not especially physically located on this patch. 
* `ValuesField` The type of the field of values (DField or [**DerivField**](classDerivField.md)). 
* `DerivDims` Derivative dimension types, it can be on the first dimension, the second or both. 



**Parameters:**


* `deriv_order` Order (k,l) of the derivatives (e.g. d\_x^k d\_y^l f). 
* `patch_values` Field of values of the function or derivative. 
* `patch_coords` ConstField of coordinates defined on the StoringPatch and where we want to evaluate the function or derivative. 
* `patches_coeffs` MultipatchType of interpolation coefficients on every patches. 




        

<hr>



### function apply\_integrate 

_Integrate the spline defined on the given patch._ 
```C++
template<class Patch>
inline void InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::apply_integrate (
    double & integral,
    CoeffOnPatch < Patch > const & coeffs
) const
```





**Template parameters:**


* [**Patch**](structPatch.md) [**Patch**](structPatch.md) type where the integration of the spline is computed. 



**Parameters:**


* `integral` Double, value of the integral of the spline on the given [**Patch**](structPatch.md). 
* `coeffs` ConstField of spline coefficients of the spline on the given [**Patch**](structPatch.md). 




        

<hr>



### function deriv [1/2]

_Differentiate 2D interpolations (described by their interpolation coefficients) at a given coordinate along a specified dimension of interest._ 
```C++
template<class... DerivDims>
inline KOKKOS_FUNCTION data_type InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::deriv (
    Idx< DerivDims... > const & deriv_order,
    Coord< CDim1, CDim2 > const & coord_eval,
    MultipatchCoeffs const & patches_coeffs
) const
```



See InterpolationEvaluator2DCollectionOperator.




**Warning:**

The derivative cannot be computed outside of the domain.




**Template parameters:**


* `DerivDims` Derivative dimension types, it can be on the first dimension, the second or both. 



**Parameters:**


* `deriv_order` Order (k,l) of the derivatives (e.g. d\_x^k d\_y^l f). 
* `coord_eval` The coordinate where the interpolation is differentiated. 
* `patches_coeffs` A MultipatchType of DField storing the 2D interpolation coefficients.



**Returns:**

The derivative of the interpolation at the desired coordinate on the right patch. 





        

<hr>



### function deriv [2/2]

_Differentiate 2D interpolations (described by their interpolation coefficients) on meshes along first dimension of interest._ 
```C++
template<class MultipatchValues, class... DerivDims>
inline void InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::deriv (
    Idx< DerivDims... > const & deriv_order,
    MultipatchValues const & patches_deriv,
    MultipatchCoordField const & patches_coords,
    MultipatchCoeffs const & patches_coeffs
) const
```



See InterpolationEvaluator2DCollectionOperator.




**Warning:**

The derivatives cannot be computed outside of the domain.




**Template parameters:**


* `DerivDims` Derivative dimension types, it can be on the first dimension, the second or both. 



**Parameters:**


* `deriv_order` Order (k,l) of the derivatives (e.g. d\_x^k d\_y^l f(x,y)). 
* `patches_deriv` A MultipatchType of DField or [**DerivField**](classDerivField.md) to store the derivatives of the interpolations at the given coordinates. 
* `patches_coords` A MultipatchType of Field of Coordinate storing the coordinates of the meshes. If patches\_deriv is a [**DerivField**](classDerivField.md), the given coordinates can be the coordinates for all function values. The operator will select the correct coordinates. 
* `patches_coeffs` A MultipatchType of DField storing the 2D interpolation coefficients. 



**Template parameters:**


* `MultipatchValues` A MultipatchField of DField or [**DerivField**](classDerivField.md). 




        

<hr>



### function integrate 

_Integration of splines (described by their spline coefficients)._ 
```C++
inline void InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::integrate (
    DKokkosView_h< n_patches > const & integrals,
    MultipatchCoeffs const & patches_coeffs
) const
```



See InterpolationEvaluator2DCollectionOperator.


This function is only available if all the evaluators are spline evaluators (ddc::SplineEvaluatorND).




**Parameters:**


* `integrals` An Kokkos::View on host containing the integrals of each spline on each patch. 
* `patches_coeffs` A MultipatchType of DField storing the 2D interpolation coefficients. 




        

<hr>



### function operator() 

_Evaluate the 2D interpolations (described by their interpolation coefficients) at a given coordinate._ 
```C++
template<class Coord>
inline KOKKOS_FUNCTION data_type InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::operator() (
    Coord const coord_eval,
    MultipatchCoeffs const & patches_coeffs
) const
```



The interpolation coefficients represent a 2D function defined on the basis of the evaluator of each patch. They can be obtained via various methods, such as using an [**InterpolationBuilder2DCollection**](classInterpolationBuilder2DCollection.md). The coordinate is defined on one patch. But it can be physically located on another patch.


 

**Parameters:**


* `coord_eval` The coordinate where the interpolation is evaluated. 
* `patches_coeffs` A MultipatchType of DField storing the 2D interpolation coefficients. 



**Returns:**

The value of the interpolation at the desired coordinate on the right patch. 





        

<hr>



### function operator() 

_Evaluate 2D interpolations (described by their interpolation coefficients) on meshes._ 
```C++
template<class MultipatchValues>
inline void InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::operator() (
    MultipatchValues const & patches_values,
    MultipatchCoordField const & patches_coords,
    MultipatchCoeffs const & patches_coeffs
) const
```



See InterpolationEvaluator2DCollectionOperator.




**Template parameters:**


* `MultipatchValues` A MultipatchField of DField or [**DerivField**](classDerivField.md). 



**Parameters:**


* `patches_values` A MultipatchType of DField or [**DerivField**](classDerivField.md) to store the values of the interpolations at the given coordinates. 
* `patches_coords` A MultipatchType of Field of Coordinate storing the coordinates of the meshes. 
* `patches_coeffs` A MultipatchType of DField storing the 2D interpolation coefficients. 




        

<hr>



### function ~InterpolationEvaluator2DCollection 

_Default destructor callable from GPU._ 
```C++
KOKKOS_DEFAULTED_FUNCTION InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::~InterpolationEvaluator2DCollection () = default
```




<hr>
## Public Static Functions Documentation




### function rank 

_The number of interpolation dimensions._ 
```C++
static inline std::size_t InterpolationEvaluator2DCollection< ddc::TypeSeq< Patches... >, ddc::TypeSeq< Evaluators... >, PatchLocator, ExtrapolationRule >::rank () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interpolation/interpolation_evaluator_2d_collection.hpp`

