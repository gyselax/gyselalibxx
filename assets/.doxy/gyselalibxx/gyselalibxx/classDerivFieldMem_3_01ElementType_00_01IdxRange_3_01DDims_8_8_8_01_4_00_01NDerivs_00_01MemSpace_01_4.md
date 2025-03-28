

# Class DerivFieldMem&lt; ElementType, IdxRange&lt; DDims... &gt;, NDerivs, MemSpace &gt;

**template &lt;class ElementType, class... DDims, int NDerivs, class MemSpace&gt;**



[**ClassList**](annotated.md) **>** [**DerivFieldMem&lt; ElementType, IdxRange&lt; DDims... &gt;, NDerivs, MemSpace &gt;**](classDerivFieldMem_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01NDerivs_00_01MemSpace_01_4.md)



_A class which holds a chunk of memory describing a field and its derivatives._ [More...](#detailed-description)

* `#include <derivative_field_mem.hpp>`



Inherits the following classes: [DerivFieldCommon](classDerivFieldCommon.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename base\_type::chunk\_type | [**chunk\_type**](#typedef-chunk_type)  <br>_The type of the field stored in the array._  |
| typedef typename base\_type::deriv\_tags | [**deriv\_tags**](#typedef-deriv_tags)  <br>_A type sequence containing all derivatives present in this object._  |
| typedef typename base\_type::discrete\_domain\_type | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The IdxRange on which the chunks in this object are defined._  |
| typedef typename base\_type::discrete\_element\_type | [**discrete\_element\_type**](#typedef-discrete_element_type)  <br>_The Idx which can be used to index this object._  |
| typedef typename base\_type::element\_type | [**element\_type**](#typedef-element_type)  <br>_The type of the elements in the chunks._  |
| typedef typename base\_type::physical\_deriv\_grids | [**physical\_deriv\_grids**](#typedef-physical_deriv_grids)  <br>_A type sequence containing all dimensions for which derivatives are present in this object._  |
| typedef typename base\_type::physical\_grids | [**physical\_grids**](#typedef-physical_grids)  <br>_A type sequence containing all the physical dimensions on which the chunks are defined._  |
| typedef typename base\_type::physical\_idx\_range\_type | [**physical\_idx\_range\_type**](#typedef-physical_idx_range_type)  <br>_The physical index range on which the field is defined._  |
| typedef [**DerivField**](classDerivField.md)&lt; ElementType, IdxRange&lt; DDims... &gt;, typename chunk\_type::memory\_space, typename chunk\_type::layout\_type &gt; | [**span\_type**](#typedef-span_type)  <br>_The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field._  |
| typedef [**DerivField**](classDerivField.md)&lt; ElementType const, IdxRange&lt; DDims... &gt;, typename chunk\_type::memory\_space, typename chunk\_type::layout\_type &gt; | [**view\_type**](#typedef-view_type)  <br>_The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field._  |








































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**DerivFieldMem**](#function-derivfieldmem-12) ([**physical\_idx\_range\_type**](classDerivFieldMem_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01NDerivs_00_01MemSpace_01_4.md#typedef-physical_idx_range_type) val\_idx\_range, [**IdxRangeSlice**](classIdxRangeSlice.md)&lt; DerivDoms &gt;... m\_deriv\_idx\_range) <br>_The constructor for_ [_**DerivFieldMem**_](classDerivFieldMem.md) _. The constructor initialises the chunks using the provided index ranges._ |
|   | [**DerivFieldMem**](#function-derivfieldmem-22) (allocator\_type allocator, [**physical\_idx\_range\_type**](classDerivFieldMem_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01NDerivs_00_01MemSpace_01_4.md#typedef-physical_idx_range_type) val\_idx\_range, [**IdxRangeSlice**](classIdxRangeSlice.md)&lt; DerivDoms &gt;... m\_deriv\_idx\_range) <br>_The constructor for_ [_**DerivFieldMem**_](classDerivFieldMem.md) _. The constructor initialises the chunks using the provided index ranges._ |
|  [**element\_type**](classDerivFieldMem_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01NDerivs_00_01MemSpace_01_4.md#typedef-element_type) & | [**operator()**](#function-operator) (DElem... elems) noexcept<br>_Get a modifiable reference to an element from a constant field. A Idx describes the element of interest. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested._  |
|  [**element\_type**](classDerivFieldMem_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01NDerivs_00_01MemSpace_01_4.md#typedef-element_type) const & | [**operator()**](#function-operator_1) (DElem... elems) noexcept const<br>_Get an element from a constant field. A Idx describes the element of interest. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested._  |
|  [**DerivFieldMem**](classDerivFieldMem.md) & | [**operator=**](#function-operator_2) ([**DerivFieldMem**](classDerivFieldMem.md) const & other) = delete<br>_Deleted copy operator._  |
|  [**DerivFieldMem**](classDerivFieldMem.md) & | [**operator=**](#function-operator_3) ([**DerivFieldMem**](classDerivFieldMem.md) && other) = default<br> |
|  [**view\_type**](classDerivFieldMem_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01NDerivs_00_01MemSpace_01_4.md#typedef-view_type) | [**span\_cview**](#function-span_cview) () const<br>_Get a constant_ [_**DerivField**_](classDerivField.md) _of this field._ |
|  [**view\_type**](classDerivFieldMem_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01NDerivs_00_01MemSpace_01_4.md#typedef-view_type) | [**span\_view**](#function-span_view-12) () const<br>_Get a constant_ [_**DerivField**_](classDerivField.md) _of this field._ |
|  [**span\_type**](classDerivFieldMem_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01NDerivs_00_01MemSpace_01_4.md#typedef-span_type) | [**span\_view**](#function-span_view-22) () <br>_Get a modifiable_ [_**DerivField**_](classDerivField.md) _of this field._ |
|   | [**~DerivFieldMem**](#function-derivfieldmem) () = default<br>_Defaulted destructor._  |
























































## Detailed Description


The values of the field and the derivatives may be defined on different index ranges, but the underlying mesh must be the same for both.


 

**Template parameters:**


* `ElementType` The type of the elements inside the chunks. 
* `IdxRange<DDims...>` The index range on which the internal fields are defined. This index range is the physical index range on which the values are defined combined with the index range of the derivatives of interest (e.g. IdxRange&lt;Deriv&lt;IDimX&gt;, IDimX, IDimY&gt;). 
* `NDerivs` The number of derivatives which are defined in the dimensions where derivatives appear. 
* `MemSpace` The memory space where the data will be saved. 




    
## Public Types Documentation




### typedef chunk\_type 

_The type of the field stored in the array._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::chunk_type =  typename base_type::chunk_type;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC FieldMem types are referred to as Chunk types and Field types are referred to as ChunkSpan/ChunkView. 


        

<hr>



### typedef deriv\_tags 

_A type sequence containing all derivatives present in this object._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::deriv_tags =  typename base_type::deriv_tags;
```




<hr>



### typedef discrete\_domain\_type 

_The IdxRange on which the chunks in this object are defined._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::discrete_domain_type =  typename base_type::discrete_domain_type;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC IdxRange types are referred to as DiscreteDomain types. 


        

<hr>



### typedef discrete\_element\_type 

_The Idx which can be used to index this object._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::discrete_element_type =  typename base_type::discrete_element_type;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC Idx types are referred to as DiscreteElement types. 


        

<hr>



### typedef element\_type 

_The type of the elements in the chunks._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::element_type =  typename base_type::element_type;
```




<hr>



### typedef physical\_deriv\_grids 

_A type sequence containing all dimensions for which derivatives are present in this object._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::physical_deriv_grids =  typename base_type::physical_deriv_grids;
```




<hr>



### typedef physical\_grids 

_A type sequence containing all the physical dimensions on which the chunks are defined._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::physical_grids =  typename base_type::physical_grids;
```




<hr>



### typedef physical\_idx\_range\_type 

_The physical index range on which the field is defined._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::physical_idx_range_type =  typename base_type::physical_idx_range_type;
```




<hr>



### typedef span\_type 

_The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::span_type =  DerivField< ElementType, IdxRange<DDims...>, typename chunk_type::memory_space, typename chunk_type::layout_type>;
```




<hr>



### typedef view\_type 

_The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field._ 
```C++
using DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::view_type =  DerivField< ElementType const, IdxRange<DDims...>, typename chunk_type::memory_space, typename chunk_type::layout_type>;
```




<hr>
## Public Functions Documentation




### function DerivFieldMem [1/2]

_The constructor for_ [_**DerivFieldMem**_](classDerivFieldMem.md) _. The constructor initialises the chunks using the provided index ranges._
```C++
template<class... DerivDoms>
inline DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::DerivFieldMem (
    physical_idx_range_type val_idx_range,
    IdxRangeSlice < DerivDoms >... m_deriv_idx_range
) 
```





**Parameters:**


* `val_idx_range` The index range on which the values of the field are defined. 
* `m_deriv_idx_range` The 1D sub-index ranges on which the derivatives of the field are defined. 




        

<hr>



### function DerivFieldMem [2/2]

_The constructor for_ [_**DerivFieldMem**_](classDerivFieldMem.md) _. The constructor initialises the chunks using the provided index ranges._
```C++
template<class... DerivDoms>
inline DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::DerivFieldMem (
    allocator_type allocator,
    physical_idx_range_type val_idx_range,
    IdxRangeSlice < DerivDoms >... m_deriv_idx_range
) 
```





**Parameters:**


* `allocator` The object which allocates the memory on the CPU or GPU. 
* `val_idx_range` The index range on which the values of the field are defined. 
* `m_deriv_idx_range` The 1D sub-index ranges on which the derivatives of the field are defined. 




        

<hr>



### function operator() 

_Get a modifiable reference to an element from a constant field. A Idx describes the element of interest. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested._ 
```C++
template<class... DElem>
inline element_type & DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::operator() (
    DElem... elems
) noexcept
```





**Parameters:**


* `elems` The element of interest.



**Returns:**

The requested element. 





        

<hr>



### function operator() 

_Get an element from a constant field. A Idx describes the element of interest. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested._ 
```C++
template<class... DElem>
inline element_type const & DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::operator() (
    DElem... elems
) noexcept const
```





**Parameters:**


* `elems` The element of interest.



**Returns:**

The requested element. 





        

<hr>



### function operator= 

_Deleted copy operator._ 
```C++
DerivFieldMem & DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::operator= (
    DerivFieldMem const & other
) = delete
```




<hr>



### function operator= 

```C++
DerivFieldMem & DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::operator= (
    DerivFieldMem && other
) = default
```



Move-assigns a new value to this field 

**Parameters:**


* `other` the FieldMem to move 



**Returns:**

\*this 





        

<hr>



### function span\_cview 

_Get a constant_ [_**DerivField**_](classDerivField.md) _of this field._
```C++
inline view_type DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::span_cview () const
```



This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_const\_field should be used.




**Returns:**

A constant span of this field. 





        

<hr>



### function span\_view [1/2]

_Get a constant_ [_**DerivField**_](classDerivField.md) _of this field._
```C++
inline view_type DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::span_view () const
```



This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_field should be used.




**Returns:**

A constant span of this field. 





        

<hr>



### function span\_view [2/2]

_Get a modifiable_ [_**DerivField**_](classDerivField.md) _of this field._
```C++
inline span_type DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::span_view () 
```



This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_field should be used.




**Returns:**

A span of this field. 





        

<hr>



### function ~DerivFieldMem 

_Defaulted destructor._ 
```C++
DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::~DerivFieldMem () = default
```




<hr>## Friends Documentation





### friend DerivField 

```C++
template<class, class, class, class>
class DerivFieldMem< ElementType, IdxRange< DDims... >, NDerivs, MemSpace >::DerivField (
    DerivField
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/derivative_field_mem.hpp`

