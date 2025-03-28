

# Class DerivField&lt; ElementType, IdxRange&lt; DDims... &gt;, MemorySpace, LayoutStridedPolicy &gt;

**template &lt;class ElementType, class... DDims, class MemorySpace, class LayoutStridedPolicy&gt;**



[**ClassList**](annotated.md) **>** [**DerivField&lt; ElementType, IdxRange&lt; DDims... &gt;, MemorySpace, LayoutStridedPolicy &gt;**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md)



_A class which holds references to chunks of memory describing a field and its derivatives._ [More...](#detailed-description)

* `#include <derivative_field.hpp>`



Inherits the following classes: [DerivFieldCommon](classDerivFieldCommon.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename base\_type::chunk\_type | [**chunk\_type**](#typedef-chunk_type)  <br>_The type of the field stored in the array._  |
| typedef typename base\_type::deriv\_tags | [**deriv\_tags**](#typedef-deriv_tags)  <br>_A type sequence containing all derivatives present in this object._  |
| typedef typename base\_type::discrete\_domain\_type | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The IdxRange on which the chunks in this object are defined._  |
| typedef typename base\_type::discrete\_element\_type | [**discrete\_element\_type**](#typedef-discrete_element_type)  <br>_The Idx which can be used to index this object._  |
| typedef typename base\_type::element\_type | [**element\_type**](#typedef-element_type)  <br>_The type of the elements in the chunks._  |
| typedef typename base\_type::index\_range\_type | [**index\_range\_type**](#typedef-index_range_type)  <br>_The IdxRange on which the fields in this object are defined._  |
| typedef typename detail::strip\_deriv\_t&lt; [**deriv\_tags**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md#typedef-deriv_tags) &gt; | [**physical\_deriv\_grids**](#typedef-physical_deriv_grids)  <br>_A type sequence containing all grid types for which derivatives are present in this object._  |
| typedef typename base\_type::physical\_grids | [**physical\_grids**](#typedef-physical_grids)  <br>_A type sequence containing all the grids on which the fields are defined._  |
| typedef typename base\_type::physical\_idx\_range\_type | [**physical\_idx\_range\_type**](#typedef-physical_idx_range_type)  <br>_The physical index range on which the field is defined._  |
| typedef [**DerivField**](classDerivField.md)&lt; ElementType, IdxRange&lt; DDims... &gt;, MemorySpace, LayoutStridedPolicy &gt; | [**span\_type**](#typedef-span_type)  <br>_The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field._  |
| typedef [**DerivField**](classDerivField.md)&lt; ElementType const, IdxRange&lt; DDims... &gt;, MemorySpace, LayoutStridedPolicy &gt; | [**view\_type**](#typedef-view_type)  <br>_The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field._  |








































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr | [**DerivField**](#function-derivfield-15) ([**DerivField**](classDerivField.md) const & other) = default<br>_Copy-construct a_ [_**DerivField**_](classDerivField.md) _._ |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr | [**DerivField**](#function-derivfield-25) ([**DerivField**](classDerivField.md) && other) = default<br>_Move-construct a_ [_**DerivField**_](classDerivField.md) _._ |
|  constexpr | [**DerivField**](#function-derivfield-35) ([**DerivFieldMem**](classDerivFieldMem.md)&lt; OElementType, [**index\_range\_type**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md#typedef-index_range_type), NDerivs, Allocator &gt; & field) <br>_Constructs a new_ [_**DerivField**_](classDerivField.md) _containing a modifiable view on the data in a_[_**DerivFieldMem**_](classDerivFieldMem.md) _._ |
|  constexpr | [**DerivField**](#function-derivfield-45) ([**DerivFieldMem**](classDerivFieldMem.md)&lt; OElementType, [**index\_range\_type**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md#typedef-index_range_type), NDerivs, Allocator &gt; const & field) <br>_Constructs a new_ [_**DerivField**_](classDerivField.md) _containing a constant view on the data in a_[_**DerivFieldMem**_](classDerivFieldMem.md) _._ |
|  KOKKOS\_FUNCTION constexpr | [**DerivField**](#function-derivfield-55) ([**DerivField**](classDerivField.md)&lt; OElementType, [**index\_range\_type**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md#typedef-index_range_type), MemorySpace, LayoutStridedPolicy &gt; const & field) <br>_Copy construct a_ [_**DerivField**_](classDerivField.md) _. The element type may be changed to a complatible type. (e.g. double -&gt; const double)._ |
|  void | [**deepcopy**](#function-deepcopy-12) ([**DerivField**](classDerivField.md)&lt; OElementType, [**index\_range\_type**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md#typedef-index_range_type), OMemorySpace, OLayoutStridedPolicy &gt; src) <br>_Copy the source_ [_**DerivField**_](classDerivField.md) _into this_[_**DerivField**_](classDerivField.md) _using Kokkos::deep\_copy._ |
|  void | [**deepcopy**](#function-deepcopy-22) (ExecSpace const & execution\_space, [**DerivField**](classDerivField.md)&lt; OElementType, [**index\_range\_type**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md#typedef-index_range_type), OMemorySpace, OLayoutStridedPolicy &gt; src) <br>_Copy the source_ [_**DerivField**_](classDerivField.md) _into this_[_**DerivField**_](classDerivField.md) _using Kokkos::deep\_copy._ |
|  KOKKOS\_FUNCTION constexpr reference | [**operator()**](#function-operator) (DElem... elems) noexcept const<br>_Get an element from a constant field. An Idx describes the element of interest. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested._  |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr [**DerivField**](classDerivField.md) & | [**operator=**](#function-operator_1) ([**DerivField**](classDerivField.md) const & other) = default<br> |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr [**DerivField**](classDerivField.md) & | [**operator=**](#function-operator_2) ([**DerivField**](classDerivField.md) && other) = default<br> |
|  KOKKOS\_FUNCTION constexpr [**view\_type**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md#typedef-view_type) | [**span\_cview**](#function-span_cview) () const<br>_Get a constant_ [_**DerivField**_](classDerivField.md) _of this field._ |
|  KOKKOS\_FUNCTION constexpr [**span\_type**](classDerivField_3_01ElementType_00_01IdxRange_3_01DDims_8_8_8_01_4_00_01MemorySpace_00_01LayoutStridedPolicy_01_4.md#typedef-span_type) | [**span\_view**](#function-span_view) () const<br>_Get a modifiable_ [_**DerivField**_](classDerivField.md) _of this field._ |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~DerivField**](#function-derivfield) () = default<br> |
























































## Detailed Description


The values of the field and the derivatives may be defined on different index ranges, but the underlying mesh must be the same for both.


 

**Template parameters:**


* `ElementType` The type of the elements inside the chunks. 
* `IdxRange<DDims...>` The index range on which the internal fields are defined. This index range is the physical index range on which the values are defined combined with the index range of the derivatives of interest (e.g. IdxRange&lt;Deriv&lt;IDimX&gt;, IDimX, IDimY&gt;). 
* `MemorySpace` The memory space where the data is saved (CPU/GPU). 
* `LayoutStridedPolicy` The way in which the memory is laid out in memory (contiguous in the leading/trailing dimension, strided, etc). 




    
## Public Types Documentation




### typedef chunk\_type 

_The type of the field stored in the array._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::chunk_type =  typename base_type::chunk_type;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC FieldMem types are referred to as Chunk types and Field types are referred to as ChunkSpan/ChunkView. 


        

<hr>



### typedef deriv\_tags 

_A type sequence containing all derivatives present in this object._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::deriv_tags =  typename base_type::deriv_tags;
```




<hr>



### typedef discrete\_domain\_type 

_The IdxRange on which the chunks in this object are defined._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::discrete_domain_type =  typename base_type::discrete_domain_type;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC IdxRange types are referred to as DiscreteDomain types. 


        

<hr>



### typedef discrete\_element\_type 

_The Idx which can be used to index this object._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::discrete_element_type =  typename base_type::discrete_element_type;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC Idx types are referred to as DiscreteElement types. 


        

<hr>



### typedef element\_type 

_The type of the elements in the chunks._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::element_type =  typename base_type::element_type;
```




<hr>



### typedef index\_range\_type 

_The IdxRange on which the fields in this object are defined._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::index_range_type =  typename base_type::index_range_type;
```




<hr>



### typedef physical\_deriv\_grids 

_A type sequence containing all grid types for which derivatives are present in this object._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::physical_deriv_grids =  typename detail::strip_deriv_t<deriv_tags>;
```




<hr>



### typedef physical\_grids 

_A type sequence containing all the grids on which the fields are defined._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::physical_grids =  typename base_type::physical_grids;
```




<hr>



### typedef physical\_idx\_range\_type 

_The physical index range on which the field is defined._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::physical_idx_range_type =  typename base_type::physical_idx_range_type;
```




<hr>



### typedef span\_type 

_The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::span_type =  DerivField<ElementType, IdxRange<DDims...>, MemorySpace, LayoutStridedPolicy>;
```




<hr>



### typedef view\_type 

_The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field._ 
```C++
using DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::view_type =  DerivField<ElementType const, IdxRange<DDims...>, MemorySpace, LayoutStridedPolicy>;
```




<hr>
## Public Functions Documentation




### function DerivField [1/5]

_Copy-construct a_ [_**DerivField**_](classDerivField.md) _._
```C++
KOKKOS_DEFAULTED_FUNCTION constexpr DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::DerivField (
    DerivField const & other
) = default
```





**Parameters:**


* `other` The [**DerivField**](classDerivField.md) being copied. 




        

<hr>



### function DerivField [2/5]

_Move-construct a_ [_**DerivField**_](classDerivField.md) _._
```C++
KOKKOS_DEFAULTED_FUNCTION constexpr DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::DerivField (
    DerivField && other
) = default
```





**Parameters:**


* `other` The [**DerivField**](classDerivField.md) being moved. 




        

<hr>



### function DerivField [3/5]

_Constructs a new_ [_**DerivField**_](classDerivField.md) _containing a modifiable view on the data in a_[_**DerivFieldMem**_](classDerivFieldMem.md) _._
```C++
template<class OElementType, int NDerivs, class Allocator, class>
inline explicit constexpr DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::DerivField (
    DerivFieldMem < OElementType, index_range_type , NDerivs, Allocator > & field
) 
```





**Parameters:**


* `field` The [**DerivFieldMem**](classDerivFieldMem.md) to view. 




        

<hr>



### function DerivField [4/5]

_Constructs a new_ [_**DerivField**_](classDerivField.md) _containing a constant view on the data in a_[_**DerivFieldMem**_](classDerivFieldMem.md) _._
```C++
template<class OElementType, class SFINAEElementType, class, int NDerivs, class Allocator, class>
inline explicit constexpr DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::DerivField (
    DerivFieldMem < OElementType, index_range_type , NDerivs, Allocator > const & field
) 
```





**Parameters:**


* `field` The [**DerivFieldMem**](classDerivFieldMem.md) to view. 




        

<hr>



### function DerivField [5/5]

_Copy construct a_ [_**DerivField**_](classDerivField.md) _. The element type may be changed to a complatible type. (e.g. double -&gt; const double)._
```C++
template<class OElementType>
inline KOKKOS_FUNCTION constexpr DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::DerivField (
    DerivField < OElementType, index_range_type , MemorySpace, LayoutStridedPolicy > const & field
) 
```





**Parameters:**


* `field` The [**DerivField**](classDerivField.md) to be copied. 




        

<hr>



### function deepcopy [1/2]

_Copy the source_ [_**DerivField**_](classDerivField.md) _into this_[_**DerivField**_](classDerivField.md) _using Kokkos::deep\_copy._
```C++
template<class OElementType, class OLayoutStridedPolicy, class OMemorySpace>
inline void DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::deepcopy (
    DerivField < OElementType, index_range_type , OMemorySpace, OLayoutStridedPolicy > src
) 
```





**Parameters:**


* `src` The [**DerivField**](classDerivField.md) containing the data to be copied. 




        

<hr>



### function deepcopy [2/2]

_Copy the source_ [_**DerivField**_](classDerivField.md) _into this_[_**DerivField**_](classDerivField.md) _using Kokkos::deep\_copy._
```C++
template<class ExecSpace, class OElementType, class OMemorySpace, class OLayoutStridedPolicy>
inline void DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::deepcopy (
    ExecSpace const & execution_space,
    DerivField < OElementType, index_range_type , OMemorySpace, OLayoutStridedPolicy > src
) 
```





**Parameters:**


* `execution_space` The execution space on which the copy will be carried out. 
* `src` The [**DerivField**](classDerivField.md) containing the data to be copied. 




        

<hr>



### function operator() 

_Get an element from a constant field. An Idx describes the element of interest. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested._ 
```C++
template<class... DElem>
inline KOKKOS_FUNCTION constexpr reference DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::operator() (
    DElem... elems
) noexcept const
```





**Parameters:**


* `elems` The element of interest.



**Returns:**

The requested element. 





        

<hr>



### function operator= 

```C++
KOKKOS_DEFAULTED_FUNCTION constexpr DerivField & DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::operator= (
    DerivField const & other
) = default
```



Copy-assigns a new value to this [**DerivField**](classDerivField.md), yields a new view to the same data 

**Parameters:**


* `other` the [**DerivField**](classDerivField.md) to copy 



**Returns:**

\*this 





        

<hr>



### function operator= 

```C++
KOKKOS_DEFAULTED_FUNCTION constexpr DerivField & DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::operator= (
    DerivField && other
) = default
```



Move-assigns a new value to this [**DerivField**](classDerivField.md) 

**Parameters:**


* `other` the [**DerivField**](classDerivField.md) to move 



**Returns:**

\*this 





        

<hr>



### function span\_cview 

_Get a constant_ [_**DerivField**_](classDerivField.md) _of this field._
```C++
inline KOKKOS_FUNCTION constexpr view_type DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::span_cview () const
```



This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_const\_field should be used.




**Returns:**

A constant span of this field. 





        

<hr>



### function span\_view 

_Get a modifiable_ [_**DerivField**_](classDerivField.md) _of this field._
```C++
inline KOKKOS_FUNCTION constexpr span_type DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::span_view () const
```



This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_field should be used.




**Returns:**

A span of this field. 





        

<hr>



### function ~DerivField 

```C++
KOKKOS_DEFAULTED_FUNCTION DerivField< ElementType, IdxRange< DDims... >, MemorySpace, LayoutStridedPolicy >::~DerivField () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/derivative_field.hpp`

