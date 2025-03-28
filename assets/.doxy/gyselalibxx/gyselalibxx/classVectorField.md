

# Class VectorField

**template &lt;class ElementType, class IdxRangeType, class VectorIndexSetType, class MemorySpace, class LayoutStridedPolicy&gt;**



[**ClassList**](annotated.md) **>** [**VectorField**](classVectorField.md)



_A class which holds multiple (scalar) fields in order to represent a vector field._ [More...](#detailed-description)

* `#include <vector_field.hpp>`



Inherits the following classes: [VectorFieldCommon](classVectorFieldCommon.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename base\_type::discrete\_domain\_type | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The type of the index range on which the field is defined. This is a DDC keyword used to make this class interchangeable with Field. In DDC IdxRange types are referred to as DiscreteDomain types._  |
| typedef typename base\_type::element\_type | [**element\_type**](#typedef-element_type)  <br>_The type of an element in one of the Fields comprising the_ [_**VectorField**_](classVectorField.md) _._ |
| typedef Field&lt; ElementType, IdxRangeType, MemorySpace, LayoutStridedPolicy &gt; | [**field\_type**](#typedef-field_type)  <br>_Type describing the object which can be extracted from this_ [_**VectorField**_](classVectorField.md) _using the get&lt;&gt; function._ |
| typedef [**discrete\_domain\_type**](classVectorField.md#typedef-discrete_domain_type) | [**index\_range\_type**](#typedef-index_range_type)  <br>_The IdxRange on which the fields in this object are defined._  |
| typedef LayoutStridedPolicy | [**layout\_type**](#typedef-layout_type)  <br>_Type describing the way in which the data is laid out in the Field memory._  |
| typedef typename field\_type::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |
| typedef [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; | [**span\_type**](#typedef-span_type)  <br>_A type which can hold a modifiable reference to a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._ |
| typedef [**VectorField**](classVectorField.md)&lt; const ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; | [**view\_type**](#typedef-view_type)  <br>_A type which can hold a constant reference to a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _. This is a DDC keyword used to make this class interchangeable with Field._ |








































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr | [**VectorField**](#function-vectorfield-513) () = default<br>_Empty_ [_**VectorField**_](classVectorField.md) _._ |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr | [**VectorField**](#function-vectorfield-613) ([**VectorField**](classVectorField.md) const & other) = default<br> |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr | [**VectorField**](#function-vectorfield-713) ([**VectorField**](classVectorField.md) && other) = default<br> |
|  KOKKOS\_FUNCTION constexpr | [**VectorField**](#function-vectorfield-813) ([**VectorFieldMem**](classVectorFieldMem.md)&lt; OElementType, IdxRangeType, VectorIndexSetType, Allocator &gt; & other) noexcept<br> |
|  KOKKOS\_FUNCTION constexpr | [**VectorField**](#function-vectorfield-913) ([**VectorFieldMem**](classVectorFieldMem.md)&lt; OElementType, IdxRangeType, VectorIndexSetType, Allocator &gt; const & other) noexcept<br> |
|   | [**VectorField**](#function-vectorfield-1013) ([**VectorFieldMem**](classVectorFieldMem.md)&lt; OElementType, IdxRangeType, VectorIndexSetType, Allocator &gt; && other) = delete<br> |
|  KOKKOS\_FUNCTION constexpr | [**VectorField**](#function-vectorfield-1113) ([**VectorField**](classVectorField.md)&lt; OElementType, [**index\_range\_type**](classVectorField.md#typedef-index_range_type), VectorIndexSetType, MemorySpace, LayoutStridedPolicy &gt; const & other) noexcept<br> |
|  KOKKOS\_FUNCTION | [**VectorField**](#function-vectorfield-1213) ([**index\_range\_type**](classVectorField.md#typedef-index_range_type) const & idx\_range, OElementType \*... ptr) <br> |
|  KOKKOS\_FUNCTION constexpr | [**VectorField**](#function-vectorfield-1313) (FieldType... fields) <br> |
|  constexpr chunk\_span\_type | [**get**](#function-get) () noexcept const<br>_Get the Field describing the component in the QueryTag direction._  |
|  KOKKOS\_FUNCTION [**element\_type**](classVectorField.md#typedef-element_type) const | [**operator()**](#function-operator) (ddc::DiscreteElement&lt; ODDims &gt; const &... delems) noexcept const<br> |
|  KOKKOS\_FUNCTION [**element\_type**](classVectorField.md#typedef-element_type) const | [**operator()**](#function-operator_1) (Idx&lt; ODDims... &gt; const & delems) noexcept const<br> |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr [**VectorField**](classVectorField.md) & | [**operator=**](#function-operator_2) ([**VectorField**](classVectorField.md) const & other) = default<br> |
|  KOKKOS\_DEFAULTED\_FUNCTION constexpr [**VectorField**](classVectorField.md) & | [**operator=**](#function-operator_3) ([**VectorField**](classVectorField.md) && other) = default<br> |
|  constexpr auto | [**operator[]**](#function-operator_4) (Idx&lt; QueryDDims... &gt; const & slice\_spec) <br>_Slice out some dimensions._  |
|  constexpr auto | [**operator[]**](#function-operator_5) (IdxRange&lt; QueryDDims... &gt; const & oidx\_range) <br>_Slice out some dimensions._  |
|  constexpr [**view\_type**](classVectorField.md#typedef-view_type) | [**span\_cview**](#function-span_cview) () const<br> |
|  constexpr [**span\_type**](classVectorField.md#typedef-span_type) | [**span\_view**](#function-span_view) () const<br> |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~VectorField**](#function-vectorfield) () = default<br>[_**VectorField**_](classVectorField.md) _destructor._ |
























































## Detailed Description


Pre-declaration of [**VectorField**](classVectorField.md).




**Template parameters:**


* `ElementType` The data type of a scalar element of the vector field. 
* `IdxRangeType` 
* `VectorIndexSetType` A VectorIndexSet describing the dimensions described by the scalar elements of a vector field element. 
* `MemorySpace` The memory space (CPU/GPU). 
* `LayoutStridedPolicy` The memory layout. See DDC. 




    
## Public Types Documentation




### typedef discrete\_domain\_type 

_The type of the index range on which the field is defined. This is a DDC keyword used to make this class interchangeable with Field. In DDC IdxRange types are referred to as DiscreteDomain types._ 
```C++
using VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy >::discrete_domain_type =  typename base_type::discrete_domain_type;
```




<hr>



### typedef element\_type 

_The type of an element in one of the Fields comprising the_ [_**VectorField**_](classVectorField.md) _._
```C++
using VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy >::element_type =  typename base_type::element_type;
```




<hr>



### typedef field\_type 

_Type describing the object which can be extracted from this_ [_**VectorField**_](classVectorField.md) _using the get&lt;&gt; function._
```C++
using VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy >::field_type =  Field<ElementType, IdxRangeType, MemorySpace, LayoutStridedPolicy>;
```




<hr>



### typedef index\_range\_type 

_The IdxRange on which the fields in this object are defined._ 
```C++
using VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy >::index_range_type =  discrete_domain_type;
```




<hr>



### typedef layout\_type 

_Type describing the way in which the data is laid out in the Field memory._ 
```C++
using VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy >::layout_type =  LayoutStridedPolicy;
```



Type describing the way in which the data is laid out in the Field memory. I.e. it describes whether it is contiguous or not. 


        

<hr>



### typedef memory\_space 

_The type of the memory space where the field is saved (CPU vs GPU)._ 
```C++
using VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy >::memory_space =  typename field_type::memory_space;
```




<hr>



### typedef span\_type 

_A type which can hold a modifiable reference to a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._
```C++
using VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy >::span_type =  VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy>;
```



A type which can hold a reference to a [**VectorFieldMem**](classVectorFieldMem.md). If this object is modifiable then so is the span type. This is a DDC keyword used to make this class interchangeable with Field. 


        

<hr>



### typedef view\_type 

_A type which can hold a constant reference to a_ [_**VectorFieldMem**_](classVectorFieldMem.md) _. This is a DDC keyword used to make this class interchangeable with Field._
```C++
using VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy >::view_type =  VectorField< const ElementType, IdxRangeType, VectorIndexSetType, MemorySpace, LayoutStridedPolicy>;
```




<hr>
## Public Functions Documentation




### function VectorField [5/13]

_Empty_ [_**VectorField**_](classVectorField.md) _._
```C++
KOKKOS_DEFAULTED_FUNCTION constexpr VectorField::VectorField () = default
```




<hr>



### function VectorField [6/13]

```C++
KOKKOS_DEFAULTED_FUNCTION constexpr VectorField::VectorField (
    VectorField const & other
) = default
```



Constructs a new [**VectorField**](classVectorField.md) by copy, yields a new view to the same data 

**Parameters:**


* `other` the [**VectorField**](classVectorField.md) to copy 




        

<hr>



### function VectorField [7/13]

```C++
KOKKOS_DEFAULTED_FUNCTION constexpr VectorField::VectorField (
    VectorField && other
) = default
```



Constructs a new [**VectorField**](classVectorField.md) by move 

**Parameters:**


* `other` the [**VectorField**](classVectorField.md) to move 




        

<hr>



### function VectorField [8/13]

```C++
template<class OElementType, class Allocator>
inline explicit KOKKOS_FUNCTION constexpr VectorField::VectorField (
    VectorFieldMem < OElementType, IdxRangeType, VectorIndexSetType, Allocator > & other
) noexcept
```



Constructs a new [**VectorField**](classVectorField.md) from a [**VectorFieldMem**](classVectorFieldMem.md), yields a new view to the same data 

**Parameters:**


* `other` the [**VectorFieldMem**](classVectorFieldMem.md) to view 




        

<hr>



### function VectorField [9/13]

```C++
template<class OElementType, class SFINAEElementType, class, class Allocator>
inline explicit KOKKOS_FUNCTION constexpr VectorField::VectorField (
    VectorFieldMem < OElementType, IdxRangeType, VectorIndexSetType, Allocator > const & other
) noexcept
```



Constructs a new [**VectorField**](classVectorField.md) from a [**VectorFieldMem**](classVectorFieldMem.md), yields a new view to the same data 

**Parameters:**


* `other` the [**VectorFieldMem**](classVectorFieldMem.md) to view 




        

<hr>



### function VectorField [10/13]

```C++
template<class OElementType, class Allocator>
VectorField::VectorField (
    VectorFieldMem < OElementType, IdxRangeType, VectorIndexSetType, Allocator > && other
) = delete
```




<hr>



### function VectorField [11/13]

```C++
template<class OElementType>
inline KOKKOS_FUNCTION constexpr VectorField::VectorField (
    VectorField < OElementType, index_range_type , VectorIndexSetType, MemorySpace, LayoutStridedPolicy > const & other
) noexcept
```



Constructs a new [**VectorField**](classVectorField.md) by copy of a chunk, yields a new view to the same data 

**Parameters:**


* `other` the [**VectorField**](classVectorField.md) to move 




        

<hr>



### function VectorField [12/13]

```C++
template<class... OElementType, class, class>
inline KOKKOS_FUNCTION VectorField::VectorField (
    index_range_type const & idx_range,
    OElementType *... ptr
) 
```



Constructs a new [**VectorField**](classVectorField.md) from scratch 

**Parameters:**


* `ptr` the allocation pointer to the data 
* `idx_range` the index range that sustains the view 




        

<hr>



### function VectorField [13/13]

```C++
template<class... FieldType, class>
inline KOKKOS_FUNCTION constexpr VectorField::VectorField (
    FieldType... fields
) 
```



Constructs a new [**VectorField**](classVectorField.md) containing references to Field. 

**Parameters:**


* `fields` The Fields. 




        

<hr>



### function get 

_Get the Field describing the component in the QueryTag direction._ 
```C++
template<class QueryTag>
inline constexpr chunk_span_type VectorField::get () noexcept const
```





**Returns:**

The field in the specified direction. 





        

<hr>



### function operator() 

```C++
template<class... ODDims>
inline KOKKOS_FUNCTION element_type const VectorField::operator() (
    ddc::DiscreteElement< ODDims > const &... delems
) noexcept const
```



Element access using a list of Idxs 

**Parameters:**


* `delems` 1D discrete coordinates 



**Returns:**

copy of this element 





        

<hr>



### function operator() 

```C++
template<class... ODDims, class>
inline KOKKOS_FUNCTION element_type const VectorField::operator() (
    Idx< ODDims... > const & delems
) noexcept const
```



Element access using a multi-dimensional Idx 

**Parameters:**


* `delems` discrete coordinates 



**Returns:**

copy of this element 





        

<hr>



### function operator= 

```C++
KOKKOS_DEFAULTED_FUNCTION constexpr VectorField & VectorField::operator= (
    VectorField const & other
) = default
```



Copy-assigns a new value to this [**VectorField**](classVectorField.md), yields a new view to the same data 

**Parameters:**


* `other` the [**VectorField**](classVectorField.md) to copy 



**Returns:**

\*this 





        

<hr>



### function operator= 

```C++
KOKKOS_DEFAULTED_FUNCTION constexpr VectorField & VectorField::operator= (
    VectorField && other
) = default
```



Move-assigns a new value to this [**VectorField**](classVectorField.md) 

**Parameters:**


* `other` the [**VectorField**](classVectorField.md) to move 



**Returns:**

\*this 





        

<hr>



### function operator[] 

_Slice out some dimensions._ 
```C++
template<class... QueryDDims>
inline constexpr auto VectorField::operator[] (
    Idx< QueryDDims... > const & slice_spec
) 
```



Get the [**VectorFieldMem**](classVectorFieldMem.md) on the reduced index range which is obtained by indexing the dimensions QueryDDims at the position slice\_spec.




**Parameters:**


* `slice_spec` The slice describing the index range of interest.



**Returns:**

A reference to the vector field on the sliced index range. 





        

<hr>



### function operator[] 

_Slice out some dimensions._ 
```C++
template<class... QueryDDims>
inline constexpr auto VectorField::operator[] (
    IdxRange< QueryDDims... > const & oidx_range
) 
```



Get the [**VectorFieldMem**](classVectorFieldMem.md) on the reduced index range passed as an argument.




**Parameters:**


* `oidx_range` The index range of interest.



**Returns:**

A reference to the vector field on the sliced index range. 





        

<hr>



### function span\_cview 

```C++
inline constexpr view_type VectorField::span_cview () const
```



Get a constant reference to the vector field referred to by this vector field span.


This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_const\_field should be used.




**Returns:**

A constant reference to the vector field. 





        

<hr>



### function span\_view 

```C++
inline constexpr span_type VectorField::span_view () const
```



Get a modifiable reference to the vector field referred to by this vector field span.


This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_field should be used.




**Returns:**

A constant reference to the vector field. 





        

<hr>



### function ~VectorField 

[_**VectorField**_](classVectorField.md) _destructor._
```C++
KOKKOS_DEFAULTED_FUNCTION VectorField::~VectorField () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/vector_field.hpp`

