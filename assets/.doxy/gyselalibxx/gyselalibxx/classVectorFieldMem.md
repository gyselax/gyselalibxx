

# Class VectorFieldMem

**template &lt;class ElementType, class IdxRangeType, class VectorIndexSetType, class MemSpace&gt;**



[**ClassList**](annotated.md) **>** [**VectorFieldMem**](classVectorFieldMem.md)



_Pre-declaration of_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._[More...](#detailed-description)

* `#include <vector_field_mem.hpp>`



Inherits the following classes: [VectorFieldCommon](classVectorFieldCommon.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::KokkosAllocator&lt; ElementType, MemSpace &gt; | [**Allocator**](#typedef-allocator)  <br>_The type of allocator that will be used to allocate the data._  |
| typedef FieldMem&lt; ElementType, IdxRangeType, MemSpace &gt; | [**chunk\_type**](#typedef-chunk_type)  <br>_Type describing the object which can be extracted from this_ [_**VectorFieldMem**_](classVectorFieldMem.md) _using the get&lt;&gt; function._ |
| typedef typename base\_type::discrete\_domain\_type | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The type of the index range on which the field is defined. This is a DDC keyword used to make this class interchangeable with Field. In DDC IdxRange types are referred to as DiscreteDomain types._  |
| typedef [**discrete\_domain\_type**](classVectorFieldMem.md#typedef-discrete_domain_type) | [**index\_range\_type**](#typedef-index_range_type)  <br>_The IdxRange on which the fields in this object are defined._  |
| typedef typename chunk\_type::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |
| typedef [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, VectorIndexSetType, MemSpace, Kokkos::layout\_right &gt; | [**span\_type**](#typedef-span_type)  <br>_A type which can hold a reference to this_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._ |
| typedef [**VectorField**](classVectorField.md)&lt; const ElementType, IdxRangeType, VectorIndexSetType, MemSpace, Kokkos::layout\_right &gt; | [**view\_type**](#typedef-view_type)  <br>_A type which can hold a constant reference to this_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._ |








































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**VectorFieldMem**](#function-vectorfieldmem-25) () = default<br>_Empty_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._ |
|   | [**VectorFieldMem**](#function-vectorfieldmem-35) ([**index\_range\_type**](classVectorFieldMem.md#typedef-index_range_type) const & idx\_range, [**Allocator**](classVectorFieldMem.md#typedef-allocator) allocator=[**Allocator**](classVectorFieldMem.md#typedef-allocator)()) <br> |
|   | [**VectorFieldMem**](#function-vectorfieldmem-45) ([**VectorFieldMem**](classVectorFieldMem.md) const & other) = delete<br>_Deleted: use deepcopy instead._  |
|   | [**VectorFieldMem**](#function-vectorfieldmem-55) ([**VectorFieldMem**](classVectorFieldMem.md) && other) = default<br> |
|  constexpr chunk\_span\_type | [**get**](#function-get-12) () noexcept<br>_Get the Field describing the component in the QueryTag direction._  |
|  constexpr chunk\_view\_type | [**get**](#function-get-22) () noexcept const<br>_Get the Field describing the component in the QueryTag direction._  |
|  element\_type const | [**operator()**](#function-operator) (ddc::DiscreteElement&lt; ODDims &gt; const &... delems) noexcept const<br> |
|  element\_type const | [**operator()**](#function-operator_1) (Idx&lt; ODDims... &gt; const & delems) noexcept const<br> |
|  [**VectorFieldMem**](classVectorFieldMem.md) & | [**operator=**](#function-operator_2) ([**VectorFieldMem**](classVectorFieldMem.md) const & other) = delete<br>_Deleted: use deepcopy instead._  |
|  [**VectorFieldMem**](classVectorFieldMem.md) & | [**operator=**](#function-operator_3) ([**VectorFieldMem**](classVectorFieldMem.md) && other) = default<br> |
|  auto | [**operator[]**](#function-operator_4) (Idx&lt; QueryDDims... &gt; const & slice\_spec) const<br>_Slice out some dimensions._  |
|  auto | [**operator[]**](#function-operator_5) (Idx&lt; QueryDDims... &gt; const & slice\_spec) <br>_Slice out some dimensions._  |
|  auto | [**operator[]**](#function-operator_6) (IdxRange&lt; QueryDDims... &gt; const & oidx\_range) const<br>_Slice out some dimensions._  |
|  auto | [**operator[]**](#function-operator_7) (IdxRange&lt; QueryDDims... &gt; const & oidx\_range) <br>_Slice out some dimensions._  |
|  [**view\_type**](classVectorFieldMem.md#typedef-view_type) | [**span\_cview**](#function-span_cview) () const<br> |
|  [**view\_type**](classVectorFieldMem.md#typedef-view_type) | [**span\_view**](#function-span_view-12) () const<br> |
|  [**span\_type**](classVectorFieldMem.md#typedef-span_type) | [**span\_view**](#function-span_view-22) () <br> |
|   | [**~VectorFieldMem**](#function-vectorfieldmem) () = default<br> |
























































## Detailed Description


A class which describes the storage for a vector field.


A class which describes the storage for a vector field. In other words a class which maps a position on an index range to a vector (x,y,z,...). This is done by storing the values at the positions in individual FieldMems.




**Template parameters:**


* `ElementType` The data type of a scalar element of the vector field. 
* `IdxRangeType` 
* `VectorIndexSetType` A VectorIndexSet describing the dimensions described by the scalar elements of a vector field element. 
* `MemSpace` The type describing where the memory is allocated. See DDC. 




    
## Public Types Documentation




### typedef Allocator 

_The type of allocator that will be used to allocate the data._ 
```C++
using VectorFieldMem< ElementType, IdxRangeType, VectorIndexSetType, MemSpace >::Allocator =  ddc::KokkosAllocator<ElementType, MemSpace>;
```




<hr>



### typedef chunk\_type 

_Type describing the object which can be extracted from this_ [_**VectorFieldMem**_](classVectorFieldMem.md) _using the get&lt;&gt; function._
```C++
using VectorFieldMem< ElementType, IdxRangeType, VectorIndexSetType, MemSpace >::chunk_type =  FieldMem<ElementType, IdxRangeType, MemSpace>;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC FieldMem types are referred to as Chunk types and Field types are referred to as ChunkSpan/ChunkView. 


        

<hr>



### typedef discrete\_domain\_type 

_The type of the index range on which the field is defined. This is a DDC keyword used to make this class interchangeable with Field. In DDC IdxRange types are referred to as DiscreteDomain types._ 
```C++
using VectorFieldMem< ElementType, IdxRangeType, VectorIndexSetType, MemSpace >::discrete_domain_type =  typename base_type::discrete_domain_type;
```




<hr>



### typedef index\_range\_type 

_The IdxRange on which the fields in this object are defined._ 
```C++
using VectorFieldMem< ElementType, IdxRangeType, VectorIndexSetType, MemSpace >::index_range_type =  discrete_domain_type;
```




<hr>



### typedef memory\_space 

_The type of the memory space where the field is saved (CPU vs GPU)._ 
```C++
using VectorFieldMem< ElementType, IdxRangeType, VectorIndexSetType, MemSpace >::memory_space =  typename chunk_type::memory_space;
```




<hr>



### typedef span\_type 

_A type which can hold a reference to this_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._
```C++
using VectorFieldMem< ElementType, IdxRangeType, VectorIndexSetType, MemSpace >::span_type =  VectorField< ElementType, IdxRangeType, VectorIndexSetType, MemSpace, Kokkos::layout_right>;
```



This is a DDC keyword used to make this class interchangeable with Field. 


        

<hr>



### typedef view\_type 

_A type which can hold a constant reference to this_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._
```C++
using VectorFieldMem< ElementType, IdxRangeType, VectorIndexSetType, MemSpace >::view_type =  VectorField< const ElementType, IdxRangeType, VectorIndexSetType, MemSpace, Kokkos::layout_right>;
```



This is a DDC keyword used to make this class interchangeable with Field. 


        

<hr>
## Public Functions Documentation




### function VectorFieldMem [2/5]

_Empty_ [_**VectorFieldMem**_](classVectorFieldMem.md) _._
```C++
VectorFieldMem::VectorFieldMem () = default
```




<hr>



### function VectorFieldMem [3/5]

```C++
inline explicit VectorFieldMem::VectorFieldMem (
    index_range_type const & idx_range,
    Allocator allocator=Allocator ()
) 
```



Construct a [**VectorFieldMem**](classVectorFieldMem.md) on an index range with uninitialized values




**Parameters:**


* `idx_range` The index range on which the chunk will be defined. 
* `allocator` An optional allocator used to create the chunks. 




        

<hr>



### function VectorFieldMem [4/5]

_Deleted: use deepcopy instead._ 
```C++
VectorFieldMem::VectorFieldMem (
    VectorFieldMem const & other
) = delete
```




<hr>



### function VectorFieldMem [5/5]

```C++
VectorFieldMem::VectorFieldMem (
    VectorFieldMem && other
) = default
```



Constructs a new [**VectorFieldMem**](classVectorFieldMem.md) by move 

**Parameters:**


* `other` the [**VectorFieldMem**](classVectorFieldMem.md) to move 




        

<hr>



### function get [1/2]

_Get the Field describing the component in the QueryTag direction._ 
```C++
template<class QueryTag>
inline constexpr chunk_span_type VectorFieldMem::get () noexcept
```





**Returns:**

The field in the specified direction. 





        

<hr>



### function get [2/2]

_Get the Field describing the component in the QueryTag direction._ 
```C++
template<class QueryTag>
inline constexpr chunk_view_type VectorFieldMem::get () noexcept const
```





**Returns:**

The constant field in the specified direction. 





        

<hr>



### function operator() 

```C++
template<class... ODDims>
inline element_type const VectorFieldMem::operator() (
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
inline element_type const VectorFieldMem::operator() (
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

_Deleted: use deepcopy instead._ 
```C++
VectorFieldMem & VectorFieldMem::operator= (
    VectorFieldMem const & other
) = delete
```




<hr>



### function operator= 

```C++
VectorFieldMem & VectorFieldMem::operator= (
    VectorFieldMem && other
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
inline auto VectorFieldMem::operator[] (
    Idx< QueryDDims... > const & slice_spec
) const
```



Get the [**VectorFieldMem**](classVectorFieldMem.md) on the reduced index range which is obtained by indexing the dimensions QueryDDims at the position slice\_spec.




**Parameters:**


* `slice_spec` The slice describing the index range of interest.



**Returns:**

A constant reference to the vector field on the sliced index range. 





        

<hr>



### function operator[] 

_Slice out some dimensions._ 
```C++
template<class... QueryDDims>
inline auto VectorFieldMem::operator[] (
    Idx< QueryDDims... > const & slice_spec
) 
```



Get the [**VectorFieldMem**](classVectorFieldMem.md) on the reduced index range which is obtained by indexing the dimensions QueryDDims at the position slice\_spec.




**Parameters:**


* `slice_spec` The slice describing the index range of interest.



**Returns:**

A modifiable reference to the vector field on the sliced index range. 





        

<hr>



### function operator[] 

_Slice out some dimensions._ 
```C++
template<class... QueryDDims>
inline auto VectorFieldMem::operator[] (
    IdxRange< QueryDDims... > const & oidx_range
) const
```



Get the [**VectorFieldMem**](classVectorFieldMem.md) on the reduced index range passed as an argument.




**Parameters:**


* `oidx_range` The index range of interest.



**Returns:**

A modifiable reference to the vector field on the sliced index range. 





        

<hr>



### function operator[] 

_Slice out some dimensions._ 
```C++
template<class... QueryDDims>
inline auto VectorFieldMem::operator[] (
    IdxRange< QueryDDims... > const & oidx_range
) 
```



Get the [**VectorFieldMem**](classVectorFieldMem.md) on the reduced index range passed as an argument.




**Parameters:**


* `oidx_range` The index range of interest.



**Returns:**

A modifiable reference to the vector field on the sliced index range. 





        

<hr>



### function span\_cview 

```C++
inline view_type VectorFieldMem::span_cview () const
```



Get a constant reference to this vector field.


This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_const\_field should be used.




**Returns:**

A constant reference to this vector field. 





        

<hr>



### function span\_view [1/2]

```C++
inline view_type VectorFieldMem::span_view () const
```



Get a constant reference to this vector field.


This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_field should be used.




**Returns:**

A constant reference to this vector field. 





        

<hr>



### function span\_view [2/2]

```C++
inline span_type VectorFieldMem::span_view () 
```



Get a modifiable reference to this vector field.


This function is designed to match the equivalent function in DDC. In Gysela it should not be called directly. Instead the global function get\_field should be used.




**Returns:**

A modifiable reference to this vector field. 





        

<hr>



### function ~VectorFieldMem 

```C++
VectorFieldMem::~VectorFieldMem () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/vector_field_mem.hpp`

