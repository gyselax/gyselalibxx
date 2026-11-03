

# Class VectorFieldCommon&lt; FieldType, ddc::TypeSeq&lt; DDims... &gt; &gt;

**template &lt;class FieldType, class... DDims&gt;**



[**ClassList**](annotated.md) **>** [**VectorFieldCommon&lt; FieldType, ddc::TypeSeq&lt; DDims... &gt; &gt;**](classVectorFieldCommon_3_01FieldType_00_01ddc_1_1TypeSeq_3_01DDims_8_8_8_01_4_01_4.md)



_A base class containing the functionalities common to_ [_**VectorField**_](classVectorField.md) _and_[_**VectorFieldMem**_](classVectorFieldMem.md) _._[More...](#detailed-description)

* `#include <vector_field_common.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::TypeSeq&lt; DDims... &gt; | [**NDTypeTag**](#typedef-ndtypetag)  <br>_The tags describing the dimensions which make up the elements of the vector._  |
| typedef typename FieldType::span\_type | [**chunk\_span\_type**](#typedef-chunk_span_type)  <br>_The type of a modifiable span of one of these field. This is a DDC keyword used to make this class interchangeable with Field._  |
| typedef typename FieldType::view\_type | [**chunk\_view\_type**](#typedef-chunk_view_type)  <br>_The type of a constant view of one of these field. This is a DDC keyword used to make this class interchangeable with Field._  |
| typedef typename FieldType::discrete\_domain\_type | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The IdxRange on which the fields in this object are defined._  |
| typedef Vector&lt; data\_type &, DDims... &gt; | [**element\_ref\_type**](#typedef-element_ref_type)  <br>_The type of the element returned when indexing this field (e.g. a tuple of doubles)._  |
| typedef Vector&lt; std::remove\_const\_t&lt; data\_type &gt;, DDims... &gt; | [**element\_type**](#typedef-element_type)  <br>_The type of the elements in the fields._  |
| typedef [**discrete\_domain\_type**](classVectorFieldCommon_3_01FieldType_00_01ddc_1_1TypeSeq_3_01DDims_8_8_8_01_4_01_4.md#typedef-discrete_domain_type) | [**index\_range\_type**](#typedef-index_range_type)  <br>_The IdxRange on which the fields in this object are defined._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**deepcopy**](#function-deepcopy) (FieldSrc const & src) <br>_Copy the values of another vector field into this vector field._  |
|  [**index\_range\_type**](classVectorFieldCommon_3_01FieldType_00_01ddc_1_1TypeSeq_3_01DDims_8_8_8_01_4_01_4.md#typedef-index_range_type) | [**idx\_range**](#function-idx_range-12) () noexcept const<br>_Provide access to the index range on which this field is defined._  |
|  IdxRange&lt; QueryDDims... &gt; | [**idx\_range**](#function-idx_range-22) () noexcept const<br>_Provide access to the sub-set of the index range on which this field is defined._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  int | [**rank**](#function-rank) () noexcept<br>_Get the number of dimensions of the index range on which the field is defined._  |






## Protected Attributes

| Type | Name |
| ---: | :--- |
|  std::array&lt; FieldType, [**NDims**](classVectorFieldCommon_3_01FieldType_00_01ddc_1_1TypeSeq_3_01DDims_8_8_8_01_4_01_4.md#variable-ndims) &gt; | [**m\_values**](#variable-m_values)  <br>_The scalar fields storing each component of the vector field._  |


## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**NDims**](#variable-ndims)   = `sizeof...(DDims)`<br>_The number of components of the vector._  |














## Protected Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**VectorFieldCommon**](#function-vectorfieldcommon-13) () = default<br>_Empty Field._  |
|   | [**VectorFieldCommon**](#function-vectorfieldcommon-23) (FieldTypes &&... fields) <br>_Construct a Field on an index range with uninitialized values._  |
|  KOKKOS\_FUNCTION | [**VectorFieldCommon**](#function-vectorfieldcommon-33) (FieldTypes &&... fields) <br>_Construct a Field on an index range with uninitialized values._  |




## Detailed Description


The vector field is stored as one scalar field per component.




**Template parameters:**


* `FieldType` The type of the scalar field used to store each component. 
* `DDims` The tags describing the dimensions which make up the elements of the vector. 




    
## Public Types Documentation




### typedef NDTypeTag 

_The tags describing the dimensions which make up the elements of the vector._ 
```C++
using VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::NDTypeTag =  ddc::TypeSeq<DDims...>;
```




<hr>



### typedef chunk\_span\_type 

_The type of a modifiable span of one of these field. This is a DDC keyword used to make this class interchangeable with Field._ 
```C++
using VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::chunk_span_type =  typename FieldType::span_type;
```




<hr>



### typedef chunk\_view\_type 

_The type of a constant view of one of these field. This is a DDC keyword used to make this class interchangeable with Field._ 
```C++
using VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::chunk_view_type =  typename FieldType::view_type;
```




<hr>



### typedef discrete\_domain\_type 

_The IdxRange on which the fields in this object are defined._ 
```C++
using VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::discrete_domain_type =  typename FieldType::discrete_domain_type;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC IdxRange types are referred to as DiscreteDomain types. 


        

<hr>



### typedef element\_ref\_type 

_The type of the element returned when indexing this field (e.g. a tuple of doubles)._ 
```C++
using VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::element_ref_type =  Vector<data_type&, DDims...>;
```




<hr>



### typedef element\_type 

_The type of the elements in the fields._ 
```C++
using VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::element_type =  Vector<std::remove_const_t<data_type>, DDims...>;
```




<hr>



### typedef index\_range\_type 

_The IdxRange on which the fields in this object are defined._ 
```C++
using VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::index_range_type =  discrete_domain_type;
```




<hr>
## Public Functions Documentation




### function deepcopy 

_Copy the values of another vector field into this vector field._ 
```C++
template<class FieldSrc>
inline void VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::deepcopy (
    FieldSrc const & src
) 
```





**Parameters:**


* `src` The vector field whose values should be copied. 




        

<hr>



### function idx\_range [1/2]

_Provide access to the index range on which this field is defined._ 
```C++
inline index_range_type VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::idx_range () noexcept const
```





**Returns:**

The index range on which this field is defined. 





        

<hr>



### function idx\_range [2/2]

_Provide access to the sub-set of the index range on which this field is defined._ 
```C++
template<class... QueryDDims>
inline IdxRange< QueryDDims... > VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::idx_range () noexcept const
```





**Template parameters:**


* `QueryDDims` The dimensions of interest. 



**Returns:**

The index range on which this field is defined. 





        

<hr>
## Public Static Functions Documentation




### function rank 

_Get the number of dimensions of the index range on which the field is defined._ 
```C++
static inline int VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::rank () noexcept
```





**Returns:**

The rank of the index range. 





        

<hr>
## Protected Attributes Documentation




### variable m\_values 

_The scalar fields storing each component of the vector field._ 
```C++
std::array<FieldType, NDims> VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::m_values;
```




<hr>
## Protected Static Attributes Documentation




### variable NDims 

_The number of components of the vector._ 
```C++
std::size_t VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::NDims;
```




<hr>
## Protected Functions Documentation




### function VectorFieldCommon [1/3]

_Empty Field._ 
```C++
KOKKOS_DEFAULTED_FUNCTION VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::VectorFieldCommon () = default
```




<hr>



### function VectorFieldCommon [2/3]

_Construct a Field on an index range with uninitialized values._ 
```C++
template<class... FieldTypes, class, std::enable_if_t< !std::conjunction_v< std::bool_constant< ddc::is_borrowed_chunk_v< FieldTypes > >... >, int >>
inline explicit VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::VectorFieldCommon (
    FieldTypes &&... fields
) 
```




<hr>



### function VectorFieldCommon [3/3]

_Construct a Field on an index range with uninitialized values._ 
```C++
template<class... FieldTypes, class, std::enable_if_t< std::conjunction_v< std::bool_constant< ddc::is_borrowed_chunk_v< FieldTypes > >... >, int >>
inline explicit KOKKOS_FUNCTION VectorFieldCommon< FieldType, ddc::TypeSeq< DDims... > >::VectorFieldCommon (
    FieldTypes &&... fields
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/vector_field_common.hpp`

