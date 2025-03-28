

# Class DerivFieldCommon&lt; FieldType, IdxRange&lt; DDims... &gt; &gt;

**template &lt;class FieldType, class... DDims&gt;**



[**ClassList**](annotated.md) **>** [**DerivFieldCommon&lt; FieldType, IdxRange&lt; DDims... &gt; &gt;**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md)



_An abstract class which holds a chunk of memory describing a field and its derivatives. This is the superclass for_ [_**DerivFieldMem**_](classDerivFieldMem.md) _and_[_**DerivField**_](classDerivField.md) _._[More...](#detailed-description)

* `#include <derivative_field_common.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef FieldType | [**chunk\_type**](#typedef-chunk_type)  <br>_The type of the field stored in the array._  |
| typedef detail::deriv\_sub\_set\_t&lt; ddc::detail::TypeSeq&lt; DDims... &gt; &gt; | [**deriv\_tags**](#typedef-deriv_tags)  <br>_A type sequence containing all derivatives present in this object._  |
| typedef IdxRange&lt; DDims... &gt; | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The IdxRange on which the fields in this object are defined._  |
| typedef Idx&lt; DDims... &gt; | [**discrete\_element\_type**](#typedef-discrete_element_type)  <br>_The Idx which can be used to index this object._  |
| typedef typename chunk\_type::element\_type | [**element\_type**](#typedef-element_type)  <br>_The type of the elements in the fields._  |
| typedef [**discrete\_domain\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-discrete_domain_type) | [**index\_range\_type**](#typedef-index_range_type)  <br>_The IdxRange on which the fields in this object are defined._  |
| typedef [**discrete\_element\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-discrete_element_type) | [**index\_type**](#typedef-index_type)  <br>_The Idx which can be used to index this object._  |
| typedef typename detail::strip\_deriv\_t&lt; [**deriv\_tags**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-deriv_tags) &gt; | [**physical\_deriv\_grids**](#typedef-physical_deriv_grids)  <br>_A type sequence containing all physical dimensions for which derivatives are present in this object._  |
| typedef ddc::type\_seq\_remove\_t&lt; ddc::detail::TypeSeq&lt; DDims... &gt;, [**deriv\_tags**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-deriv_tags) &gt; | [**physical\_grids**](#typedef-physical_grids)  <br>_A type sequence containing all the physical dimensions on which the fields are defined._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  auto | [**get\_mdspan**](#function-get_mdspan-12) (IdxRange&lt; ODims... &gt; provided\_deriv\_idx\_range) <br>_Get one of the mdspans from the internal array internal\_fields. This function takes index ranges on the derivative directions. Where derivatives are missing it is assumed that the 0-th order derivative is requested. This dimension is stripped from the resulting field. This is the recommended way to access the internal fields._  |
|  auto | [**get\_mdspan**](#function-get_mdspan-22) () <br>_Get the mdspan holding the values of the function from the internal array internal\_fields._  |
|  auto | [**get\_values\_field**](#function-get_values_field-12) () <br>_Get the Field which holds the values of the function._  |
|  auto | [**get\_values\_field**](#function-get_values_field-22) () const<br>_Get the Field which holds the values of the function._  |
|  constexpr auto | [**operator[]**](#function-operator) (Idx&lt; QueryDDims... &gt; const & slice\_spec) const<br>_Get a ConstField describing a subset of the data._  |
|  constexpr auto | [**operator[]**](#function-operator_1) (Idx&lt; QueryDDims... &gt; const & slice\_spec) <br>_Get a Field describing a subset of the data._  |
|  KOKKOS\_FUNCTION constexpr auto | [**operator[]**](#function-operator_2) (IdxRange&lt; QueryDDims... &gt; const & oidx\_range) <br>_Get a Field describing a subset of the data. This function allows a slice to be obtained however it is designed to return a Field. It is therefore not possible to request data from multiple fields (e.g. derivatives from 0 to 3)._  |
|  KOKKOS\_FUNCTION constexpr auto | [**operator[]**](#function-operator_3) (IdxRange&lt; QueryDDims... &gt; const & oidx\_range) const<br>_Get a ConstField describing a subset of the data. This function allows a slice to be obtained however it is designed to return a ConstField. It is therefore not possible to request data from multiple fields (e.g. derivatives from 0 to 3)._  |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~DerivFieldCommon**](#function-derivfieldcommon) () = default<br> |




## Protected Types

| Type | Name |
| ---: | :--- |
| typedef typename FieldType::span\_type | [**chunk\_span**](#typedef-chunk_span)  <br>_The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field._  |
| typedef typename FieldType::view\_type | [**chunk\_view**](#typedef-chunk_view)  <br>_The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field._  |
| typedef typename ddc::detail::convert\_type\_seq\_to\_discrete\_domain\_t&lt; [**deriv\_tags**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-deriv_tags) &gt; | [**discrete\_deriv\_idx\_range\_type**](#typedef-discrete_deriv_idx_range_type)  <br>_The IdxRange which describes the derivatives present on each field._  |
| typedef typename discrete\_deriv\_idx\_range\_type::discrete\_element\_type | [**discrete\_deriv\_index\_type**](#typedef-discrete_deriv_index_type)  <br>_The Idx which describes the order of the derivatives in each dimension. (e.g. second-order derivative)._  |
| typedef typename discrete\_deriv\_idx\_range\_type::discrete\_vector\_type | [**discrete\_deriv\_vector\_type**](#typedef-discrete_deriv_vector_type)  <br>_The Idx which describes the order of the derivatives in each dimension. (e.g. second-order derivative)._  |
| typedef Kokkos::mdspan&lt; [**element\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-element_type), Kokkos::dextents&lt; std::size\_t, sizeof...(DDims)&gt;, Kokkos::layout\_stride &gt; | [**internal\_mdspan\_type**](#typedef-internal_mdspan_type)  <br>_The type of the memory block stored in the array internal\_fields._  |
| typedef Kokkos::mdspan&lt; const [**element\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-element_type), Kokkos::dextents&lt; std::size\_t, sizeof...(DDims)&gt;, Kokkos::layout\_stride &gt; | [**internal\_mdview\_type**](#typedef-internal_mdview_type)  <br>_The type of a constant view on the memory block stored in the array internal\_fields._  |
| typedef typename ddc::detail::convert\_type\_seq\_to\_discrete\_domain\_t&lt; [**physical\_grids**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-physical_grids) &gt; | [**physical\_idx\_range\_type**](#typedef-physical_idx_range_type)  <br>_The index range for the field excluding derivatives._  |
| typedef typename physical\_idx\_range\_type::discrete\_element\_type | [**physical\_index\_type**](#typedef-physical_index_type)  <br>_The Idx which describes the physical position where values are defined._  |




## Protected Attributes

| Type | Name |
| ---: | :--- |
|  std::array&lt; [**internal\_mdspan\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-internal_mdspan_type), [**n\_fields**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#variable-n_fields) &gt; | [**internal\_fields**](#variable-internal_fields)  <br>_The internal fields describing the values and derivatives._  |
|  to\_subidx\_range\_collection&lt; [**physical\_deriv\_grids**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-physical_deriv_grids) &gt; | [**m\_cross\_derivative\_idx\_range**](#variable-m_cross_derivative_idx_range)  <br>_The physical index ranges on which the derivatives are defined._  |
|  [**discrete\_deriv\_idx\_range\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-discrete_deriv_idx_range_type) | [**m\_deriv\_idx\_range**](#variable-m_deriv_idx_range)  <br>_The index range of available derivatives._  |
|  [**physical\_idx\_range\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-physical_idx_range_type) | [**m\_physical\_idx\_range**](#variable-m_physical_idx_range)  <br>_The physical index range on which the values are defined._  |


## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr int | [**n\_fields**](#variable-n_fields)   = `1 &lt;&lt; ddc::type\_seq\_size\_v&lt;[**deriv\_tags**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-deriv_tags)&gt;`<br>_The number of fields which must be created to describe this object._  |














## Protected Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**DerivFieldCommon**](#function-derivfieldcommon) ([**physical\_idx\_range\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-physical_idx_range_type) physical\_idx\_range, [**discrete\_deriv\_idx\_range\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-discrete_deriv_idx_range_type) deriv\_idx\_range, to\_subidx\_range\_collection&lt; [**physical\_deriv\_grids**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-physical_deriv_grids) &gt; cross\_derivative\_idx\_range) <br>_Protected constructor to be used by subclasses to initialise index ranges._  |
|  KOKKOS\_FUNCTION int | [**get\_array\_index**](#function-get_array_index) (Idx&lt; Tag... &gt; idx) const<br>_An internal function which provides the index of a field inside the internal\_fields array. An Idx describes the derivatives of interest. n-th order derivatives are stored in the same field for all n!=0 so it is sufficient to provide any valid element from the derivatives._  |
|  KOKKOS\_FUNCTION std::pair&lt; int, [**index\_type**](classDerivFieldCommon_3_01FieldType_00_01IdxRange_3_01DDims_8_8_8_01_4_01_4.md#typedef-index_type) &gt; | [**get\_index**](#function-get_index) (DElem elem) const<br>_An internal function which provides the index of an element inside the internal\_fields. An Idx describes the element of interest. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested._  |
|  KOKKOS\_FUNCTION auto | [**get\_internal\_field**](#function-get_internal_field-12) (IdxRange&lt; ODims... &gt; idx\_range) const<br> |
|  KOKKOS\_FUNCTION auto | [**get\_internal\_field**](#function-get_internal_field-22) (Idx&lt; ODims... &gt; elem) const<br> |
|  KOKKOS\_FUNCTION constexpr auto | [**get\_slicer\_for**](#function-get_slicer_for-12) (Idx&lt; ODDims... &gt; const & slice\_idx, int array\_idx) const<br>_Get an object which can be used to slice an mdspan._  |
|  KOKKOS\_FUNCTION constexpr auto | [**get\_slicer\_for**](#function-get_slicer_for-22) (IdxRange&lt; ODDims... &gt; const & slice\_idx\_range, int array\_idx) const<br>_Get an object which can be used to slice an mdspan._  |




## Detailed Description




**Template parameters:**


* `FieldType` The type of the object stored in the internal\_fields array. 
* `IdxRange<DDims...>` The index range on which the internal fields are defined. This index range is the physical index range on which the values are defined combined with the index range of the derivatives of interest (e.g. IdxRange&lt;Deriv&lt;IDimX&gt;, IDimX, IDimY&gt;). 




    
## Public Types Documentation




### typedef chunk\_type 

_The type of the field stored in the array._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::chunk_type =  FieldType;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC FieldMem types are referred to as Chunk types and Field types are referred to as ChunkSpan/ChunkView. 


        

<hr>



### typedef deriv\_tags 

_A type sequence containing all derivatives present in this object._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::deriv_tags =  detail::deriv_sub_set_t<ddc::detail::TypeSeq<DDims...> >;
```




<hr>



### typedef discrete\_domain\_type 

_The IdxRange on which the fields in this object are defined._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::discrete_domain_type =  IdxRange<DDims...>;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC IdxRange types are referred to as DiscreteDomain types. 


        

<hr>



### typedef discrete\_element\_type 

_The Idx which can be used to index this object._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::discrete_element_type =  Idx<DDims...>;
```



This is a DDC keyword used to make this class interchangeable with Field. In DDC Idx types are referred to as DiscreteElement types. 


        

<hr>



### typedef element\_type 

_The type of the elements in the fields._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::element_type =  typename chunk_type::element_type;
```




<hr>



### typedef index\_range\_type 

_The IdxRange on which the fields in this object are defined._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::index_range_type =  discrete_domain_type;
```




<hr>



### typedef index\_type 

_The Idx which can be used to index this object._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::index_type =  discrete_element_type;
```




<hr>



### typedef physical\_deriv\_grids 

_A type sequence containing all physical dimensions for which derivatives are present in this object._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::physical_deriv_grids =  typename detail::strip_deriv_t<deriv_tags>;
```




<hr>



### typedef physical\_grids 

_A type sequence containing all the physical dimensions on which the fields are defined._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::physical_grids =  ddc::type_seq_remove_t<ddc::detail::TypeSeq<DDims...>, deriv_tags>;
```




<hr>
## Public Functions Documentation




### function get\_mdspan [1/2]

_Get one of the mdspans from the internal array internal\_fields. This function takes index ranges on the derivative directions. Where derivatives are missing it is assumed that the 0-th order derivative is requested. This dimension is stripped from the resulting field. This is the recommended way to access the internal fields._ 
```C++
template<class... ODims>
inline auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_mdspan (
    IdxRange< ODims... > provided_deriv_idx_range
) 
```





**Parameters:**


* `provided_deriv_idx_range` The derivative index range which should be retained.



**Returns:**

Field The field on the physical index range and the requested index ranges. 





        

<hr>



### function get\_mdspan [2/2]

_Get the mdspan holding the values of the function from the internal array internal\_fields._ 
```C++
inline auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_mdspan () 
```





**Returns:**

Field The field on the physical index range and the requested index ranges. 





        

<hr>



### function get\_values\_field [1/2]

_Get the Field which holds the values of the function._ 
```C++
inline auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_values_field () 
```



This function is equivalent to calling operator[] with a 0D IdxRange.




**Returns:**

Field The field on the physical index range. 





        

<hr>



### function get\_values\_field [2/2]

_Get the Field which holds the values of the function._ 
```C++
inline auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_values_field () const
```



This function is equivalent to calling operator[] with a 0D IdxRange.




**Returns:**

Field The field on the physical index range. 





        

<hr>



### function operator[] 

_Get a ConstField describing a subset of the data._ 
```C++
template<class... QueryDDims>
inline constexpr auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::operator[] (
    Idx< QueryDDims... > const & slice_spec
) const
```





**Parameters:**


* `slice_spec` A discrete element describing the position at which these dimensions should be indexed. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested.



**Returns:**

ConstField A subset of the data. 





        

<hr>



### function operator[] 

_Get a Field describing a subset of the data._ 
```C++
template<class... QueryDDims>
inline constexpr auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::operator[] (
    Idx< QueryDDims... > const & slice_spec
) 
```





**Parameters:**


* `slice_spec` A discrete element describing the position at which these dimensions should be indexed. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested.



**Returns:**

Field A subset of the data. 





        

<hr>



### function operator[] 

_Get a Field describing a subset of the data. This function allows a slice to be obtained however it is designed to return a Field. It is therefore not possible to request data from multiple fields (e.g. derivatives from 0 to 3)._ 
```C++
template<class... QueryDDims>
inline KOKKOS_FUNCTION constexpr auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::operator[] (
    IdxRange< QueryDDims... > const & oidx_range
) 
```





**Parameters:**


* `oidx_range` A discrete index range describing the position at which these dimensions should be indexed. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested.



**Returns:**

Field A subset of the data. 





        

<hr>



### function operator[] 

_Get a ConstField describing a subset of the data. This function allows a slice to be obtained however it is designed to return a ConstField. It is therefore not possible to request data from multiple fields (e.g. derivatives from 0 to 3)._ 
```C++
template<class... QueryDDims>
inline KOKKOS_FUNCTION constexpr auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::operator[] (
    IdxRange< QueryDDims... > const & oidx_range
) const
```





**Parameters:**


* `oidx_range` A discrete index range describing the position at which these dimensions should be indexed. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested.



**Returns:**

ConstField A subset of the data. 





        

<hr>



### function ~DerivFieldCommon 

```C++
KOKKOS_DEFAULTED_FUNCTION DerivFieldCommon< FieldType, IdxRange< DDims... > >::~DerivFieldCommon () = default
```




<hr>
## Protected Types Documentation




### typedef chunk\_span 

_The type of a modifiable span of this field. This is a DDC keyword used to make this class interchangeable with Field._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::chunk_span =  typename FieldType::span_type;
```




<hr>



### typedef chunk\_view 

_The type of a constant view of this field. This is a DDC keyword used to make this class interchangeable with Field._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::chunk_view =  typename FieldType::view_type;
```




<hr>



### typedef discrete\_deriv\_idx\_range\_type 

_The IdxRange which describes the derivatives present on each field._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::discrete_deriv_idx_range_type =  typename ddc::detail::convert_type_seq_to_discrete_domain_t<deriv_tags>;
```




<hr>



### typedef discrete\_deriv\_index\_type 

_The Idx which describes the order of the derivatives in each dimension. (e.g. second-order derivative)._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::discrete_deriv_index_type =  typename discrete_deriv_idx_range_type::discrete_element_type;
```




<hr>



### typedef discrete\_deriv\_vector\_type 

_The Idx which describes the order of the derivatives in each dimension. (e.g. second-order derivative)._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::discrete_deriv_vector_type =  typename discrete_deriv_idx_range_type::discrete_vector_type;
```




<hr>



### typedef internal\_mdspan\_type 

_The type of the memory block stored in the array internal\_fields._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::internal_mdspan_type =  Kokkos::mdspan< element_type, Kokkos::dextents<std::size_t, sizeof...(DDims)>, Kokkos::layout_stride>;
```




<hr>



### typedef internal\_mdview\_type 

_The type of a constant view on the memory block stored in the array internal\_fields._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::internal_mdview_type =  Kokkos::mdspan< const element_type, Kokkos::dextents<std::size_t, sizeof...(DDims)>, Kokkos::layout_stride>;
```




<hr>



### typedef physical\_idx\_range\_type 

_The index range for the field excluding derivatives._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::physical_idx_range_type =  typename ddc::detail::convert_type_seq_to_discrete_domain_t<physical_grids>;
```




<hr>



### typedef physical\_index\_type 

_The Idx which describes the physical position where values are defined._ 
```C++
using DerivFieldCommon< FieldType, IdxRange< DDims... > >::physical_index_type =  typename physical_idx_range_type::discrete_element_type;
```




<hr>
## Protected Attributes Documentation




### variable internal\_fields 

_The internal fields describing the values and derivatives._ 
```C++
std::array<internal_mdspan_type, n_fields> DerivFieldCommon< FieldType, IdxRange< DDims... > >::internal_fields;
```



The fields which contain the values have different index ranges to the fields containing derivatives so a DDC object cannot be used directly. E.g. for a 2D field ([**X**](structX.md),[**Y**](structY.md)) with derivatives provided in both directions the elements of internal\_fields have the type : DFieldMem&lt;IdxRange&lt;Deriv&lt;IDimX&gt;, Deriv&lt;IDimY&gt;, IDimX, IDimY&gt; The derivative index ranges are then defined such that the elements of internal\_fields represent: 0 :  1 :  2 :  3 :  


        

<hr>



### variable m\_cross\_derivative\_idx\_range 

_The physical index ranges on which the derivatives are defined._ 
```C++
to_subidx_range_collection<physical_deriv_grids> DerivFieldCommon< FieldType, IdxRange< DDims... > >::m_cross_derivative_idx_range;
```




<hr>



### variable m\_deriv\_idx\_range 

_The index range of available derivatives._ 
```C++
discrete_deriv_idx_range_type DerivFieldCommon< FieldType, IdxRange< DDims... > >::m_deriv_idx_range;
```




<hr>



### variable m\_physical\_idx\_range 

_The physical index range on which the values are defined._ 
```C++
physical_idx_range_type DerivFieldCommon< FieldType, IdxRange< DDims... > >::m_physical_idx_range;
```




<hr>
## Protected Static Attributes Documentation




### variable n\_fields 

_The number of fields which must be created to describe this object._ 
```C++
constexpr int DerivFieldCommon< FieldType, IdxRange< DDims... > >::n_fields;
```




<hr>
## Protected Functions Documentation




### function DerivFieldCommon 

_Protected constructor to be used by subclasses to initialise index ranges._ 
```C++
inline KOKKOS_FUNCTION DerivFieldCommon< FieldType, IdxRange< DDims... > >::DerivFieldCommon (
    physical_idx_range_type physical_idx_range,
    discrete_deriv_idx_range_type deriv_idx_range,
    to_subidx_range_collection< physical_deriv_grids > cross_derivative_idx_range
) 
```





**Parameters:**


* `physical_idx_range` The index range on which the values of the function are defined. 
* `deriv_idx_range` The index range of the provided derivatives. 
* `cross_derivative_idx_range` The cross product of the index ranges on which the derivatives of the function are defined. 




        

<hr>



### function get\_array\_index 

_An internal function which provides the index of a field inside the internal\_fields array. An Idx describes the derivatives of interest. n-th order derivatives are stored in the same field for all n!=0 so it is sufficient to provide any valid element from the derivatives._ 
```C++
template<class... Tag>
inline KOKKOS_FUNCTION int DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_array_index (
    Idx< Tag... > idx
) const
```



The index is calculated using a bit mask. Mathematically the equation which determines the index in `internal_fields` is:  where  is the index of the dimension  in the tags.




**Parameters:**


* `idx` The derivatives of interest.



**Returns:**

int The index of the internal field inside the array internal\_fields. 




**Returns:**

discrete\_deriv\_idx\_range\_type The index range of the derivatives at the field at the index. 





        

<hr>



### function get\_index 

_An internal function which provides the index of an element inside the internal\_fields. An Idx describes the element of interest. If information about the derivatives is missing then it is assumed that the 0-th order derivative is requested._ 
```C++
template<class DElem>
inline KOKKOS_FUNCTION std::pair< int, index_type > DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_index (
    DElem elem
) const
```





**Parameters:**


* `elem` The element of interest.



**Returns:**

int The index of the internal field inside the array internal\_fields. 




**Returns:**

index\_type The index of the element of interest inside the field of interest. 





        

<hr>



### function get\_internal\_field [1/2]

```C++
template<class... ODims>
inline KOKKOS_FUNCTION auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_internal_field (
    IdxRange< ODims... > idx_range
) const
```



Get a Field from a subset of one of the mdspans in internal\_fields. The provided IdxRange is used to slice the mdspan in such a way that the resulting mdspan can be saved in a Field. This means that all information about derivatives must be provided. If information about a derivative is missing then it is assumed that the 0-th order derivative is requested.




**Parameters:**


* `idx_range` The index range used to slice the mdspan.



**Returns:**

Field The subset of the internal mdspan. 





        

<hr>



### function get\_internal\_field [2/2]

```C++
template<class... ODims>
inline KOKKOS_FUNCTION auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_internal_field (
    Idx< ODims... > elem
) const
```



Get a Field from a subset of one of the mdspans in internal\_fields. The provided Idx is used to slice the mdspan in such a way that the resulting mdspan can be saved in a Field. This means that all information about derivatives must be provided. If information about a derivative is missing then it is assumed that the 0-th order derivative is requested.




**Parameters:**


* `elem` The element used to slice the mdspan.



**Returns:**

Field The subset of the internal mdspan. 





        

<hr>



### function get\_slicer\_for [1/2]

_Get an object which can be used to slice an mdspan._ 
```C++
template<class QueryDDim, class... ODDims>
inline KOKKOS_FUNCTION constexpr auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_slicer_for (
    Idx< ODDims... > const & slice_idx,
    int array_idx
) const
```





**Parameters:**


* `slice_idx` The DDC element which should be used to slice the field. 
* `array_idx` The index of the mdspan in internal\_fields that will be sliced.



**Template parameters:**


* `QueryDDim` The dimension along which we want to slice.



**Returns:**

An index or a slice which can be used to slice an mdspan. 





        

<hr>



### function get\_slicer\_for [2/2]

_Get an object which can be used to slice an mdspan._ 
```C++
template<class QueryDDim, class... ODDims>
inline KOKKOS_FUNCTION constexpr auto DerivFieldCommon< FieldType, IdxRange< DDims... > >::get_slicer_for (
    IdxRange< ODDims... > const & slice_idx_range,
    int array_idx
) const
```





**Parameters:**


* `slice_idx_range` The DDC index range which should be used to slice the field. 
* `array_idx` The index of the mdspan in internal\_fields that will be sliced.



**Template parameters:**


* `QueryDDim` The dimension along which the we want to slice.



**Returns:**

A slice (often in the form of a (start, end) pair) which can be used to slice an mdspan. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/derivative_field_common.hpp`

