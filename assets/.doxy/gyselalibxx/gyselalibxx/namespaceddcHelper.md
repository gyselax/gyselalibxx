

# Namespace ddcHelper



[**Namespace List**](namespaces.md) **>** [**ddcHelper**](namespaceddcHelper.md)




















## Classes

| Type | Name |
| ---: | :--- |
| class | [**NonUniformInterpolationPoints**](classddcHelper_1_1NonUniformInterpolationPoints.md) &lt;class BSplines, BcXmin, BcXmax&gt;<br>_Helper class for the initialisation of the mesh of interpolation points._  |
| struct | [**is\_non\_uniform\_interpolation\_points**](structddcHelper_1_1is__non__uniform__interpolation__points.md) &lt;class [**T**](structT.md)&gt;<br> |
| struct | [**is\_non\_uniform\_interpolation\_points&lt; NonUniformInterpolationPoints&lt; BSplines, BcXmin, BcXmax &gt; &gt;**](structddcHelper_1_1is__non__uniform__interpolation__points_3_01NonUniformInterpolationPoints_3_047d1c8570873e3c052e2e394afcf9270.md) &lt;class BSplines, BcXmin, BcXmax&gt;<br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename detail::ApplyTemplateToTypeSeq&lt; Templ, TypeSeq &gt;::type | [**apply\_template\_to\_type\_seq\_t**](#typedef-apply_template_to_type_seq_t)  <br>_A helper to get a type sequence by applying a template to a sequence of type tags._  |
| typedef typename detail::TypeSeqIntersection&lt; TypeSeq1, TypeSeq2, ddc::detail::TypeSeq&lt;&gt; &gt;::type | [**type\_seq\_intersection\_t**](#typedef-type_seq_intersection_t)  <br>_A helper to find all types which are found in both TypeSeq1 and TypeSeq2._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**is\_non\_uniform\_interpolation\_points\_v**](#variable-is_non_uniform_interpolation_points_v)   = `[**is\_non\_uniform\_interpolation\_points**](structddcHelper_1_1is__non__uniform__interpolation__points.md)&lt;[**T**](structT.md)&gt;::value`<br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION void | [**assign\_vector\_field\_element**](#function-assign_vector_field_element) ([**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, VectorIndexSet&lt; Dims... &gt;, MemorySpace, LayoutStridedPolicy &gt; field, typename IdxRangeType::discrete\_element\_type idx, [**Vector**](classTensor.md)&lt; ElementType, Dims... &gt; vector) <br>_Copy the elements of a vector into a vector field at a given index._  |
|  auto | [**create\_mirror\_view\_and\_copy**](#function-create_mirror_view_and_copy) (ExecSpace exec\_space, [**VectorField**](classVectorField.md)&lt; ElementType, IdxRangeType, VectorIndexSet&lt; Dims... &gt;, MemorySpace, LayoutStridedPolicy &gt; field) <br> |
|  auto | [**create\_transpose\_mirror**](#function-create_transpose_mirror) (ExecSpace const & execution\_space, Field&lt; ElementType, Domain, MemSpace, FieldLayoutType &gt; src) <br>_Create a data object in the requested dimension ordering using as allocations as possible. This function does not copy data._  |
|  auto | [**create\_transpose\_mirror\_view\_and\_copy**](#function-create_transpose_mirror_view_and_copy) (ExecSpace const & execution\_space, Field&lt; ElementType, Domain, MemSpace, FieldLayoutType &gt; src) <br>_If necessary transpose data into the requested dimension ordering._  |
|  auto | [**deepcopy**](#function-deepcopy) (FieldDst && dst, FieldSrc && src) <br>_Copy the contents of one_ [_**DerivField**_](classDerivField.md) _into another._ |
|  auto | [**deepcopy**](#function-deepcopy) (ExecSpace const & execution\_space, FieldDst && dst, FieldSrc && src) <br>_Copy the contents of one_ [_**DerivField**_](classDerivField.md) _into another._ |
|  void | [**deepcopy**](#function-deepcopy) ([**MultipatchField**](classMultipatchField.md)&lt; T1, Patches... &gt; dst, [**MultipatchField**](classMultipatchField.md)&lt; T2, Patches... &gt; src) <br>_Copy the data from one_ [_**MultipatchField**_](classMultipatchField.md) _into another._ |
|  void | [**dump\_coordinates**](#function-dump_coordinates) (ExecSpace exec\_space, DField&lt; IdxRange&lt; Grid1D &gt;, Layout, MemorySpace &gt; dump\_coord) <br>_Dump the coordinates of a field into the field._  |
|  void | [**dump\_coordinates**](#function-dump_coordinates) (ExecSpace exec\_space, Field&lt; Coord&lt; typename Grid1D::continuous\_dimension\_type &gt;, IdxRange&lt; Grid1D &gt;, Layout, MemorySpace &gt; dump\_coord) <br>_Dump the coordinates of a field into the field._  |
|  KOKKOS\_INLINE\_FUNCTION ElementType & | [**get**](#function-get) ([**Tensor**](classTensor.md)&lt; ElementType, ValidIndexSet... &gt; & tensor) <br>_A helper function to get a modifiable reference to an element of the tensor._  |
|  KOKKOS\_INLINE\_FUNCTION ElementType const & | [**get**](#function-get) ([**Tensor**](classTensor.md)&lt; ElementType, ValidIndexSet... &gt; const & tensor) <br>_A helper function to get an element of the tensor._  |
|  constexpr auto | [**get**](#function-get) (VectorFieldType const & field) noexcept<br> |
|  constexpr auto | [**get**](#function-get) (VectorFieldType & field) noexcept<br> |
|  double | [**maximum\_distance\_between\_adjacent\_points**](#function-maximum_distance_between_adjacent_points) (IdxRange&lt; GridDim &gt; const & idx\_range) <br>_Computes the maximum distance between two adjacent points within an IdxRange._  |
|  KOKKOS\_INLINE\_FUNCTION void | [**restrict\_to\_bspline\_domain**](#function-restrict_to_bspline_domain) (Coord&lt; typename BSpline::continuous\_dimension\_type &gt; & coord) <br>_Calculate the Coordinate inside the domain._  |
|  constexpr std::enable\_if\_t&lt; IDim::continuous\_dimension\_type::PERIODIC, Coord&lt; typename IDim::continuous\_dimension\_type &gt; &gt; | [**restrict\_to\_idx\_range**](#function-restrict_to_idx_range) (Coord&lt; typename IDim::continuous\_dimension\_type &gt; coord, IdxRange&lt; IDim &gt; const & idx\_range) <br>_Calculate the Coordinate inside the domain._  |
|  KOKKOS\_INLINE\_FUNCTION Coord&lt; Dims... &gt; | [**to\_coord**](#function-to_coord) ([**Vector**](classTensor.md)&lt; ElementType, Dims... &gt; const & tensor) <br>_A helper function to convert a vector to a coordinate. This is useful in order to add a Vector to a coordinate to obtain a new coordinate (e.g. when calculating the foot of a characteristic._  |
|  constexpr std::enable\_if\_t&lt;!IDim::continuous\_dimension\_type::PERIODIC, double &gt; | [**total\_interval\_length**](#function-total_interval_length) (IdxRange&lt; IDim &gt; const & idx\_range) <br> |
|  constexpr std::enable\_if\_t&lt; IDim::continuous\_dimension\_type::PERIODIC &&ddc::is\_uniform\_point\_sampling\_v&lt; IDim &gt;, double &gt; | [**total\_interval\_length**](#function-total_interval_length) (IdxRange&lt; IDim &gt; const & idx\_range) <br> |
|  constexpr std::enable\_if\_t&lt; IDim::continuous\_dimension\_type::PERIODIC &&ddc::is\_non\_uniform\_point\_sampling\_v&lt; IDim &gt;, double &gt; | [**total\_interval\_length**](#function-total_interval_length) (IdxRange&lt; IDim &gt; const & idx\_range) <br> |




























## Public Types Documentation




### typedef apply\_template\_to\_type\_seq\_t 

_A helper to get a type sequence by applying a template to a sequence of type tags._ 
```C++
using ddcHelper::apply_template_to_type_seq_t = typedef typename detail::ApplyTemplateToTypeSeq<Templ, TypeSeq>::type;
```




<hr>



### typedef type\_seq\_intersection\_t 

_A helper to find all types which are found in both TypeSeq1 and TypeSeq2._ 
```C++
using ddcHelper::type_seq_intersection_t = typedef typename detail::TypeSeqIntersection<TypeSeq1, TypeSeq2, ddc::detail::TypeSeq<> >::type;
```




<hr>
## Public Attributes Documentation




### variable is\_non\_uniform\_interpolation\_points\_v 

```C++
constexpr bool ddcHelper::is_non_uniform_interpolation_points_v;
```




<hr>
## Public Functions Documentation




### function assign\_vector\_field\_element 

_Copy the elements of a vector into a vector field at a given index._ 
```C++
template<class ElementType, class IdxRangeType, class... Dims, class MemorySpace, class LayoutStridedPolicy>
KOKKOS_INLINE_FUNCTION void ddcHelper::assign_vector_field_element (
    VectorField < ElementType, IdxRangeType, VectorIndexSet< Dims... >, MemorySpace, LayoutStridedPolicy > field,
    typename IdxRangeType::discrete_element_type idx,
    Vector < ElementType, Dims... > vector
) 
```





**Parameters:**


* `field` On output, a [**VectorField**](classVectorField.md) containing the values of vector at the index idx. 
* `idx` An index to specify where the values of the Vector should be copied in the [**VectorField**](classVectorField.md). 
* `vector` A Vector to be copied in the [**VectorField**](classVectorField.md) at the index idx. 




        

<hr>



### function create\_mirror\_view\_and\_copy 

```C++
template<class ExecSpace, class ElementType, class IdxRangeType, class... Dims, class MemorySpace, class LayoutStridedPolicy>
auto ddcHelper::create_mirror_view_and_copy (
    ExecSpace exec_space,
    VectorField < ElementType, IdxRangeType, VectorIndexSet< Dims... >, MemorySpace, LayoutStridedPolicy > field
) 
```




<hr>



### function create\_transpose\_mirror 

_Create a data object in the requested dimension ordering using as allocations as possible. This function does not copy data._ 
```C++
template<class TargetDomain, class ElementType, class Domain, class ExecSpace, class MemSpace, class FieldLayoutType>
auto ddcHelper::create_transpose_mirror (
    ExecSpace const & execution_space,
    Field< ElementType, Domain, MemSpace, FieldLayoutType > src
) 
```





**Parameters:**


* `execution_space` The execution space (Host/Device) where the code will run. 
* `src` The object to be transposed.



**Returns:**

If src is already in the correct dimension ordering, return a view on src. Otherwise return a chunk with the correct dimension ordering. 





        

<hr>



### function create\_transpose\_mirror\_view\_and\_copy 

_If necessary transpose data into the requested dimension ordering._ 
```C++
template<class TargetDomain, class ElementType, class Domain, class ExecSpace, class MemSpace, class FieldLayoutType>
auto ddcHelper::create_transpose_mirror_view_and_copy (
    ExecSpace const & execution_space,
    Field< ElementType, Domain, MemSpace, FieldLayoutType > src
) 
```





**Parameters:**


* `execution_space` The execution space (Host/Device) where the code will run. 
* `src` The object to be transposed.



**Returns:**

If src is already in the correct dimension ordering, return a view on src. Otherwise return a chunk with the correct dimension ordering in which the data from src has been copied. 





        

<hr>



### function deepcopy 

_Copy the contents of one_ [_**DerivField**_](classDerivField.md) _into another._
```C++
template<class FieldDst, class FieldSrc, std::enable_if_t< is_borrowed_deriv_field_v< FieldDst > &&is_borrowed_deriv_field_v< FieldSrc >, bool >>
auto ddcHelper::deepcopy (
    FieldDst && dst,
    FieldSrc && src
) 
```





**Parameters:**


* `dst` The [**DerivField**](classDerivField.md) where the data will be saved. 
* `src` The [**DerivField**](classDerivField.md) whose data will be copied. 




        

<hr>



### function deepcopy 

_Copy the contents of one_ [_**DerivField**_](classDerivField.md) _into another._
```C++
template<class ExecSpace, class FieldDst, class FieldSrc>
auto ddcHelper::deepcopy (
    ExecSpace const & execution_space,
    FieldDst && dst,
    FieldSrc && src
) 
```





**Parameters:**


* `execution_space` The Kokkos execution space where the copy will be carried out. 
* `dst` The [**DerivField**](classDerivField.md) where the data will be saved. 
* `src` The [**DerivField**](classDerivField.md) whose data will be copied. 




        

<hr>



### function deepcopy 

_Copy the data from one_ [_**MultipatchField**_](classMultipatchField.md) _into another._
```C++
template<template< typename P > typename T1, template< typename P > typename T2, class... Patches>
void ddcHelper::deepcopy (
    MultipatchField < T1, Patches... > dst,
    MultipatchField < T2, Patches... > src
) 
```





**Parameters:**


* `dst` The [**MultipatchField**](classMultipatchField.md) that the data will be copied to. 
* `src` The [**MultipatchField**](classMultipatchField.md) that the data will be copied from. 




        

<hr>



### function dump\_coordinates 

_Dump the coordinates of a field into the field._ 
```C++
template<class ExecSpace, class Grid1D, class Layout, class MemorySpace>
inline void ddcHelper::dump_coordinates (
    ExecSpace exec_space,
    DField< IdxRange< Grid1D >, Layout, MemorySpace > dump_coord
) 
```





**Parameters:**


* `exec_space` The execution space on which the code will run. 
* `dump_coord` The field which will contain the coordinates. 




        

<hr>



### function dump\_coordinates 

_Dump the coordinates of a field into the field._ 
```C++
template<class ExecSpace, class Grid1D, class Layout, class MemorySpace>
inline void ddcHelper::dump_coordinates (
    ExecSpace exec_space,
    Field< Coord< typename Grid1D::continuous_dimension_type >, IdxRange< Grid1D >, Layout, MemorySpace > dump_coord
) 
```





**Parameters:**


* `exec_space` The execution space on which the code will run. 
* `dump_coord` The field which will contain the coordinates. 




        

<hr>



### function get 

_A helper function to get a modifiable reference to an element of the tensor._ 
```C++
template<class... QueryIndexTag, class ElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType & ddcHelper::get (
    Tensor < ElementType, ValidIndexSet... > & tensor
) 
```





**Template parameters:**


* `QueryIndexTag` A type describing the relevant index. 



**Parameters:**


* `tensor` The tensor whose elements are examined. 



**Returns:**

The relevant element of the tensor. 





        

<hr>



### function get 

_A helper function to get an element of the tensor._ 
```C++
template<class... QueryIndexTag, class ElementType, class... ValidIndexSet>
KOKKOS_INLINE_FUNCTION ElementType const & ddcHelper::get (
    Tensor < ElementType, ValidIndexSet... > const & tensor
) 
```





**Template parameters:**


* `QueryIndexTag` A type describing the relevant index. 



**Parameters:**


* `tensor` The tensor whose elements are examined. 



**Returns:**

The relevant element of the tensor. 





        

<hr>



### function get 

```C++
template<class QueryTag, class VectorFieldType>
inline constexpr auto ddcHelper::get (
    VectorFieldType const & field
) noexcept
```




<hr>



### function get 

```C++
template<class QueryTag, class VectorFieldType>
inline constexpr auto ddcHelper::get (
    VectorFieldType & field
) noexcept
```




<hr>



### function maximum\_distance\_between\_adjacent\_points 

_Computes the maximum distance between two adjacent points within an IdxRange._ 
```C++
template<class GridDim>
double ddcHelper::maximum_distance_between_adjacent_points (
    IdxRange< GridDim > const & idx_range
) 
```





**Parameters:**


* `idx_range` The domain on which the distance should be calculated.



**Returns:**

The maximum distance between two adjacent points. 





        

<hr>



### function restrict\_to\_bspline\_domain 

_Calculate the Coordinate inside the domain._ 
```C++
template<class BSpline>
KOKKOS_INLINE_FUNCTION void ddcHelper::restrict_to_bspline_domain (
    Coord< typename BSpline::continuous_dimension_type > & coord
) 
```



In the case of a periodic domain, a Coordinate can sometimes be found outside the domain. In this case it is useful to be able to find the equivalent coordinate inside the domain. This function makes that possible.




**Parameters:**


* `coord` The 1D coordinate we want to compute inside the domain.



**Returns:**

The equivalent coordinate inside the domain. 





        

<hr>



### function restrict\_to\_idx\_range 

_Calculate the Coordinate inside the domain._ 
```C++
template<class IDim>
constexpr std::enable_if_t< IDim::continuous_dimension_type::PERIODIC, Coord< typename IDim::continuous_dimension_type > > ddcHelper::restrict_to_idx_range (
    Coord< typename IDim::continuous_dimension_type > coord,
    IdxRange< IDim > const & idx_range
) 
```



In the case of a periodic domain, a Coordinate can sometimes be found outside the domain. In this case it is useful to be able to find the equivalent coordinate inside the domain. This function makes that possible.




**Parameters:**


* `coord` The 1D coordinate we want to compute inside the domain. 
* `idx_range` The domain where the coordinate is defined.



**Returns:**

The equivalent coordinate inside the domain. 





        

<hr>



### function to\_coord 

_A helper function to convert a vector to a coordinate. This is useful in order to add a Vector to a coordinate to obtain a new coordinate (e.g. when calculating the foot of a characteristic._ 
```C++
template<class ElementType, class... Dims>
KOKKOS_INLINE_FUNCTION Coord< Dims... > ddcHelper::to_coord (
    Vector < ElementType, Dims... > const & tensor
) 
```





**Parameters:**


* `tensor` The tensor to be converted. 



**Returns:**

The new coordinate. 





        

<hr>



### function total\_interval\_length 

```C++
template<class IDim>
constexpr std::enable_if_t<!IDim::continuous_dimension_type::PERIODIC, double > ddcHelper::total_interval_length (
    IdxRange< IDim > const & idx_range
) 
```



Calculate the total length of a non-periodic domain.




**Parameters:**


* `idx_range` The domain on which the length should be calculated.



**Returns:**

The length of the domain. 





        

<hr>



### function total\_interval\_length 

```C++
template<class IDim>
constexpr std::enable_if_t< IDim::continuous_dimension_type::PERIODIC &&ddc::is_uniform_point_sampling_v< IDim >, double > ddcHelper::total_interval_length (
    IdxRange< IDim > const & idx_range
) 
```



Calculate the total length of a uniform periodic domain.




**Parameters:**


* `idx_range` The domain on which the length should be calculated.



**Returns:**

The length of the domain. 





        

<hr>



### function total\_interval\_length 

```C++
template<class IDim>
constexpr std::enable_if_t< IDim::continuous_dimension_type::PERIODIC &&ddc::is_non_uniform_point_sampling_v< IDim >, double > ddcHelper::total_interval_length (
    IdxRange< IDim > const & idx_range
) 
```



Calculate the total length of a non-uniform periodic domain.




**Parameters:**


* `idx_range` The domain on which the length should be calculated.



**Returns:**

The length of the domain. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/derivative_field.hpp`

