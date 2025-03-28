

# Class MultipatchField

**template &lt;template&lt; typename P &gt; typename T, class... Patches&gt;**



[**ClassList**](annotated.md) **>** [**MultipatchField**](classMultipatchField.md)



_A class to store field objects on patches._ [More...](#detailed-description)

* `#include <multipatch_field.hpp>`



Inherits the following classes: [MultipatchType](classMultipatchType.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**T**](structT.md)&lt; [**Patch**](structPatch.md) &gt;[**::view\_type**](classMultipatchField.md#typedef-view_type) | [**InternalConstFieldOnPatch**](#typedef-internalconstfieldonpatch)  <br>_An internal type alias that is only instantiated if the get\_const\_field method is called._  |
| typedef typename [**T**](structT.md)&lt; [**Patch**](structPatch.md) &gt;[**::span\_type**](classMultipatchField.md#typedef-span_type) | [**InternalFieldOnPatch**](#typedef-internalfieldonpatch)  <br>_An internal type alias that is only instantiated if the get\_const\_field method is called._  |
| typedef typename [**T**](structT.md)&lt; [**Patch**](structPatch.md) &gt;[**::discrete\_domain\_type**](classMultipatchField.md#typedef-discrete_domain_type) | [**InternalIdxRangeOnPatch**](#typedef-internalidxrangeonpatch)  <br>_An internal type alias that is only instantiated if the idx\_range method is called._  |
| typedef [**MultipatchType**](classMultipatchType.md)&lt; [**T**](structT.md), Patches... &gt; | [**base\_type**](#typedef-base_type)  <br>_The_ [_**MultipatchType**_](classMultipatchType.md) _from which this class inherits._ |
| typedef [**MultipatchType**](classMultipatchType.md)&lt; [**InternalIdxRangeOnPatch**](classMultipatchField.md#typedef-internalidxrangeonpatch), Patches... &gt; | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The type of the index ranges that can be used to access this field._  |
| typedef typename base\_type::example\_element::element\_type | [**element\_type**](#typedef-element_type)  <br>_The type of the elements inside the field._  |
| typedef typename [**InternalIdxRangeOnPatch**](classMultipatchField.md#typedef-internalidxrangeonpatch)&lt; [**Patch**](structPatch.md) &gt;::discrete\_element\_type | [**idx\_type**](#typedef-idx_type)  <br>_The type used to index the field on the specified patch._  |
| typedef typename base\_type::example\_element::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The memory space (CPU/GPU) where the data is saved._  |
| typedef [**MultipatchField**](classMultipatchField.md)&lt; [**InternalFieldOnPatch**](classMultipatchField.md#typedef-internalfieldonpatch), Patches... &gt; | [**span\_type**](#typedef-span_type)  <br>_The type of a modifiable reference to this multipatch field._  |
| typedef [**MultipatchField**](classMultipatchField.md)&lt; [**InternalConstFieldOnPatch**](classMultipatchField.md#typedef-internalconstfieldonpatch), Patches... &gt; | [**view\_type**](#typedef-view_type)  <br>_The type of a constant reference to this multipatch field._  |


## Public Types inherited from MultipatchType

See [MultipatchType](classMultipatchType.md)

| Type | Name |
| ---: | :--- |
| typedef ddc::detail::TypeSeq&lt; Patches... &gt; | [**PatchOrdering**](classMultipatchType.md#typedef-patchordering)  <br>_A tag storing the order of Patches in this_ [_**MultipatchType**_](classMultipatchType.md) _._ |
| typedef [**T**](structT.md)&lt; ddc::type\_seq\_element\_t&lt; 0, [**PatchOrdering**](classMultipatchType.md#typedef-patchordering) &gt; &gt; | [**example\_element**](classMultipatchType.md#typedef-example_element)  <br>_The type of one of the elements of the_ [_**MultipatchType**_](classMultipatchType.md) _. This can be used to check that types are as expected using functions such as ddc::is\_chunk\_v._ |






































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**MultipatchField**](#function-multipatchfield-25) ([**T**](structT.md)&lt; Patches &gt;... args) <br> |
|  KOKKOS\_FUNCTION | [**MultipatchField**](#function-multipatchfield-35) (MultipatchObj & other) <br> |
|   | [**MultipatchField**](#function-multipatchfield-45) (MultipatchObj & other) <br> |
|   | [**MultipatchField**](#function-multipatchfield-55) ([**MultipatchField**](classMultipatchField.md)&lt; OtherType, OPatches... &gt; && other) <br> |
|  KOKKOS\_FUNCTION auto | [**get**](#function-get) () const<br> |
|  KOKKOS\_FUNCTION auto | [**get\_const\_field**](#function-get_const_field) () const<br>_Get a_ [_**MultipatchField**_](classMultipatchField.md) _containing constant fields so the values cannot be modified._ |
|  KOKKOS\_FUNCTION auto | [**get\_field**](#function-get_field) () <br>_Get a_ [_**MultipatchField**_](classMultipatchField.md) _containing modifiable fields._ |
|  auto | [**idx\_range**](#function-idx_range) () const<br>_Get a_ [_**MultipatchType**_](classMultipatchType.md) _containing the index ranges on which the fields are defined._ |
|  KOKKOS\_FUNCTION auto | [**span\_cview**](#function-span_cview) () <br>_Get a_ [_**MultipatchField**_](classMultipatchField.md) _containing constant fields so the values cannot be modified. This function matches the DDC name to allow the global get\_const\_field to be defined._ |
|  KOKKOS\_FUNCTION auto | [**span\_view**](#function-span_view) () <br>_Get a_ [_**MultipatchField**_](classMultipatchField.md) _containing modifiable fields. This function matches the DDC name to allow the global get\_const\_field to be defined._ |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~MultipatchField**](#function-multipatchfield) () noexcept<br> |


## Public Functions inherited from MultipatchType

See [MultipatchType](classMultipatchType.md)

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**MultipatchType**](classMultipatchType.md#function-multipatchtype-35) ([**T**](structT.md)&lt; Patches &gt;... args) <br> |
|  KOKKOS\_FUNCTION | [**MultipatchType**](classMultipatchType.md#function-multipatchtype-45) ([**MultipatchType**](classMultipatchType.md)&lt; OtherType, OPatches... &gt; const & other) <br> |
|   | [**MultipatchType**](classMultipatchType.md#function-multipatchtype-55) ([**MultipatchType**](classMultipatchType.md)&lt; OtherType, OPatches... &gt; && other) <br> |
|  KOKKOS\_FUNCTION [**T**](structT.md)&lt; [**Patch**](structPatch.md) &gt; | [**get**](classMultipatchType.md#function-get) () const<br> |
|  KOKKOS\_FUNCTION std::tuple&lt; [**T**](structT.md)&lt; Patches &gt;... &gt; const & | [**get\_tuple**](classMultipatchType.md#function-get_tuple) () const<br>_Get a constant reference to the tuple of objects stored inside this_ [_**MultipatchType**_](classMultipatchType.md) _._ |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~MultipatchType**](classMultipatchType.md#function-multipatchtype) () noexcept<br> |




## Public Static Functions inherited from MultipatchType

See [MultipatchType](classMultipatchType.md)

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**size**](classMultipatchType.md#function-size) () <br>_Get the number of objects stored inside the class. This is equal to the number of patches._  |












## Protected Attributes inherited from MultipatchType

See [MultipatchType](classMultipatchType.md)

| Type | Name |
| ---: | :--- |
|  std::tuple&lt; [**T**](structT.md)&lt; Patches &gt;... &gt; | [**m\_tuple**](classMultipatchType.md#variable-m_tuple)  <br>_The internal tuple containing the data._  |
































## Protected Functions inherited from MultipatchType

See [MultipatchType](classMultipatchType.md)

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**MultipatchType**](classMultipatchType.md#function-multipatchtype-25) (std::tuple&lt; [**T**](structT.md)&lt; Patches &gt;... &gt; && tuple) <br> |






## Detailed Description


On a multipatch domain when we have objects and types defined on different patches, e.g. fields. They can be stored in this class and then be accessed by the patch they are defined on.




**Template parameters:**


* [**T**](structT.md) The type of the fields/derivative fields/vector fields that are stored on the given patches. 
* `Patches` The patches of the objects in the same order of the patches that the given objects are defined on.



**Warning:**

The objects have to be defined on different patches. Otherwise retrieving them by their patch is ill-defined. 





    
## Public Types Documentation




### typedef InternalConstFieldOnPatch 

_An internal type alias that is only instantiated if the get\_const\_field method is called._ 
```C++
using MultipatchField< T, Patches >::InternalConstFieldOnPatch =  typename T<Patch>::view_type;
```




<hr>



### typedef InternalFieldOnPatch 

_An internal type alias that is only instantiated if the get\_const\_field method is called._ 
```C++
using MultipatchField< T, Patches >::InternalFieldOnPatch =  typename T<Patch>::span_type;
```




<hr>



### typedef InternalIdxRangeOnPatch 

_An internal type alias that is only instantiated if the idx\_range method is called._ 
```C++
using MultipatchField< T, Patches >::InternalIdxRangeOnPatch =  typename T<Patch>::discrete_domain_type;
```




<hr>



### typedef base\_type 

_The_ [_**MultipatchType**_](classMultipatchType.md) _from which this class inherits._
```C++
using MultipatchField< T, Patches >::base_type =  MultipatchType<T, Patches...>;
```




<hr>



### typedef discrete\_domain\_type 

_The type of the index ranges that can be used to access this field._ 
```C++
using MultipatchField< T, Patches >::discrete_domain_type =  MultipatchType<InternalIdxRangeOnPatch, Patches...>;
```




<hr>



### typedef element\_type 

_The type of the elements inside the field._ 
```C++
using MultipatchField< T, Patches >::element_type =  typename base_type::example_element::element_type;
```




<hr>



### typedef idx\_type 

_The type used to index the field on the specified patch._ 
```C++
using MultipatchField< T, Patches >::idx_type =  typename InternalIdxRangeOnPatch<Patch>::discrete_element_type;
```




<hr>



### typedef memory\_space 

_The memory space (CPU/GPU) where the data is saved._ 
```C++
using MultipatchField< T, Patches >::memory_space =  typename base_type::example_element::memory_space;
```




<hr>



### typedef span\_type 

_The type of a modifiable reference to this multipatch field._ 
```C++
using MultipatchField< T, Patches >::span_type =  MultipatchField<InternalFieldOnPatch, Patches...>;
```




<hr>



### typedef view\_type 

_The type of a constant reference to this multipatch field._ 
```C++
using MultipatchField< T, Patches >::view_type =  MultipatchField<InternalConstFieldOnPatch, Patches...>;
```




<hr>
## Public Functions Documentation




### function MultipatchField [2/5]

```C++
inline explicit KOKKOS_FUNCTION MultipatchField::MultipatchField (
    T < Patches >... args
) 
```



Instantiate the [**MultipatchField**](classMultipatchField.md) class from an arbitrary number of objects.




**Parameters:**


* `args` The objects to be stored in the class. 




        

<hr>



### function MultipatchField [3/5]

```C++
template<class MultipatchObj, std::enable_if_t<!is_mem_type_v< MultipatchObj >, bool >>
inline KOKKOS_FUNCTION MultipatchField::MultipatchField (
    MultipatchObj & other
) 
```



Create a [**MultipatchField**](classMultipatchField.md) class by copying an instance of another compatible [**MultipatchField**](classMultipatchField.md).


A compatible [**MultipatchField**](classMultipatchField.md) is one which uses all the patches used by this class. The object being copied may include more patches than this [**MultipatchField**](classMultipatchField.md). Further the original [**MultipatchField**](classMultipatchField.md) must store objects of the correct type (the type template may be different but return the same type depending on how it is designed).


This function is not explicit as it is helpful to be able to change between equivalent multipatch definitions if the internal type is the same but the definition comes from different locations in the code.




**Parameters:**


* `other` The equivalent [**MultipatchField**](classMultipatchField.md) being copied. 




        

<hr>



### function MultipatchField [4/5]

```C++
template<class MultipatchObj, std::enable_if_t< is_mem_type_v< MultipatchObj >, bool >>
inline explicit MultipatchField::MultipatchField (
    MultipatchObj & other
) 
```



Create a [**MultipatchField**](classMultipatchField.md) class from a compatible [**MultipatchFieldMem**](classMultipatchFieldMem.md).


A compatible [**MultipatchField**](classMultipatchField.md) is one which uses all the patches used by this class. The object being copied may include more patches than this [**MultipatchField**](classMultipatchField.md). Further the original [**MultipatchField**](classMultipatchField.md) must store objects of the correct type.




**Parameters:**


* `other` The [**MultipatchFieldMem**](classMultipatchFieldMem.md) being accessed. 




        

<hr>



### function MultipatchField [5/5]

```C++
template<template< typename P > typename OtherType, class... OPatches>
inline MultipatchField::MultipatchField (
    MultipatchField < OtherType, OPatches... > && other
) 
```



Create a [**MultipatchField**](classMultipatchField.md) class from an r-value (temporary) instance of another [**MultipatchField**](classMultipatchField.md) which uses the same type for the internal tuple.




**Parameters:**


* `other` The equivalent [**MultipatchField**](classMultipatchField.md) being copied. 




        

<hr>



### function get 

```C++
template<class Patch>
inline KOKKOS_FUNCTION auto MultipatchField::get () const
```



Retrieve an object from the patch that it is defined on.




**Template parameters:**


* [**Patch**](structPatch.md) The patch of the object to be returned. 



**Returns:**

The object on the given patch. 





        

<hr>



### function get\_const\_field 

_Get a_ [_**MultipatchField**_](classMultipatchField.md) _containing constant fields so the values cannot be modified._
```C++
inline KOKKOS_FUNCTION auto MultipatchField::get_const_field () const
```





**Returns:**

A set of constant fields providing access to the fields stored in this class. 





        

<hr>



### function get\_field 

_Get a_ [_**MultipatchField**_](classMultipatchField.md) _containing modifiable fields._
```C++
inline KOKKOS_FUNCTION auto MultipatchField::get_field () 
```





**Returns:**

A set of modifiable fields providing access to the fields stored in this class. 





        

<hr>



### function idx\_range 

_Get a_ [_**MultipatchType**_](classMultipatchType.md) _containing the index ranges on which the fields are defined._
```C++
inline auto MultipatchField::idx_range () const
```





**Returns:**

The set of index ranges on which the set of fields stored in this class are defined. 





        

<hr>



### function span\_cview 

_Get a_ [_**MultipatchField**_](classMultipatchField.md) _containing constant fields so the values cannot be modified. This function matches the DDC name to allow the global get\_const\_field to be defined._
```C++
inline KOKKOS_FUNCTION auto MultipatchField::span_cview () 
```





**Returns:**

A set of constant fields providing access to the fields stored in this class. 





        

<hr>



### function span\_view 

_Get a_ [_**MultipatchField**_](classMultipatchField.md) _containing modifiable fields. This function matches the DDC name to allow the global get\_const\_field to be defined._
```C++
inline KOKKOS_FUNCTION auto MultipatchField::span_view () 
```





**Returns:**

A set of modifiable fields providing access to the fields stored in this class. 





        

<hr>



### function ~MultipatchField 

```C++
KOKKOS_DEFAULTED_FUNCTION MultipatchField::~MultipatchField () noexcept
```




<hr>## Friends Documentation





### friend MultipatchField [1/5]

```C++
template<template< typename P > typename OtherType, class... OPatches>
class MultipatchField::MultipatchField (
    MultipatchField
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/data_types/multipatch_field.hpp`

