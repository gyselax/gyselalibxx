

# Class MultipatchType

**template &lt;template&lt; typename P &gt; typename T, class... Patches&gt;**



[**ClassList**](annotated.md) **>** [**MultipatchType**](classMultipatchType.md)



_A class to store several objects that are of a type which is templated by the patch._ [More...](#detailed-description)

* `#include <multipatch_type.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::detail::TypeSeq&lt; Patches... &gt; | [**PatchOrdering**](#typedef-patchordering)  <br>_A tag storing the order of Patches in this_ [_**MultipatchType**_](classMultipatchType.md) _._ |
| typedef [**T**](structT.md)&lt; ddc::type\_seq\_element\_t&lt; 0, [**PatchOrdering**](classMultipatchType.md#typedef-patchordering) &gt; &gt; | [**example\_element**](#typedef-example_element)  <br>_The type of one of the elements of the_ [_**MultipatchType**_](classMultipatchType.md) _. This can be used to check that types are as expected using functions such as ddc::is\_chunk\_v._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**MultipatchType**](#function-multipatchtype-35) ([**T**](structT.md)&lt; Patches &gt;... args) <br> |
|  KOKKOS\_FUNCTION | [**MultipatchType**](#function-multipatchtype-45) ([**MultipatchType**](classMultipatchType.md)&lt; OtherType, OPatches... &gt; const & other) <br> |
|   | [**MultipatchType**](#function-multipatchtype-55) ([**MultipatchType**](classMultipatchType.md)&lt; OtherType, OPatches... &gt; && other) <br> |
|  KOKKOS\_FUNCTION [**T**](structT.md)&lt; [**Patch**](structPatch.md) &gt; | [**get**](#function-get) () const<br> |
|  KOKKOS\_FUNCTION std::tuple&lt; [**T**](structT.md)&lt; Patches &gt;... &gt; const & | [**get\_tuple**](#function-get_tuple) () const<br>_Get a constant reference to the tuple of objects stored inside this_ [_**MultipatchType**_](classMultipatchType.md) _._ |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**~MultipatchType**](#function-multipatchtype) () noexcept<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**size**](#function-size) () <br>_Get the number of objects stored inside the class. This is equal to the number of patches._  |






## Protected Attributes

| Type | Name |
| ---: | :--- |
|  std::tuple&lt; [**T**](structT.md)&lt; Patches &gt;... &gt; | [**m\_tuple**](#variable-m_tuple)  <br>_The internal tuple containing the data._  |
















## Protected Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**MultipatchType**](#function-multipatchtype-25) (std::tuple&lt; [**T**](structT.md)&lt; Patches &gt;... &gt; && tuple) <br> |




## Detailed Description


On a multipatch domain when we have objects and types defined on different patches, e.g. fields. They can be stored in this class and then be accessed by the patch they are defined on.




**Template parameters:**


* [**T**](structT.md) The type of the objects that are stored on the given patches. 
* `Patches` The patches of the objects in the same order of the patches that the given objects are defined on.



**Warning:**

The objects have to be defined on different patches. Otherwise retrieving them by their patch is ill-defined. 





    
## Public Types Documentation




### typedef PatchOrdering 

_A tag storing the order of Patches in this_ [_**MultipatchType**_](classMultipatchType.md) _._
```C++
using MultipatchType< T, Patches >::PatchOrdering =  ddc::detail::TypeSeq<Patches...>;
```




<hr>



### typedef example\_element 

_The type of one of the elements of the_ [_**MultipatchType**_](classMultipatchType.md) _. This can be used to check that types are as expected using functions such as ddc::is\_chunk\_v._
```C++
using MultipatchType< T, Patches >::example_element =  T<ddc::type_seq_element_t<0, PatchOrdering> >;
```




<hr>
## Public Functions Documentation




### function MultipatchType [3/5]

```C++
inline explicit KOKKOS_FUNCTION MultipatchType::MultipatchType (
    T < Patches >... args
) 
```



Instantiate the [**MultipatchType**](classMultipatchType.md) class from an arbitrary number of objects.




**Parameters:**


* `args` The objects to be stored in the class. 




        

<hr>



### function MultipatchType [4/5]

```C++
template<template< typename P > typename OtherType, class... OPatches>
inline KOKKOS_FUNCTION MultipatchType::MultipatchType (
    MultipatchType < OtherType, OPatches... > const & other
) 
```



Create a [**MultipatchType**](classMultipatchType.md) class by copying an instance of another compatible [**MultipatchType**](classMultipatchType.md).


A compatible [**MultipatchType**](classMultipatchType.md) is one which uses all the patches used by this class. The object being copied may include more patches than this [**MultipatchType**](classMultipatchType.md). Further the original [**MultipatchType**](classMultipatchType.md) must store objects of the correct type (the type template may be different but return the same type depending on how it is designed.


This function is not explicit as it is helpful to be able to change between equivalent multipatch definitions if the internal type is the same but the definition comes from different locations in the code.




**Parameters:**


* `other` The equivalent [**MultipatchType**](classMultipatchType.md) being copied. 




        

<hr>



### function MultipatchType [5/5]

```C++
template<template< typename P > typename OtherType, class... OPatches>
inline MultipatchType::MultipatchType (
    MultipatchType < OtherType, OPatches... > && other
) 
```



Create a [**MultipatchType**](classMultipatchType.md) class from an r-value (temporary) instance of another [**MultipatchType**](classMultipatchType.md) which uses the same type for the internal tuple.




**Parameters:**


* `other` The equivalent [**MultipatchType**](classMultipatchType.md) being copied. 




        

<hr>



### function get 

```C++
template<class Patch, std::enable_if_t<!has_data_access_methods_v< T < Patch > >, bool >>
inline KOKKOS_FUNCTION T < Patch > MultipatchType::get () const
```



Retrieve an object from the patch that it is defined on.




**Template parameters:**


* [**Patch**](structPatch.md) The patch of the object to be returned. 



**Returns:**

The object on the given patch. 





        

<hr>



### function get\_tuple 

_Get a constant reference to the tuple of objects stored inside this_ [_**MultipatchType**_](classMultipatchType.md) _._
```C++
inline KOKKOS_FUNCTION std::tuple< T < Patches >... > const & MultipatchType::get_tuple () const
```





**Returns:**

A constant reference to the tuple of objects stored inside this [**MultipatchType**](classMultipatchType.md). 





        

<hr>



### function ~MultipatchType 

```C++
KOKKOS_DEFAULTED_FUNCTION MultipatchType::~MultipatchType () noexcept
```




<hr>
## Public Static Functions Documentation




### function size 

_Get the number of objects stored inside the class. This is equal to the number of patches._ 
```C++
static inline constexpr std::size_t MultipatchType::size () 
```





**Returns:**

Number of elements stored in the tuple of the class. 





        

<hr>
## Protected Attributes Documentation




### variable m\_tuple 

_The internal tuple containing the data._ 
```C++
std::tuple<T<Patches>...> MultipatchType< T, Patches >::m_tuple;
```




<hr>
## Protected Functions Documentation




### function MultipatchType [2/5]

```C++
inline explicit KOKKOS_FUNCTION MultipatchType::MultipatchType (
    std::tuple< T < Patches >... > && tuple
) 
```



A constructor for sub-classes which can build the necessary tuple directly following their own rules.




**Parameters:**


* `tuple` The internal tuple. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/data_types/multipatch_type.hpp`

