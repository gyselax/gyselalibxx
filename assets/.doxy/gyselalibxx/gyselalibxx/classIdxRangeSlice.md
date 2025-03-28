

# Class IdxRangeSlice

**template &lt;class... Dims&gt;**



[**ClassList**](annotated.md) **>** [**IdxRangeSlice**](classIdxRangeSlice.md)



_A class which describes a collection of equally spaced Idxs which form a index range._ [More...](#detailed-description)

* `#include <idx_range_slice.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**IdxRangeSlice**](#function-idxrangeslice-24) () = default<br>_Default constructor for_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _creating an empty index range._ |
|  KOKKOS\_FUNCTION | [**IdxRangeSlice**](#function-idxrangeslice-34) (Idx&lt; Dims... &gt; front, IdxStep&lt; Dims... &gt; size, IdxStep&lt; Dims... &gt; stride) <br>_Build a_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _from vectors of valid Idxs in each dimension._ |
|  KOKKOS\_FUNCTION | [**IdxRangeSlice**](#function-idxrangeslice-44) (DDoms const &... valid\_indices) <br>_Build a_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _from a set of 1D IdxRangeSlices._ |
|  KOKKOS\_FUNCTION constexpr Idx&lt; Dims... &gt; | [**back**](#function-back) () noexcept const<br>_Get the last element in the index range slice._  |
|  KOKKOS\_FUNCTION auto | [**begin**](#function-begin) () const<br>_Get the iterator to the first element of the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _. The elements are pairs of Idxs and indices._ |
|  KOKKOS\_FUNCTION bool | [**contains**](#function-contains-12) (Idx&lt; DDims... &gt; elem) const<br>_Check if the specified Idx is found in this index range slice._  |
|  KOKKOS\_FUNCTION bool | [**contains**](#function-contains-22) (IdxRange&lt; DDims... &gt; idx\_range) const<br>_Check if all elements of the specified IdxRange are found in this index range slice._  |
|  KOKKOS\_FUNCTION auto | [**end**](#function-end) () const<br>_Get the iterator to the end of the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _. The elements are pairs of Idxs and indices._ |
|  KOKKOS\_FUNCTION constexpr IdxStep&lt; QueryDDim &gt; | [**extent**](#function-extent) () noexcept const<br>_Get the size of the index range slice in the specified dimension._  |
|  KOKKOS\_FUNCTION constexpr IdxStep&lt; Dims... &gt; | [**extents**](#function-extents) () noexcept const<br>_Get the size of the index range slice in each dimension._  |
|  KOKKOS\_FUNCTION constexpr Idx&lt; Dims... &gt; | [**front**](#function-front) () noexcept const<br>_Get the first element in the index range slice._  |
|  KOKKOS\_FUNCTION std::size\_t | [**get\_index**](#function-get_index) (Idx&lt; Dim &gt; elem) const<br>_Get the index of the Idx within the index range slice. This function is particularly useful to index an mdspan over the index range slice._  |
|  KOKKOS\_FUNCTION constexpr std::size\_t | [**size**](#function-size) () const<br>_Get the total number of elements in the index range slice._  |
|  KOKKOS\_FUNCTION constexpr IdxStep&lt; QueryDim &gt; | [**stride**](#function-stride) () noexcept const<br>_Get the stride from one element of the index range slice to another._  |
|  KOKKOS\_FUNCTION constexpr IdxStep&lt; Dims... &gt; | [**strides**](#function-strides) () noexcept const<br>_Get the strides from one element of the index range slice to another._  |




























## Detailed Description


This class should eventually be replaced by a DDC functionality when this becomes available. 


    
## Public Functions Documentation




### function IdxRangeSlice [2/4]

_Default constructor for_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _creating an empty index range._
```C++
IdxRangeSlice::IdxRangeSlice () = default
```




<hr>



### function IdxRangeSlice [3/4]

_Build a_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _from vectors of valid Idxs in each dimension._
```C++
inline KOKKOS_FUNCTION IdxRangeSlice::IdxRangeSlice (
    Idx< Dims... > front,
    IdxStep< Dims... > size,
    IdxStep< Dims... > stride
) 
```





**Parameters:**


* `front` The Idx describing the first point in the index range. 
* `size` The number of elements along each direction. 
* `stride` The IdxStep distance between subsequent elements of the index range. 




        

<hr>



### function IdxRangeSlice [4/4]

_Build a_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _from a set of 1D IdxRangeSlices._
```C++
template<class... DDoms, class>
inline KOKKOS_FUNCTION IdxRangeSlice::IdxRangeSlice (
    DDoms const &... valid_indices
) 
```





**Parameters:**


* `valid_indices` The IdxRangeSlices which comprise this [**IdxRangeSlice**](classIdxRangeSlice.md). 




        

<hr>



### function back 

_Get the last element in the index range slice._ 
```C++
inline KOKKOS_FUNCTION constexpr Idx< Dims... > IdxRangeSlice::back () noexcept const
```





**Returns:**

The last element in the index range slice. 





        

<hr>



### function begin 

_Get the iterator to the first element of the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _. The elements are pairs of Idxs and indices._
```C++
inline KOKKOS_FUNCTION auto IdxRangeSlice::begin () const
```





**Returns:**

The iterator to the first element of the [**IdxRangeSlice**](classIdxRangeSlice.md). 





        

<hr>



### function contains [1/2]

_Check if the specified Idx is found in this index range slice._ 
```C++
template<class... DDims>
inline KOKKOS_FUNCTION bool IdxRangeSlice::contains (
    Idx< DDims... > elem
) const
```





**Parameters:**


* `elem` The element which may or may not be in this index range slice.



**Returns:**

bool True if the element is found in this index range slice, False otherwise. 





        

<hr>



### function contains [2/2]

_Check if all elements of the specified IdxRange are found in this index range slice._ 
```C++
template<class... DDims>
inline KOKKOS_FUNCTION bool IdxRangeSlice::contains (
    IdxRange< DDims... > idx_range
) const
```





**Parameters:**


* `idx_range` The index range which may or may not be in this index range slice.



**Returns:**

bool True if the all elements of the specified IdxRange are found in this index range slice. False otherwise. 





        

<hr>



### function end 

_Get the iterator to the end of the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _. The elements are pairs of Idxs and indices._
```C++
inline KOKKOS_FUNCTION auto IdxRangeSlice::end () const
```





**Returns:**

The iterator to the end of the [**IdxRangeSlice**](classIdxRangeSlice.md). 





        

<hr>



### function extent 

_Get the size of the index range slice in the specified dimension._ 
```C++
template<class QueryDDim>
inline KOKKOS_FUNCTION constexpr IdxStep< QueryDDim > IdxRangeSlice::extent () noexcept const
```





**Template parameters:**


* `QueryDDim` The dimension in which the size is requested.



**Returns:**

A IdxStep describing the size of the index range slice in the specified dimension. 





        

<hr>



### function extents 

_Get the size of the index range slice in each dimension._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxStep< Dims... > IdxRangeSlice::extents () noexcept const
```





**Returns:**

A IdxStep describing the size of the index range slice in each dimension. 





        

<hr>



### function front 

_Get the first element in the index range slice._ 
```C++
inline KOKKOS_FUNCTION constexpr Idx< Dims... > IdxRangeSlice::front () noexcept const
```





**Returns:**

The first element in the index range slice. 





        

<hr>



### function get\_index 

_Get the index of the Idx within the index range slice. This function is particularly useful to index an mdspan over the index range slice._ 
```C++
template<class Dim, class>
inline KOKKOS_FUNCTION std::size_t IdxRangeSlice::get_index (
    Idx< Dim > elem
) const
```





**Parameters:**


* `elem` A 1D Idx which is inside the index range slice.



**Returns:**

The index of the element. 





        

<hr>



### function size 

_Get the total number of elements in the index range slice._ 
```C++
inline KOKKOS_FUNCTION constexpr std::size_t IdxRangeSlice::size () const
```





**Returns:**

The total number of elements in the index range slice. 





        

<hr>



### function stride 

_Get the stride from one element of the index range slice to another._ 
```C++
template<class QueryDim>
inline KOKKOS_FUNCTION constexpr IdxStep< QueryDim > IdxRangeSlice::stride () noexcept const
```





**Template parameters:**


* `QueryDim` The dimension being queried



**Returns:**

A 1D IdxStep describing a stride from one element to another. 





        

<hr>



### function strides 

_Get the strides from one element of the index range slice to another._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxStep< Dims... > IdxRangeSlice::strides () noexcept const
```





**Returns:**

A IdxStep describing the strides from one element to another. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/idx_range_slice.hpp`

