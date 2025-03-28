

# Struct IdxRangeSliceIterator

**template &lt;class Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**IdxRangeSliceIterator**](structIdxRangeSliceIterator.md)



_An iterator type for the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _._

* `#include <idx_range_slice.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::ptrdiff\_t | [**difference\_type**](#typedef-difference_type)  <br>_The type that can be used to increment the iterator._  |
| typedef std::random\_access\_iterator\_tag | [**iterator\_category**](#typedef-iterator_category)  <br>_The type of iterator._  |
| typedef IdxStep&lt; Grid1D &gt; | [**stride\_type**](#typedef-stride_type)  <br>_The type of the stride between values._  |
| typedef Idx&lt; Grid1D &gt; | [**value\_type**](#typedef-value_type)  <br>_The type of the values stored in the iterator._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**IdxRangeSliceIterator**](#function-idxrangesliceiterator-12) () = default<br> |
|  KOKKOS\_FUNCTION constexpr | [**IdxRangeSliceIterator**](#function-idxrangesliceiterator-22) (Idx&lt; Grid1D &gt; value, IdxStep&lt; Grid1D &gt; stride) <br>_Build an iterator using the current value and the distance to the following element._  |
|  KOKKOS\_FUNCTION constexpr Idx&lt; Grid1D &gt; | [**operator\***](#function-operator) () noexcept const<br>_Get the value referred to by the iterator._  |
|  KOKKOS\_FUNCTION constexpr [**IdxRangeSliceIterator**](structIdxRangeSliceIterator.md) & | [**operator++**](#function-operator_1) () <br>_The prefix increment operator._  |
|  KOKKOS\_FUNCTION constexpr [**IdxRangeSliceIterator**](structIdxRangeSliceIterator.md) | [**operator++**](#function-operator_2) (int) <br>_The postfix increment operator._  |
|  KOKKOS\_FUNCTION constexpr [**IdxRangeSliceIterator**](structIdxRangeSliceIterator.md) & | [**operator+=**](#function-operator_3) ([**difference\_type**](structIdxRangeSliceIterator.md#typedef-difference_type) n) <br>_Increment the current iterator by n elements._  |
|  KOKKOS\_FUNCTION constexpr [**IdxRangeSliceIterator**](structIdxRangeSliceIterator.md) & | [**operator--**](#function-operator-) () <br>_The prefix decrement operator._  |
|  KOKKOS\_FUNCTION constexpr [**IdxRangeSliceIterator**](structIdxRangeSliceIterator.md) | [**operator--**](#function-operator-_1) (int) <br>_The postfix decrement operator._  |
|  KOKKOS\_FUNCTION constexpr [**IdxRangeSliceIterator**](structIdxRangeSliceIterator.md) & | [**operator-=**](#function-operator-) ([**difference\_type**](structIdxRangeSliceIterator.md#typedef-difference_type) n) <br>_Decrement the current iterator by n elements._  |
|  KOKKOS\_FUNCTION constexpr Idx&lt; Grid1D &gt; | [**operator[]**](#function-operator_4) ([**difference\_type**](structIdxRangeSliceIterator.md#typedef-difference_type) n) const<br>_Access the n-th following element._  |




























## Public Types Documentation




### typedef difference\_type 

_The type that can be used to increment the iterator._ 
```C++
using IdxRangeSliceIterator< Grid1D >::difference_type =  std::ptrdiff_t;
```




<hr>



### typedef iterator\_category 

_The type of iterator._ 
```C++
using IdxRangeSliceIterator< Grid1D >::iterator_category =  std::random_access_iterator_tag;
```




<hr>



### typedef stride\_type 

_The type of the stride between values._ 
```C++
using IdxRangeSliceIterator< Grid1D >::stride_type =  IdxStep<Grid1D>;
```




<hr>



### typedef value\_type 

_The type of the values stored in the iterator._ 
```C++
using IdxRangeSliceIterator< Grid1D >::value_type =  Idx<Grid1D>;
```




<hr>
## Public Functions Documentation




### function IdxRangeSliceIterator [1/2]

```C++
KOKKOS_DEFAULTED_FUNCTION IdxRangeSliceIterator::IdxRangeSliceIterator () = default
```




<hr>



### function IdxRangeSliceIterator [2/2]

_Build an iterator using the current value and the distance to the following element._ 
```C++
inline explicit KOKKOS_FUNCTION constexpr IdxRangeSliceIterator::IdxRangeSliceIterator (
    Idx< Grid1D > value,
    IdxStep< Grid1D > stride
) 
```





**Parameters:**


* `value` The value of the discrete sub-index range element. 
* `stride` The stride between consecutive sub-index range elements. 




        

<hr>



### function operator\* 

_Get the value referred to by the iterator._ 
```C++
inline KOKKOS_FUNCTION constexpr Idx< Grid1D > IdxRangeSliceIterator::operator* () noexcept const
```





**Returns:**

The value referred to by the iterator. 





        

<hr>



### function operator++ 

_The prefix increment operator._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator & IdxRangeSliceIterator::operator++ () 
```





**Returns:**

A reference to the current incremented iterator. 





        

<hr>



### function operator++ 

_The postfix increment operator._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator IdxRangeSliceIterator::operator++ (
    int
) 
```





**Returns:**

A reference to the non-incremented iterator. 





        

<hr>



### function operator+= 

_Increment the current iterator by n elements._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator & IdxRangeSliceIterator::operator+= (
    difference_type n
) 
```





**Parameters:**


* `n` The number of strides by which the iterator should be incremented. 



**Returns:**

A reference to the current incremented iterator. 





        

<hr>



### function operator-- 

_The prefix decrement operator._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator & IdxRangeSliceIterator::operator-- () 
```





**Returns:**

A reference to the current decremented iterator. 





        

<hr>



### function operator-- 

_The postfix decrement operator._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator IdxRangeSliceIterator::operator-- (
    int
) 
```





**Returns:**

A reference to the non-decremented iterator. 





        

<hr>



### function operator-= 

_Decrement the current iterator by n elements._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator & IdxRangeSliceIterator::operator-= (
    difference_type n
) 
```





**Parameters:**


* `n` The number of strides by which the iterator should be decremented. 



**Returns:**

A reference to the current decremented iterator. 





        

<hr>



### function operator[] 

_Access the n-th following element._ 
```C++
inline KOKKOS_FUNCTION constexpr Idx< Grid1D > IdxRangeSliceIterator::operator[] (
    difference_type n
) const
```





**Parameters:**


* `n` The index of the element. 



**Returns:**

The Idx n-th position in the sub-index range following the value indicated by this iterator. 





        

<hr>## Friends Documentation





### friend operator!= 

_Compare iterator non-equality._ 
```C++
inline KOKKOS_FUNCTION constexpr bool IdxRangeSliceIterator::operator!= (
    IdxRangeSliceIterator const & xx,
    IdxRangeSliceIterator const & yy
) 
```




<hr>



### friend operator+ 

_Increment an iterator by n elements._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator IdxRangeSliceIterator::operator+ (
    IdxRangeSliceIterator i,
    difference_type n
) 
```




<hr>



### friend operator+ 

_Increment an iterator by n elements._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator IdxRangeSliceIterator::operator+ (
    difference_type n,
    IdxRangeSliceIterator i
) 
```




<hr>



### friend operator- 

_Decrement an iterator by n elements._ 
```C++
inline KOKKOS_FUNCTION constexpr IdxRangeSliceIterator IdxRangeSliceIterator::operator- (
    IdxRangeSliceIterator i,
    difference_type n
) 
```




<hr>



### friend operator- 

_Decrement an iterator by n elements._ 
```C++
inline KOKKOS_FUNCTION constexpr difference_type IdxRangeSliceIterator::operator- (
    IdxRangeSliceIterator const & xx,
    IdxRangeSliceIterator const & yy
) 
```




<hr>



### friend operator&lt; 

_Compare the order of iterators._ 
```C++
inline KOKKOS_FUNCTION constexpr bool IdxRangeSliceIterator::operator< (
    IdxRangeSliceIterator const & xx,
    IdxRangeSliceIterator const & yy
) 
```




<hr>



### friend operator&lt;= 

_Compare the order of iterators._ 
```C++
inline KOKKOS_FUNCTION constexpr bool IdxRangeSliceIterator::operator<= (
    IdxRangeSliceIterator const & xx,
    IdxRangeSliceIterator const & yy
) 
```




<hr>



### friend operator== 

_Compare iterator equality._ 
```C++
inline KOKKOS_FUNCTION constexpr bool IdxRangeSliceIterator::operator== (
    IdxRangeSliceIterator const & xx,
    IdxRangeSliceIterator const & yy
) 
```




<hr>



### friend operator&gt; 

_Compare the order of iterators._ 
```C++
inline KOKKOS_FUNCTION constexpr bool IdxRangeSliceIterator::operator> (
    IdxRangeSliceIterator const & xx,
    IdxRangeSliceIterator const & yy
) 
```




<hr>



### friend operator&gt;= 

_Compare the order of iterators._ 
```C++
inline KOKKOS_FUNCTION constexpr bool IdxRangeSliceIterator::operator>= (
    IdxRangeSliceIterator const & xx,
    IdxRangeSliceIterator const & yy
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/idx_range_slice.hpp`

