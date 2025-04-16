

# File math\_tools.hpp



[**FileList**](files.md) **>** [**math\_tools**](dir_3ced5d1c6eac490d7704c2e023d148d8.md) **>** [**math\_tools.hpp**](math__tools_8hpp.md)

[Go to the source code of this file](math__tools_8hpp_source.md)



* `#include <algorithm>`
* `#include <cmath>`
* `#include <ddc/ddc.hpp>`
* `#include <Kokkos_Core.hpp>`
* `#include "indexed_tensor.hpp"`
* `#include "tensor.hpp"`
* `#include "vector_field.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION double | [**determinant**](#function-determinant) ([**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; RowDim1, RowDim2 &gt;, VectorIndexSet&lt; ColDim1, ColDim2 &gt; &gt; arr) <br> |
|  KOKKOS\_INLINE\_FUNCTION double | [**determinant**](#function-determinant) ([**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; RowDim1, RowDim2, RowDim3 &gt;, VectorIndexSet&lt; ColDim1, ColDim2, ColDim3 &gt; &gt; arr) <br> |
|  std::size\_t | [**factorial**](#function-factorial) (std::size\_t f) <br> |
|  KOKKOS\_INLINE\_FUNCTION [**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; typename ColDim1::Dual, typename ColDim2::Dual &gt;, VectorIndexSet&lt; typename RowDim1::Dual, typename RowDim2::Dual &gt; &gt; | [**inverse**](#function-inverse) ([**DTensor**](classTensor.md)&lt; VectorIndexSet&lt; RowDim1, RowDim2 &gt;, VectorIndexSet&lt; ColDim1, ColDim2 &gt; &gt; arr) <br> |
|  KOKKOS\_INLINE\_FUNCTION constexpr double | [**ipow**](#function-ipow) (double a, std::size\_t i) <br> |
|  KOKKOS\_INLINE\_FUNCTION double | [**ipow**](#function-ipow) (double a, int i) <br> |
|  [**T**](structT.md) | [**max**](#function-max) ([**T**](structT.md) x, [**T**](structT.md) y) <br> |
|  [**T**](structT.md) | [**min**](#function-min) ([**T**](structT.md) x, [**T**](structT.md) y) <br> |
|  [**T**](structT.md) | [**modulo**](#function-modulo) ([**T**](structT.md) x, [**T**](structT.md) y) <br> |
|  KOKKOS\_INLINE\_FUNCTION ElementType | [**norm**](#function-norm) ([**Tensor**](classTensor.md)&lt; ElementType, VectorIndexSetType, VectorIndexSetType &gt; const & metric, [**Tensor**](classTensor.md)&lt; ElementType, vector\_index\_set\_dual\_t&lt; VectorIndexSetType &gt; &gt; const & vec) <br> |
|  void | [**norm**](#function-norm) (ExecSpace exec\_space, DField&lt; IdxRangeType, typename ExecSpace::memory\_space &gt; norm\_vals, [**MetricTensorEvaluator**](classMetricTensorEvaluator.md) const & get\_metric, [**DVectorConstField**](classVectorField.md)&lt; IdxRangeType, VectorIndexSetType, typename ExecSpace::memory\_space &gt; vals) <br> |
|  KOKKOS\_INLINE\_FUNCTION [**T**](structT.md) | [**sum**](#function-sum) (const [**T**](structT.md) \* array, int size) <br> |
|  KOKKOS\_INLINE\_FUNCTION ElementType | [**sum**](#function-sum) (Kokkos::mdspan&lt; ElementType, Kokkos::extents&lt; std::size\_t, Ext &gt;, LayoutPolicy, AccessorPolicy &gt; const & array) <br> |
|  KOKKOS\_INLINE\_FUNCTION ElementType | [**sum**](#function-sum) (Kokkos::mdspan&lt; ElementType, Kokkos::extents&lt; std::size\_t, Ext &gt;, LayoutPolicy, AccessorPolicy &gt; const & array, int start, int end) <br> |




























## Public Functions Documentation




### function determinant 

```C++
template<class RowDim1, class RowDim2, class ColDim1, class ColDim2>
KOKKOS_INLINE_FUNCTION double determinant (
    DTensor < VectorIndexSet< RowDim1, RowDim2 >, VectorIndexSet< ColDim1, ColDim2 > > arr
) 
```




<hr>



### function determinant 

```C++
template<class RowDim1, class RowDim2, class RowDim3, class ColDim1, class ColDim2, class ColDim3>
KOKKOS_INLINE_FUNCTION double determinant (
    DTensor < VectorIndexSet< RowDim1, RowDim2, RowDim3 >, VectorIndexSet< ColDim1, ColDim2, ColDim3 > > arr
) 
```




<hr>



### function factorial 

```C++
inline std::size_t factorial (
    std::size_t f
) 
```




<hr>



### function inverse 

```C++
template<class RowDim1, class RowDim2, class ColDim1, class ColDim2>
KOKKOS_INLINE_FUNCTION DTensor < VectorIndexSet< typename ColDim1::Dual, typename ColDim2::Dual >, VectorIndexSet< typename RowDim1::Dual, typename RowDim2::Dual > > inverse (
    DTensor < VectorIndexSet< RowDim1, RowDim2 >, VectorIndexSet< ColDim1, ColDim2 > > arr
) 
```




<hr>



### function ipow 

```C++
KOKKOS_INLINE_FUNCTION constexpr double ipow (
    double a,
    std::size_t i
) 
```




<hr>



### function ipow 

```C++
KOKKOS_INLINE_FUNCTION double ipow (
    double a,
    int i
) 
```




<hr>



### function max 

```C++
template<typename T>
inline T max (
    T x,
    T y
) 
```




<hr>



### function min 

```C++
template<typename T>
inline T min (
    T x,
    T y
) 
```




<hr>



### function modulo 

```C++
template<typename T>
inline T modulo (
    T x,
    T y
) 
```




<hr>



### function norm 

```C++
template<class ElementType, class VectorIndexSetType>
KOKKOS_INLINE_FUNCTION ElementType norm (
    Tensor < ElementType, VectorIndexSetType, VectorIndexSetType > const & metric,
    Tensor < ElementType, vector_index_set_dual_t< VectorIndexSetType > > const & vec
) 
```




<hr>



### function norm 

```C++
template<class ExecSpace, class IdxRangeType, class MetricTensorEvaluator, class VectorIndexSetType>
void norm (
    ExecSpace exec_space,
    DField< IdxRangeType, typename ExecSpace::memory_space > norm_vals,
    MetricTensorEvaluator const & get_metric,
    DVectorConstField < IdxRangeType, VectorIndexSetType, typename ExecSpace::memory_space > vals
) 
```




<hr>



### function sum 

```C++
template<typename T>
KOKKOS_INLINE_FUNCTION T sum (
    const T * array,
    int size
) 
```




<hr>



### function sum 

```C++
template<class ElementType, class LayoutPolicy, class AccessorPolicy, std::size_t Ext>
KOKKOS_INLINE_FUNCTION ElementType sum (
    Kokkos::mdspan< ElementType, Kokkos::extents< std::size_t, Ext >, LayoutPolicy, AccessorPolicy > const & array
) 
```




<hr>



### function sum 

```C++
template<class ElementType, class LayoutPolicy, class AccessorPolicy, std::size_t Ext>
KOKKOS_INLINE_FUNCTION ElementType sum (
    Kokkos::mdspan< ElementType, Kokkos::extents< std::size_t, Ext >, LayoutPolicy, AccessorPolicy > const & array,
    int start,
    int end
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/math_tools.hpp`

