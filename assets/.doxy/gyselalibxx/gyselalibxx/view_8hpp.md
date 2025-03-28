

# File view.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**view.hpp**](view_8hpp.md)

[Go to the source code of this file](view_8hpp_source.md)



* `#include <array>`
* `#include <ostream>`
* `#include <utility>`
* `#include <Kokkos_Core.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Span1D&lt; double const  &gt; | [**CDSpan1D**](#typedef-cdspan1d)  <br> |
| typedef Span2D&lt; double const  &gt; | [**CDSpan2D**](#typedef-cdspan2d)  <br> |
| typedef Kokkos::View&lt; double[N]&gt; | [**DKokkosView**](#typedef-dkokkosview)  <br> |
| typedef Kokkos::View&lt; double[N], Kokkos::HostSpace &gt; | [**DKokkosView\_h**](#typedef-dkokkosview_h)  <br> |
| typedef Span1D&lt; double &gt; | [**DSpan1D**](#typedef-dspan1d)  <br> |
| typedef Span2D&lt; double &gt; | [**DSpan2D**](#typedef-dspan2d)  <br> |
| typedef View1D&lt; double &gt; | [**DView1D**](#typedef-dview1d)  <br> |
| typedef View2D&lt; double &gt; | [**DView2D**](#typedef-dview2d)  <br> |
| typedef std::array&lt; std::array&lt; double, 2 &gt;, 2 &gt; | [**Matrix\_2x2**](#typedef-matrix_2x2)  <br>_The type of the Jacobian matrix and its inverse._  |
| typedef SpanND&lt; 1, ElementType &gt; | [**Span1D**](#typedef-span1d)  <br> |
| typedef SpanND&lt; 2, ElementType &gt; | [**Span2D**](#typedef-span2d)  <br> |
| typedef Kokkos::mdspan&lt; ElementType, Kokkos::dextents&lt; std::size\_t, N &gt; &gt; | [**SpanND**](#typedef-spannd)  <br> |
| typedef ViewND&lt; 1, ElementType &gt; | [**View1D**](#typedef-view1d)  <br> |
| typedef ViewND&lt; 2, ElementType &gt; | [**View2D**](#typedef-view2d)  <br> |
| typedef SpanND&lt; N, ElementType const  &gt; | [**ViewND**](#typedef-viewnd)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION Span1D&lt; ElementType &gt; | [**as\_span**](#function-as_span) (std::array&lt; ElementType, N &gt; & arr) noexcept<br> |
|  KOKKOS\_FUNCTION Span1D&lt; const ElementType &gt; | [**as\_span**](#function-as_span) (std::array&lt; ElementType, N &gt; const & arr) noexcept<br> |
|  std::ostream & | [**operator&lt;&lt;**](#function-operator) (std::ostream & os, Kokkos::mdspan&lt; ElementType, Extents, Layout, Accessor &gt; const & s) <br> |




























## Public Types Documentation




### typedef CDSpan1D 

```C++
using CDSpan1D =  Span1D<double const>;
```




<hr>



### typedef CDSpan2D 

```C++
using CDSpan2D =  Span2D<double const>;
```




<hr>



### typedef DKokkosView 

```C++
using DKokkosView =  Kokkos::View<double[N]>;
```




<hr>



### typedef DKokkosView\_h 

```C++
using DKokkosView_h =  Kokkos::View<double[N], Kokkos::HostSpace>;
```




<hr>



### typedef DSpan1D 

```C++
using DSpan1D =  Span1D<double>;
```




<hr>



### typedef DSpan2D 

```C++
using DSpan2D =  Span2D<double>;
```




<hr>



### typedef DView1D 

```C++
using DView1D =  View1D<double>;
```




<hr>



### typedef DView2D 

```C++
using DView2D =  View2D<double>;
```




<hr>



### typedef Matrix\_2x2 

_The type of the Jacobian matrix and its inverse._ 
```C++
using Matrix_2x2 =  std::array<std::array<double, 2>, 2>;
```




<hr>



### typedef Span1D 

```C++
using Span1D =  SpanND<1, ElementType>;
```




<hr>



### typedef Span2D 

```C++
using Span2D =  SpanND<2, ElementType>;
```




<hr>



### typedef SpanND 

```C++
using SpanND =  Kokkos::mdspan<ElementType, Kokkos::dextents<std::size_t, N> >;
```




<hr>



### typedef View1D 

```C++
using View1D =  ViewND<1, ElementType>;
```




<hr>



### typedef View2D 

```C++
using View2D =  ViewND<2, ElementType>;
```




<hr>



### typedef ViewND 

```C++
using ViewND =  SpanND<N, ElementType const>;
```




<hr>
## Public Functions Documentation




### function as\_span 

```C++
template<class ElementType, std::size_t N>
KOKKOS_FUNCTION Span1D< ElementType > as_span (
    std::array< ElementType, N > & arr
) noexcept
```




<hr>



### function as\_span 

```C++
template<class ElementType, std::size_t N>
KOKKOS_FUNCTION Span1D< const ElementType > as_span (
    std::array< ElementType, N > const & arr
) noexcept
```




<hr>



### function operator&lt;&lt; 

```C++
template<class ElementType, class Extents, class Layout, class Accessor>
std::ostream & operator<< (
    std::ostream & os,
    Kokkos::mdspan< ElementType, Extents, Layout, Accessor > const & s
) 
```



Convenient function to dump a mdspan, it recursively prints all dimensions. Disclaimer: use with caution for large arrays 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/view.hpp`

