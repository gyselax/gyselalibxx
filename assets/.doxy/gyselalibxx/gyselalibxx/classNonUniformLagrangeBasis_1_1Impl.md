

# Class NonUniformLagrangeBasis::Impl

**template &lt;class DDim, class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**NonUniformLagrangeBasis**](classNonUniformLagrangeBasis.md) **>** [**Impl**](classNonUniformLagrangeBasis_1_1Impl.md)



_Storage class of the static attributes of the discrete dimension._ [More...](#detailed-description)

* `#include <lagrange_basis_non_uniform.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**NonUniformLagrangeKnots**](structNonUniformLagrangeKnots.md)&lt; DDim &gt; | [**knot\_grid**](#typedef-knot_grid)  <br>_The type of the knots defining the B-splines._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Impl**](#function-impl-24) () = default<br> |
|   | [**Impl**](#function-impl-34) (IdxRange&lt; Grid1D &gt; break\_point\_domain) <br>_Initialise the possible Lagrange bases._  |
|   | [**Impl**](#function-impl-44) ([**Impl**](classNonUniformLagrangeBasis_1_1Impl.md)&lt; DDim, OriginMemorySpace &gt; const & impl) <br>_Copy-constructs from another_ [_**Impl**_](classNonUniformLagrangeBasis_1_1Impl.md) _with a different Kokkos memory space._ |
|  KOKKOS\_INLINE\_FUNCTION IdxRange&lt; [**knot\_grid**](classNonUniformLagrangeBasis_1_1Impl.md#typedef-knot_grid) &gt; | [**break\_point\_domain**](#function-break_point_domain) () const<br>_Returns the index range of the break points._  |
|  KOKKOS\_INLINE\_FUNCTION void | [**eval\_basis**](#function-eval_basis) (Span1D&lt; DataType &gt; values, [**coord\_type**](classNonUniformLagrangeBasis.md#typedef-coord_type) const & x, Idx&lt; [**knot\_grid**](classNonUniformLagrangeBasis_1_1Impl.md#typedef-knot_grid) &gt; poly\_start) const<br>_Evaluate the selected set of bases at the coordinate._  |
|  KOKKOS\_INLINE\_FUNCTION IdxRange&lt; [**knot\_grid**](classNonUniformLagrangeBasis_1_1Impl.md#typedef-knot_grid) &gt; | [**full\_domain**](#function-full_domain) () const<br>_Returns the index range including eventual duplicate values for the periodic case._  |
|  KOKKOS\_INLINE\_FUNCTION DataType | [**length**](#function-length) () noexcept const<br>_Returns the length of the domain._  |
|  KOKKOS\_INLINE\_FUNCTION [**coord\_type**](classNonUniformLagrangeBasis.md#typedef-coord_type) | [**rmax**](#function-rmax) () noexcept const<br>_Returns the coordinate of the upper bound of the domain on which the B-splines are defined._  |
|  KOKKOS\_INLINE\_FUNCTION [**coord\_type**](classNonUniformLagrangeBasis.md#typedef-coord_type) | [**rmin**](#function-rmin) () noexcept const<br>_Returns the coordinate of the lower bound of the domain on which the B-splines are defined._  |




























## Detailed Description




**Template parameters:**


* `DDim` The name of the discrete dimension. 
* `MemorySpace` The Kokkos memory space where the attributes are being stored. 




    
## Public Types Documentation




### typedef knot\_grid 

_The type of the knots defining the B-splines._ 
```C++
using NonUniformLagrangeBasis< Dim, D, DataType >::Impl< DDim, MemorySpace >::knot_grid =  NonUniformLagrangeKnots<DDim>;
```




<hr>
## Public Functions Documentation




### function Impl [2/4]

```C++
NonUniformLagrangeBasis::Impl::Impl () = default
```




<hr>



### function Impl [3/4]

_Initialise the possible Lagrange bases._ 
```C++
template<class Grid1D>
inline explicit NonUniformLagrangeBasis::Impl::Impl (
    IdxRange< Grid1D > break_point_domain
) 
```



Initialise the class such that Lagrange bases can be evaluated on domains derived from the break point domain. 


        

<hr>



### function Impl [4/4]

_Copy-constructs from another_ [_**Impl**_](classNonUniformLagrangeBasis_1_1Impl.md) _with a different Kokkos memory space._
```C++
template<class OriginMemorySpace>
inline explicit NonUniformLagrangeBasis::Impl::Impl (
    Impl < DDim, OriginMemorySpace > const & impl
) 
```





**Parameters:**


* `impl` A reference to the other [**Impl**](classNonUniformLagrangeBasis_1_1Impl.md). 




        

<hr>



### function break\_point\_domain 

_Returns the index range of the break points._ 
```C++
inline KOKKOS_INLINE_FUNCTION IdxRange< knot_grid > NonUniformLagrangeBasis::Impl::break_point_domain () const
```





**Returns:**

The index range describing the break points. 





        

<hr>



### function eval\_basis 

_Evaluate the selected set of bases at the coordinate._ 
```C++
inline KOKKOS_INLINE_FUNCTION void NonUniformLagrangeBasis::Impl::eval_basis (
    Span1D< DataType > values,
    coord_type const & x,
    Idx< knot_grid > poly_start
) const
```



Evaluate all d+1 bases which span the domain [coordinate(poly\_start), coordinate(poly\_start+d)] at the coordinate x.




**Parameters:**


* `values` The values of each basis at the coordinate x. 
* `x` The coordinate where the bases are evaluated. 
* `poly_start` The index of the first of the d+1 knots describing the set of bases to be evaluated. 




        

<hr>



### function full\_domain 

_Returns the index range including eventual duplicate values for the periodic case._ 
```C++
inline KOKKOS_INLINE_FUNCTION IdxRange< knot_grid > NonUniformLagrangeBasis::Impl::full_domain () const
```





**Returns:**

The index range including eventual duplicate values. 





        

<hr>



### function length 

_Returns the length of the domain._ 
```C++
inline KOKKOS_INLINE_FUNCTION DataType NonUniformLagrangeBasis::Impl::length () noexcept const
```





**Returns:**

The length of the domain. 





        

<hr>



### function rmax 

_Returns the coordinate of the upper bound of the domain on which the B-splines are defined._ 
```C++
inline KOKKOS_INLINE_FUNCTION coord_type NonUniformLagrangeBasis::Impl::rmax () noexcept const
```





**Returns:**

Coordinate of the upper bound of the domain. 





        

<hr>



### function rmin 

_Returns the coordinate of the lower bound of the domain on which the B-splines are defined._ 
```C++
inline KOKKOS_INLINE_FUNCTION coord_type NonUniformLagrangeBasis::Impl::rmin () noexcept const
```





**Returns:**

Coordinate of the lower bound of the domain. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/lagrange_basis_non_uniform.hpp`

