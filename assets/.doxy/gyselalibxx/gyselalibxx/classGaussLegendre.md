

# Class GaussLegendre

**template &lt;class GLGrid, std::size\_t NPoints&gt;**



[**ClassList**](annotated.md) **>** [**GaussLegendre**](classGaussLegendre.md)



_An operator for constructing a Gauss-Legendre quadrature._ [More...](#detailed-description)

* `#include <gauss_legendre_integration.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef GLGrid | [**Grid1D**](#typedef-grid1d)  <br>_The grid on which the quadrature scheme is defined._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**GaussLegendre**](#function-gausslegendre-13) (InputIt mesh\_edges\_begin, InputIt mesh\_edges\_end) <br>_A constructor of the_ [_**GaussLegendre**_](classGaussLegendre.md) _class._ |
|   | [**GaussLegendre**](#function-gausslegendre-23) (std::initializer\_list&lt; Coord&lt; Dim &gt; &gt; const mesh\_edges) <br>_A constructor of the_ [_**GaussLegendre**_](classGaussLegendre.md) _class._ |
|   | [**GaussLegendre**](#function-gausslegendre-33) (ConstField&lt; Coord&lt; Dim &gt;, IdxRange&lt; [**Grid1D**](classGaussLegendre.md#typedef-grid1d) &gt;, Kokkos::HostSpace &gt; mesh\_edges) <br>_A constructor of the_ [_**GaussLegendre**_](classGaussLegendre.md) _class._ |
|  DFieldMem&lt; IdxRange&lt; GLGrid &gt;, typename ExecSpace::memory\_space &gt; | [**gauss\_legendre\_coefficients**](#function-gauss_legendre_coefficients) () const<br>_Get a FieldMem containing the coefficients for the Gauss-Legendre quadrature._  |
|  IdxRange&lt; GLGrid &gt; | [**get\_idx\_range**](#function-get_idx_range) () const<br>_Get the index range of the points of the Gauss-Legendre quadrature._  |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**order**](#function-order) () <br>_The order of the quadrature scheme._  |


























## Detailed Description




**Template parameters:**


* `GLGrid` The grid describing the Gauss-Legendre points. 
* `NPoints` The number of points in the Gauss-Legendre scheme 




    
## Public Types Documentation




### typedef Grid1D 

_The grid on which the quadrature scheme is defined._ 
```C++
using GaussLegendre< GLGrid, NPoints >::Grid1D =  GLGrid;
```




<hr>
## Public Functions Documentation




### function GaussLegendre [1/3]

_A constructor of the_ [_**GaussLegendre**_](classGaussLegendre.md) _class._
```C++
template<class InputIt>
inline GaussLegendre::GaussLegendre (
    InputIt mesh_edges_begin,
    InputIt mesh_edges_end
) 
```





**Parameters:**


* `mesh_edges_begin` An iterator pointing to the first element in an iterable decscribing the edges of the cells on which the Gauss-Legendre quadrature is calculated. 
* `mesh_edges_end` An iterator pointing to the end of an iterable decscribing the edges of the cells on which the Gauss-Legendre quadrature is calculated. 




        

<hr>



### function GaussLegendre [2/3]

_A constructor of the_ [_**GaussLegendre**_](classGaussLegendre.md) _class._
```C++
inline explicit GaussLegendre::GaussLegendre (
    std::initializer_list< Coord< Dim > > const mesh_edges
) 
```





**Parameters:**


* `mesh_edges` An initialiser list containing the edges of the cells on which the Gauss-Legendre quadrature is calculated. 




        

<hr>



### function GaussLegendre [3/3]

_A constructor of the_ [_**GaussLegendre**_](classGaussLegendre.md) _class._
```C++
template<class Grid1D>
inline explicit GaussLegendre::GaussLegendre (
    ConstField< Coord< Dim >, IdxRange< Grid1D >, Kokkos::HostSpace > mesh_edges
) 
```





**Parameters:**


* `mesh_edges` A constant Field containing the edges of the cells on which the Gauss-Legendre quadrature is calculated. 




        

<hr>



### function gauss\_legendre\_coefficients 

_Get a FieldMem containing the coefficients for the Gauss-Legendre quadrature._ 
```C++
template<class ExecSpace>
inline DFieldMem< IdxRange< GLGrid >, typename ExecSpace::memory_space > GaussLegendre::gauss_legendre_coefficients () const
```





**Returns:**

The Gauss-Legendre quadrature. 





        

<hr>



### function get\_idx\_range 

_Get the index range of the points of the Gauss-Legendre quadrature._ 
```C++
inline IdxRange< GLGrid > GaussLegendre::get_idx_range () const
```





**Returns:**

The index range where functions should be evaluated. 





        

<hr>
## Public Static Functions Documentation




### function order 

_The order of the quadrature scheme._ 
```C++
static inline constexpr std::size_t GaussLegendre::order () 
```





**Returns:**

The order of the quadrature scheme. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/gauss_legendre_integration.hpp`

