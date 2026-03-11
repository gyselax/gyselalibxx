

# Class NonUniformLagrangeBasis

**template &lt;class Dim, std::size\_t D, class DataType&gt;**



[**ClassList**](annotated.md) **>** [**NonUniformLagrangeBasis**](classNonUniformLagrangeBasis.md)



_Class describing Lagrange polynomials on a uniform grid._ [More...](#detailed-description)

* `#include <lagrange_basis_non_uniform.hpp>`



Inherits the following classes: detail::NonUniformLagrangeBasisBase












## Classes

| Type | Name |
| ---: | :--- |
| class | [**Impl**](classNonUniformLagrangeBasis_1_1Impl.md) &lt;class DDim, class MemorySpace&gt;<br>_Storage class of the static attributes of the discrete dimension._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef Dim | [**continuous\_dimension\_type**](#typedef-continuous_dimension_type)  <br>_The tag identifying the continuous dimension on which the Lagrange polynomials are defined._  |
| typedef Coord&lt; Dim &gt; | [**coord\_type**](#typedef-coord_type)  <br>_The type of the coordinates on which the Lagrange polynomials can be evaluated._  |
| typedef [**NonUniformLagrangeBasis**](classNonUniformLagrangeBasis.md) | [**discrete\_dimension\_type**](#typedef-discrete_dimension_type)  <br>_The discrete dimension representing B-splines._  |






















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**degree**](#function-degree) () noexcept<br>_The degree of the Lagrange polynomials._  |
|  constexpr bool | [**is\_periodic**](#function-is_periodic) () noexcept<br>_Indicates if the Lagrange polynomials are periodic or not._  |
|  constexpr bool | [**is\_uniform**](#function-is_uniform) () noexcept<br>_Indicates if the Lagrange polynomials are uniform or not (this is the case here)._  |


























## Detailed Description


This class uses the second barycentric formulation to evaluate the polynomials. This formula is used for stability. It is described in Barycentric Lagrange Interpolation Jean-Paul Berrut and Lloyd N. Trefethen SIAM Review 2004 46:3, 501-517




**Template parameters:**


* `Dim` The dimension on which the Lagrange polynomials are defined. 
* `D` The degree of the polynomials, equal to the number of cells over which the Lagrange polynomials are defined. 
* `DataType` The data type used for the calculations. Double by default. 




    
## Public Types Documentation




### typedef continuous\_dimension\_type 

_The tag identifying the continuous dimension on which the Lagrange polynomials are defined._ 
```C++
using NonUniformLagrangeBasis< Dim, D, DataType >::continuous_dimension_type =  Dim;
```




<hr>



### typedef coord\_type 

_The type of the coordinates on which the Lagrange polynomials can be evaluated._ 
```C++
using NonUniformLagrangeBasis< Dim, D, DataType >::coord_type =  Coord<Dim>;
```




<hr>



### typedef discrete\_dimension\_type 

_The discrete dimension representing B-splines._ 
```C++
using NonUniformLagrangeBasis< Dim, D, DataType >::discrete_dimension_type =  NonUniformLagrangeBasis;
```




<hr>
## Public Static Functions Documentation




### function degree 

_The degree of the Lagrange polynomials._ 
```C++
static inline constexpr std::size_t NonUniformLagrangeBasis::degree () noexcept
```





**Returns:**

The degree. 





        

<hr>



### function is\_periodic 

_Indicates if the Lagrange polynomials are periodic or not._ 
```C++
static inline constexpr bool NonUniformLagrangeBasis::is_periodic () noexcept
```





**Returns:**

A boolean indicating if the Lagrange polynomials are periodic or not. 





        

<hr>



### function is\_uniform 

_Indicates if the Lagrange polynomials are uniform or not (this is the case here)._ 
```C++
static inline constexpr bool NonUniformLagrangeBasis::is_uniform () noexcept
```





**Returns:**

A boolean indicating if the Lagrange polynomials are uniform or not. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/lagrange_basis_non_uniform.hpp`

