

# Class TriangularBernsteinPolynomialBasis

**template &lt;class [**X**](structX.md), class [**Y**](structY.md), class Corner1Tag, class Corner2Tag, class Corner3Tag, std::size\_t D&gt;**



[**ClassList**](annotated.md) **>** [**TriangularBernsteinPolynomialBasis**](classTriangularBernsteinPolynomialBasis.md)



_A class which evaluates the triangular Bernstein polynomials._ [More...](#detailed-description)

* `#include <bernstein.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**Impl**](classTriangularBernsteinPolynomialBasis_1_1Impl.md) &lt;class DDim, class MemorySpace&gt;<br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**TriangularBernsteinPolynomialBasis**](classTriangularBernsteinPolynomialBasis.md) | [**discrete\_dimension\_type**](#typedef-discrete_dimension_type)  <br>_A tag for DDC to recognise the discrete dimension type._  |






















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**degree**](#function-degree) () noexcept<br>_The degree of the triangular Bernstein polynomials._  |
|  constexpr std::size\_t | [**nbasis**](#function-nbasis) () <br>_The number of basis elements._  |
|  constexpr std::size\_t | [**rank**](#function-rank) () <br>_The rank of the system of equations. This is equal to the number of dimensions in the coordinates on which the polynomials are evaluated._  |


























## Detailed Description


Triangular Bernstein polynomials of degree  are defined as: 


Where , , and  are the barycentric coordinates of 


c.f. Toshniwal et al. 2017 Multi-degree smooth polar splines: A framework for geometric modelling and isogeometric analysis [https://doi.org/10.1016/j.cma.2016.11.009](https://doi.org/10.1016/j.cma.2016.11.009)




**Template parameters:**


* [**X**](structX.md) The first dimension of the Cartesian coordinate system. 
* [**Y**](structY.md) The second dimension of the Cartesian coordinate system. 
* `Corner1Tag` A class to identify the first corner. 
* `Corner2Tag` A class to identify the second corner. 
* `Corner3Tag` A class to identify the third corner. 
* `D` The degree of the polynomial. 




    
## Public Types Documentation




### typedef discrete\_dimension\_type 

_A tag for DDC to recognise the discrete dimension type._ 
```C++
using TriangularBernsteinPolynomialBasis< X, Y, Corner1Tag, Corner2Tag, Corner3Tag, D >::discrete_dimension_type =  TriangularBernsteinPolynomialBasis;
```




<hr>
## Public Static Functions Documentation




### function degree 

_The degree of the triangular Bernstein polynomials._ 
```C++
static inline constexpr std::size_t TriangularBernsteinPolynomialBasis::degree () noexcept
```





**Returns:**

The degree of the triangular Bernstein polynomials. 





        

<hr>



### function nbasis 

_The number of basis elements._ 
```C++
static inline constexpr std::size_t TriangularBernsteinPolynomialBasis::nbasis () 
```





**Returns:**

The number of basis elements. 





        

<hr>



### function rank 

_The rank of the system of equations. This is equal to the number of dimensions in the coordinates on which the polynomials are evaluated._ 
```C++
static inline constexpr std::size_t TriangularBernsteinPolynomialBasis::rank () 
```





**Returns:**

The rank of the system. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/bernstein.hpp`

