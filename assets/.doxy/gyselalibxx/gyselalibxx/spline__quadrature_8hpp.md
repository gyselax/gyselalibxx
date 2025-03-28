

# File spline\_quadrature.hpp



[**FileList**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**spline\_quadrature.hpp**](spline__quadrature_8hpp.md)

[Go to the source code of this file](spline__quadrature_8hpp_source.md)

[More...](#detailed-description)

* `#include <cassert>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  DFieldMem&lt; IdxRange&lt; DDims... &gt;, typename ExecSpace::memory\_space &gt; | [**spline\_quadrature\_coefficients**](#function-spline_quadrature_coefficients) (IdxRange&lt; DDims... &gt; const & idx\_range, SplineBuilders const &... builders) <br>_Get the spline quadrature coefficients in ND from N 1D quadrature coefficient._  |
|  host\_t&lt; DFieldMem&lt; IdxRange&lt; Grid1D &gt; &gt; &gt; | [**spline\_quadrature\_coefficients\_1d**](#function-spline_quadrature_coefficients_1d) (IdxRange&lt; Grid1D &gt; const & idx\_range, SplineBuilder const & builder) <br>_Get the spline quadrature coefficients._  |




























## Detailed Description


File providing quadrature coefficients via a spline quadrature. 


    
## Public Functions Documentation




### function spline\_quadrature\_coefficients 

_Get the spline quadrature coefficients in ND from N 1D quadrature coefficient._ 
```C++
template<class ExecSpace, class... DDims, class... SplineBuilders>
DFieldMem< IdxRange< DDims... >, typename ExecSpace::memory_space > spline_quadrature_coefficients (
    IdxRange< DDims... > const & idx_range,
    SplineBuilders const &... builders
) 
```



Calculate the quadrature coefficients for the spline quadrature method defined on the provided index range.




**Parameters:**


* `idx_range` The index range on which the coefficients will be defined. 
* `builders` The spline builder used for the quadrature coefficients in the different dimensions.



**Returns:**

The coefficients which define the spline quadrature method in ND. 





        

<hr>



### function spline\_quadrature\_coefficients\_1d 

_Get the spline quadrature coefficients._ 
```C++
template<class Grid1D, class SplineBuilder>
host_t< DFieldMem< IdxRange< Grid1D > > > spline_quadrature_coefficients_1d (
    IdxRange< Grid1D > const & idx_range,
    SplineBuilder const & builder
) 
```



To integrate a function with a spline quadrature, we use:


,


which rewritten gives


,


with
*  the values of the function at the interpolation points;
*  the quadrature coefficients we compute thanks to ,
  * with  the matrix of B-splines ,
  * and  the integrated B-splines.






More details are given in Emily Bourne's thesis "Non-Uniform Numerical Schemes for the Modelling of Turbulence
in the 5D GYSELA Code". December 2022.




**Parameters:**


* `idx_range` The index range where the functions we want to integrate are defined. 
* `builder` The spline builder describing the way in which splines would be constructed.



**Returns:**

A chunk with the quadrature coefficients . 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/spline_quadrature.hpp`

