

# File neumann\_spline\_quadrature.hpp



[**FileList**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**neumann\_spline\_quadrature.hpp**](neumann__spline__quadrature_8hpp.md)

[Go to the source code of this file](neumann__spline__quadrature_8hpp_source.md)

[More...](#detailed-description)

* `#include <cassert>`
* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`
* `#include "ddc_aliases.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  DFieldMem&lt; IdxRange&lt; DDims... &gt;, typename ExecSpace::memory\_space &gt; | [**neumann\_spline\_quadrature\_coefficients**](#function-neumann_spline_quadrature_coefficients) (IdxRange&lt; DDims... &gt; const & idx\_range, SplineBuilders const &... builders) <br>_Get the spline quadrature coefficients in ND from N 1D quadrature coefficient._  |
|  DFieldMem&lt; IdxRange&lt; Grid1D &gt;, typename ExecSpace::memory\_space &gt; | [**neumann\_spline\_quadrature\_coefficients\_1d**](#function-neumann_spline_quadrature_coefficients_1d) (IdxRange&lt; Grid1D &gt; const & idx\_range, SplineBuilder const & builder) <br>_Get the spline quadrature coefficients in 1D._  |




























## Detailed Description


File providing quadrature coefficients via a spline quadrature. 


    
## Public Functions Documentation




### function neumann\_spline\_quadrature\_coefficients 

_Get the spline quadrature coefficients in ND from N 1D quadrature coefficient._ 
```C++
template<class ExecSpace, class... DDims, class... SplineBuilders>
DFieldMem< IdxRange< DDims... >, typename ExecSpace::memory_space > neumann_spline_quadrature_coefficients (
    IdxRange< DDims... > const & idx_range,
    SplineBuilders const &... builders
) 
```



This function calculates the quadrature coefficients which define a quadrature equivalent to calculating and integrating a spline approximation of a function. The spline approximation would be calculated with homogeneous Neumann boundary conditions.




**Parameters:**


* `idx_range` The index range on which the coefficients will be defined. 
* `builders` The spline builder used for the quadrature coefficients in the different dimensions.



**Returns:**

The coefficients which define the spline quadrature method in ND. 





        

<hr>



### function neumann\_spline\_quadrature\_coefficients\_1d 

_Get the spline quadrature coefficients in 1D._ 
```C++
template<class ExecSpace, class Grid1D, class SplineBuilder>
DFieldMem< IdxRange< Grid1D >, typename ExecSpace::memory_space > neumann_spline_quadrature_coefficients_1d (
    IdxRange< Grid1D > const & idx_range,
    SplineBuilder const & builder
) 
```



This function calculates the quadrature coefficients which define a quadrature equivalent to calculating and integrating a spline approximation of a function. The spline approximation would be calculated with homogeneous Neumann boundary conditions. This method of defining quadrature coefficients is described in section Emily Bourne's thesis[1].


[1] Non-Uniform Numerical Schemes for the Modelling of Turbulence in the 5D GYSELA Code Emily Bourne, December 2022-




**Parameters:**


* `idx_range` The index range on which the splines quadrature will be carried out. 
* `builder` The spline builder used for the quadrature coefficients.



**Returns:**

The quadrature coefficients for the method defined on the provided index range. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/neumann_spline_quadrature.hpp`

