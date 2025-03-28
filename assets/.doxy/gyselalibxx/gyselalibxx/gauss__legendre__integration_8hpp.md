

# File gauss\_legendre\_integration.hpp



[**FileList**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**gauss\_legendre\_integration.hpp**](gauss__legendre__integration_8hpp.md)

[Go to the source code of this file](gauss__legendre__integration_8hpp_source.md)



* `#include <array>`
* `#include <vector>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**GaussLegendre**](classGaussLegendre.md) &lt;class GLGrid, NPoints&gt;<br>_An operator for constructing a Gauss-Legendre quadrature._  |
| struct | [**GaussLegendreCoefficients**](structGaussLegendreCoefficients.md) &lt;NPoints&gt;<br>_A structure containing the weights and positions associated with a Gauss-Legendre quadrature using NPoints points._  |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  DFieldMem&lt; IdxRange&lt; typename GaussLegendreQuad::Grid1D... &gt;, typename ExecSpace::memory\_space &gt; | [**gauss\_legendre\_quadrature\_coefficients**](#function-gauss_legendre_quadrature_coefficients) (GaussLegendreQuad const &... gl) <br>_Get the spline quadrature coefficients in ND from N 1D quadrature coefficient._  |




























## Public Functions Documentation




### function gauss\_legendre\_quadrature\_coefficients 

_Get the spline quadrature coefficients in ND from N 1D quadrature coefficient._ 
```C++
template<class ExecSpace, class... GaussLegendreQuad>
DFieldMem< IdxRange< typename GaussLegendreQuad::Grid1D... >, typename ExecSpace::memory_space > gauss_legendre_quadrature_coefficients (
    GaussLegendreQuad const &... gl
) 
```



Calculate the quadrature coefficients for the spline quadrature method defined on the provided index range.




**Parameters:**


* `gl` The Gauss-Legendre quadrature objects describing each dimension.



**Returns:**

The coefficients which define the spline quadrature method in ND. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/gauss_legendre_integration.hpp`

