

# File volume\_quadrature\_nd.hpp



[**FileList**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**volume\_quadrature\_nd.hpp**](volume__quadrature__nd_8hpp.md)

[Go to the source code of this file](volume__quadrature__nd_8hpp_source.md)

[More...](#detailed-description)

* `#include <ddc/ddc.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "mapping_tools.hpp"`
* `#include "quadrature.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  DFieldMem&lt; IdxRangeCoeffs, typename ExecSpace::memory\_space &gt; | [**compute\_coeffs\_on\_mapping**](#function-compute_coeffs_on_mapping) (ExecSpace exec\_space, Mapping & mapping, DFieldMem&lt; IdxRangeCoeffs, typename ExecSpace::memory\_space &gt; && coefficients\_alloc) <br>_Add the Jacobian determinant to the coefficients._  |




























## Detailed Description


File providing functions to adapt quadrature coefficients so they can integrate over a N-D volume. A N-D volume is a surface in 2D and a volume in 3D. 


    
## Public Functions Documentation




### function compute\_coeffs\_on\_mapping 

_Add the Jacobian determinant to the coefficients._ 
```C++
template<class Mapping, class IdxRangeCoeffs, class ExecSpace>
DFieldMem< IdxRangeCoeffs, typename ExecSpace::memory_space > compute_coeffs_on_mapping (
    ExecSpace exec_space,
    Mapping & mapping,
    DFieldMem< IdxRangeCoeffs, typename ExecSpace::memory_space > && coefficients_alloc
) 
```



For polar integration, we can add the Jacobian determinant to the quadrature coefficients.



* 2D example: (but it is implemented for ND)







This function uses rvalues. It means that coefficients is a temporary input parameter and it returns a temporary coefficient object. The [**Quadrature**](classQuadrature.md) object can only be instantiate with rvalues.




**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `mapping` The mapping function from the logical index range  to the physical index range . 
* `coefficients_alloc` The quadrature coefficients .



**Returns:**

A rvalue FieldMem to the modified coefficients . 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/volume_quadrature_nd.hpp`

