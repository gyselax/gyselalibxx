

# File quadrature\_coeffs\_nd.hpp



[**FileList**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**quadrature\_coeffs\_nd.hpp**](quadrature__coeffs__nd_8hpp.md)

[Go to the source code of this file](quadrature__coeffs__nd_8hpp_source.md)

[More...](#detailed-description)

* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  DFieldMem&lt; IdxRange&lt; DDims... &gt;, typename ExecSpace::memory\_space &gt; | [**quadrature\_coeffs\_nd**](#function-quadrature_coeffs_nd) (IdxRange&lt; DDims... &gt; const & idx\_range, std::function&lt; DFieldMem&lt; IdxRange&lt; DDims &gt;, typename ExecSpace::memory\_space &gt;(IdxRange&lt; DDims &gt;)&gt;... funcs) <br>_Helper function which creates ND dimensions from N 1D quadrature coefficient functions._  |




























## Detailed Description


File providing helper functions for defining multi-dimensional quadrature methods. 


    
## Public Functions Documentation




### function quadrature\_coeffs\_nd 

_Helper function which creates ND dimensions from N 1D quadrature coefficient functions._ 
```C++
template<class ExecSpace, class... DDims>
DFieldMem< IdxRange< DDims... >, typename ExecSpace::memory_space > quadrature_coeffs_nd (
    IdxRange< DDims... > const & idx_range,
    std::function< DFieldMem< IdxRange< DDims >, typename ExecSpace::memory_space >(IdxRange< DDims >)>... funcs
) 
```





**Parameters:**


* `idx_range` The index range on which the coefficients will be defined. 
* `funcs` The functions which define quadrature coefficients in the different dimensions.



**Returns:**

The coefficients which define the quadrature method in ND. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/quadrature_coeffs_nd.hpp`

