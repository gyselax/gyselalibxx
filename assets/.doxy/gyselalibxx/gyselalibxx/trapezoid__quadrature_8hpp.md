

# File trapezoid\_quadrature.hpp



[**FileList**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**trapezoid\_quadrature.hpp**](trapezoid__quadrature_8hpp.md)

[Go to the source code of this file](trapezoid__quadrature_8hpp_source.md)

[More...](#detailed-description)

* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "quadrature_coeffs_nd.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  DFieldMem&lt; IdxRange&lt; ODims... &gt;, typename ExecSpace::memory\_space &gt; | [**trapezoid\_quadrature\_coefficients**](#function-trapezoid_quadrature_coefficients) (IdxRange&lt; ODims... &gt; const & idx\_range) <br>_Get the trapezoid coefficients in ND._  |
|  DFieldMem&lt; IdxRange&lt; Grid1D &gt;, typename ExecSpace::memory\_space &gt; | [**trapezoid\_quadrature\_coefficients\_1d**](#function-trapezoid_quadrature_coefficients_1d) (IdxRange&lt; Grid1D &gt; const & idx\_range) <br>_Get the trapezoid coefficients in 1D._  |




























## Detailed Description


File providing quadrature coefficients via the trapezoidal method. 


    
## Public Functions Documentation




### function trapezoid\_quadrature\_coefficients 

_Get the trapezoid coefficients in ND._ 
```C++
template<class ExecSpace, class... ODims>
DFieldMem< IdxRange< ODims... >, typename ExecSpace::memory_space > trapezoid_quadrature_coefficients (
    IdxRange< ODims... > const & idx_range
) 
```



Calculate the quadrature coefficients for the trapezoid method defined on the provided index range.




**Template parameters:**


* `ExecSpace` Execution space, depends on Kokkos.



**Parameters:**


* `idx_range` The idx\_range on which the quadrature will be carried out.



**Returns:**

The quadrature coefficients for the trapezoid method defined on the provided idx\_range. The allocation place (host or device ) will depend on the ExecSpace. 





        

<hr>



### function trapezoid\_quadrature\_coefficients\_1d 

_Get the trapezoid coefficients in 1D._ 
```C++
template<class ExecSpace, class Grid1D>
DFieldMem< IdxRange< Grid1D >, typename ExecSpace::memory_space > trapezoid_quadrature_coefficients_1d (
    IdxRange< Grid1D > const & idx_range
) 
```



Calculate the quadrature coefficients for the trapezoid method defined on the provided index range.




**Template parameters:**


* `ExecSpace` Execution space, depends on Kokkos.



**Parameters:**


* `idx_range` The idx\_range on which the quadrature will be carried out.



**Returns:**

The quadrature coefficients for the trapezoid method defined on the provided index range. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/trapezoid_quadrature.hpp`

