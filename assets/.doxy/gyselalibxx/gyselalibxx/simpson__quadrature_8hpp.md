

# File simpson\_quadrature.hpp



[**FileList**](files.md) **>** [**quadrature**](dir_264321be3574e3b1cf375050e213576e.md) **>** [**simpson\_quadrature.hpp**](simpson__quadrature_8hpp.md)

[Go to the source code of this file](simpson__quadrature_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "geometry_descriptors.hpp"`
* `#include "quadrature_coeffs_nd.hpp"`
* `#include "trapezoid_quadrature.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**fill\_simpson\_quadrature\_coefficients\_1d**](#function-fill_simpson_quadrature_coefficients_1d) (DField&lt; IdxRange&lt; Grid1D &gt;, typename ExecSpace::memory\_space &gt; coefficients) <br>_Get the Simpson coefficients in 1D._  |
|  DFieldMem&lt; IdxRange&lt; ODims... &gt;, typename ExecSpace::memory\_space &gt; | [**simpson\_quadrature\_coefficients**](#function-simpson_quadrature_coefficients) (IdxRange&lt; ODims... &gt; const & idx\_range) <br>_Get the simpson coefficients in ND._  |
|  DFieldMem&lt; IdxRange&lt; Grid1D &gt;, typename ExecSpace::memory\_space &gt; | [**simpson\_quadrature\_coefficients\_1d**](#function-simpson_quadrature_coefficients_1d) (IdxRange&lt; Grid1D &gt; const & idx\_range) <br>_Get the Simpson coefficients in 1D._  |
|  DFieldMem&lt; IdxRange&lt; Grid1D &gt;, typename ExecSpace::memory\_space &gt; | [**simpson\_trapezoid\_quadrature\_coefficients\_1d**](#function-simpson_trapezoid_quadrature_coefficients_1d) (IdxRange&lt; Grid1D &gt; const & idx\_range, Extremity trapezoid\_extremity) <br>_Get the Simpson coefficients in 1D._  |




























## Public Functions Documentation




### function fill\_simpson\_quadrature\_coefficients\_1d 

_Get the Simpson coefficients in 1D._ 
```C++
template<class ExecSpace, class Grid1D>
void fill_simpson_quadrature_coefficients_1d (
    DField< IdxRange< Grid1D >, typename ExecSpace::memory_space > coefficients
) 
```



Calculate the quadrature coefficients for the Simpson method defined on the provided index range. The non-uniform form of the Simpson quadrature is:


(x\_3-x\_1)(2-(x\_3-x\_2)/(x\_2-x\_1)) / 6 f(x\_1) + (x\_3-x\_1)^3 / 6 / (x\_3-x\_2) / (x\_2-x\_1) f(x\_2) + (x\_3-x\_1)(2-(x\_2-x\_1)/(x\_3-x\_2)) / 6 f(x\_3)




**Parameters:**


* `coefficients` The field where the quadrature coefficients should be saved. 




        

<hr>



### function simpson\_quadrature\_coefficients 

_Get the simpson coefficients in ND._ 
```C++
template<class ExecSpace, class... ODims>
DFieldMem< IdxRange< ODims... >, typename ExecSpace::memory_space > simpson_quadrature_coefficients (
    IdxRange< ODims... > const & idx_range
) 
```



Calculate the quadrature coefficients for the simpson method defined on the provided index range.




**Template parameters:**


* `ExecSpace` Execution space, depends on Kokkos.



**Parameters:**


* `idx_range` The idx\_range on which the quadrature will be carried out.



**Returns:**

The quadrature coefficients for the trapezoid method defined on the provided idx\_range. The allocation place (host or device ) will depend on the ExecSpace. 





        

<hr>



### function simpson\_quadrature\_coefficients\_1d 

_Get the Simpson coefficients in 1D._ 
```C++
template<class ExecSpace, class Grid1D>
DFieldMem< IdxRange< Grid1D >, typename ExecSpace::memory_space > simpson_quadrature_coefficients_1d (
    IdxRange< Grid1D > const & idx_range
) 
```



Calculate the quadrature coefficients for the Simpson method defined on the provided index range. The non-uniform form of the Simpson quadrature is:


(x\_3-x\_1)(2-(x\_3-x\_2)/(x\_2-x\_1)) / 6 f(x\_1) + (x\_3-x\_1)^3 / 6 / (x\_3-x\_2) / (x\_2-x\_1) f(x\_2) + (x\_3-x\_1)(2-(x\_2-x\_1)/(x\_3-x\_2)) / 6 f(x\_3)




**Parameters:**


* `idx_range` The index range on which the quadrature will be carried out.



**Returns:**

The quadrature coefficients for the Simpson method defined on the provided index range. 





        

<hr>



### function simpson\_trapezoid\_quadrature\_coefficients\_1d 

_Get the Simpson coefficients in 1D._ 
```C++
template<class ExecSpace, class Grid1D>
DFieldMem< IdxRange< Grid1D >, typename ExecSpace::memory_space > simpson_trapezoid_quadrature_coefficients_1d (
    IdxRange< Grid1D > const & idx_range,
    Extremity trapezoid_extremity
) 
```



If the number of grid points is not compatible with the Simpson quadrature scheme then use a trapezoid formula over one cell at the specified extremity.




**Parameters:**


* `idx_range` The index range on which the quadrature will be carried out. 
* `trapezoid_extremity` The extremity where the trapezoid quadrature may be used.



**Returns:**

The quadrature coefficients for the Simpson method defined on the provided index range. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/quadrature/simpson_quadrature.hpp`

