

# Class Lagrange

**template &lt;class Execspace, class GridInterp, BCond BcMin, BCond BcMax&gt;**



[**ClassList**](annotated.md) **>** [**Lagrange**](classLagrange.md)



_A class which implements_ [_**Lagrange**_](classLagrange.md) _polynomials._[More...](#detailed-description)

* `#include <Lagrange.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**Lagrange**](#function-lagrange) (int degree, Field&lt; double, IdxRangeInterp, typename Execspace::memory\_space &gt; x\_nodes\_fnodes, IdxRangeInterp idx\_range, IdxStepInterp ghost) <br>_Usual Constructor._  |
|  KOKKOS\_FUNCTION double | [**evaluate**](#function-evaluate) (CoordDimI x\_interp) const<br>_Evaluates the approximated value of a function on a point current values at a known set of interpolation points._  |




























## Detailed Description


A simple class which provides the possibility to evaluate an interpolation of a function known only over a restricted set of nodes. 


    
## Public Functions Documentation




### function Lagrange 

_Usual Constructor._ 
```C++
inline KOKKOS_FUNCTION Lagrange::Lagrange (
    int degree,
    Field< double, IdxRangeInterp, typename Execspace::memory_space > x_nodes_fnodes,
    IdxRangeInterp idx_range,
    IdxStepInterp ghost
) 
```





**Parameters:**


* `degree` integer which correspond to the degree of interpolation. 
* `x_nodes_fnodes` Field of nodes and associated values of the function. 
* `idx_range` along interest direction, usedful in periodic case 
* `ghost` DiscretVector which gives the number of ghosted points 




        

<hr>



### function evaluate 

_Evaluates the approximated value of a function on a point current values at a known set of interpolation points._ 
```C++
KOKKOS_FUNCTION double Lagrange::evaluate (
    CoordDimI x_interp
) const
```





**Parameters:**


* `x_interp` a node where we want to evaluate the function.



**Returns:**

The evaluation of [**Lagrange**](classLagrange.md) interpolation at the point x\_intercept. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/Lagrange.hpp`

