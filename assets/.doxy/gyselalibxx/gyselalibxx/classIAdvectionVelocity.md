

# Class IAdvectionVelocity

**template &lt;class Geometry, class GridV&gt;**



[**ClassList**](annotated.md) **>** [**IAdvectionVelocity**](classIAdvectionVelocity.md)



_A class which provides an advection operator._ [More...](#detailed-description)

* `#include <iadvectionvx.hpp>`





Inherited by the following classes: [BslAdvectionVelocity](classBslAdvectionVelocity.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DField&lt; typename Geometry::IdxRangeFdistribu &gt; | [**operator()**](#function-operator) (DField&lt; typename Geometry::IdxRangeFdistribu &gt; allfdistribu, DConstField&lt; typename Geometry::IdxRangeSpatial &gt; electric\_field, double dt) const = 0<br>_operates a transport of the distribution function._  |
| virtual  | [**~IAdvectionVelocity**](#function-iadvectionvelocity) () = default<br> |




























## Detailed Description


An abstract class which implements a function that applies the transport along a velocity direction of the phase space.


A generic class for a velocity advection 


    
## Public Functions Documentation




### function operator() 

_operates a transport of the distribution function._ 
```C++
virtual DField< typename Geometry::IdxRangeFdistribu > IAdvectionVelocity::operator() (
    DField< typename Geometry::IdxRangeFdistribu > allfdistribu,
    DConstField< typename Geometry::IdxRangeSpatial > electric_field,
    double dt
) const = 0
```





**Parameters:**


* `allfdistribu` Reference to an array containing the value of the distribution function. 
* `electric_field` The electric field which derives from electrostatic potential and is the advection speed. 
* `dt` Time step.



**Returns:**

A reference to an array containing the value of distribution the function at the updated time t+dt. 





        

<hr>



### function ~IAdvectionVelocity 

```C++
virtual IAdvectionVelocity::~IAdvectionVelocity () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/iadvectionvx.hpp`

