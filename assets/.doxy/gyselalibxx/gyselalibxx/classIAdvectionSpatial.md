

# Class IAdvectionSpatial

**template &lt;class Geometry, class [**GridX**](structGridX.md)&gt;**



[**ClassList**](annotated.md) **>** [**IAdvectionSpatial**](classIAdvectionSpatial.md)



_A class which provides an advection operator._ [More...](#detailed-description)

* `#include <iadvectionx.hpp>`





Inherited by the following classes: [BslAdvectionSpatial](classBslAdvectionSpatial.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DField&lt; typename Geometry::IdxRangeFdistribu &gt; | [**operator()**](#function-operator) (DField&lt; typename Geometry::IdxRangeFdistribu &gt; allfdistribu, double dt) const = 0<br>_operates a transport of the distribution function._  |
| virtual  | [**~IAdvectionSpatial**](#function-iadvectionspatial) () = default<br> |




























## Detailed Description


An abstract class which implements a function that applies the transport along a physical space direction of the phase space.


A generic class for a spatial advection 


    
## Public Functions Documentation




### function operator() 

_operates a transport of the distribution function._ 
```C++
virtual DField< typename Geometry::IdxRangeFdistribu > IAdvectionSpatial::operator() (
    DField< typename Geometry::IdxRangeFdistribu > allfdistribu,
    double dt
) const = 0
```





**Parameters:**


* `allfdistribu` reference to an array containing the value of distribution the function. 
* `dt` time step.



**Returns:**

A reference to an array containing the value of distribution the function at the updated time T+dt. 





        

<hr>



### function ~IAdvectionSpatial 

```C++
virtual IAdvectionSpatial::~IAdvectionSpatial () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/iadvectionx.hpp`

