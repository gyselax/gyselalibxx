

# Class BslAdvectionVelocity

**template &lt;class Geometry, class GridV&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvectionVelocity**](classBslAdvectionVelocity.md)



_A class which computes the velocity advection along the dimension of interest GridV. Working for every Cartesian geometry._ 

* `#include <bsl_advection_vx.hpp>`



Inherits the following classes: [IAdvectionVelocity](classIAdvectionVelocity.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvectionVelocity**](#function-bsladvectionvelocity) (PreallocatableInterpolatorType const & interpolator\_v) <br>_Constructor._  |
|  Field&lt; double, IdxRangeFdistribu &gt; | [**operator()**](#function-operator) (Field&lt; double, IdxRangeFdistribu &gt; const allfdistribu, Field&lt; const double, IdxRangeSpatial &gt; const electric\_field, double const dt) override const<br>_Advects fdistribu along GridV for a duration dt._  |
|   | [**~BslAdvectionVelocity**](#function-bsladvectionvelocity) () override<br> |


## Public Functions inherited from IAdvectionVelocity

See [IAdvectionVelocity](classIAdvectionVelocity.md)

| Type | Name |
| ---: | :--- |
| virtual DField&lt; typename Geometry::IdxRangeFdistribu &gt; | [**operator()**](classIAdvectionVelocity.md#function-operator) (DField&lt; typename Geometry::IdxRangeFdistribu &gt; allfdistribu, DConstField&lt; typename Geometry::IdxRangeSpatial &gt; electric\_field, double dt) const = 0<br>_operates a transport of the distribution function._  |
| virtual  | [**~IAdvectionVelocity**](classIAdvectionVelocity.md#function-iadvectionvelocity) () = default<br> |






















































## Public Functions Documentation




### function BslAdvectionVelocity 

_Constructor._ 
```C++
inline explicit BslAdvectionVelocity::BslAdvectionVelocity (
    PreallocatableInterpolatorType const & interpolator_v
) 
```





**Parameters:**


* `interpolator_v` interpolator along the GridV direction which refers to the velocity space. 
 




        

<hr>



### function operator() 

_Advects fdistribu along GridV for a duration dt._ 
```C++
inline Field< double, IdxRangeFdistribu > BslAdvectionVelocity::operator() (
    Field< double, IdxRangeFdistribu > const allfdistribu,
    Field< const double, IdxRangeSpatial > const electric_field,
    double const dt
) override const
```





**Parameters:**


* `allfdistribu` Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration). 
* `electric_field` Reference to the electric field which derives from electrostatic potential, allocated on the device. 
* `dt` Time step 



**Returns:**

A reference to the allfdistribu array containing the value of the function at the coordinates. 





        

<hr>



### function ~BslAdvectionVelocity 

```C++
BslAdvectionVelocity::~BslAdvectionVelocity () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/bsl_advection_vx.hpp`

