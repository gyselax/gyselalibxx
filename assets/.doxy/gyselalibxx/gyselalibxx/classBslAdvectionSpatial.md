

# Class BslAdvectionSpatial

**template &lt;class Geometry, class [**GridX**](structGridX.md)&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvectionSpatial**](classBslAdvectionSpatial.md)



_A class which computes the spatial advection along the dimension of interest_ [_**GridX**_](structGridX.md) _. Working for every Cartesian geometry._

* `#include <bsl_advection_x.hpp>`



Inherits the following classes: [IAdvectionSpatial](classIAdvectionSpatial.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvectionSpatial**](#function-bsladvectionspatial) (PreallocatableInterpolatorType const & interpolator\_x) <br>_Constructor_  __ |
|  Field&lt; double, IdxRangeFdistrib &gt; | [**operator()**](#function-operator) (Field&lt; double, IdxRangeFdistrib &gt; const allfdistribu, double const dt) override const<br>_Advects fdistribu along_ [_**GridX**_](structGridX.md) _for a duration dt._ |
|   | [**~BslAdvectionSpatial**](#function-bsladvectionspatial) () override<br> |


## Public Functions inherited from IAdvectionSpatial

See [IAdvectionSpatial](classIAdvectionSpatial.md)

| Type | Name |
| ---: | :--- |
| virtual DField&lt; typename Geometry::IdxRangeFdistribu &gt; | [**operator()**](classIAdvectionSpatial.md#function-operator) (DField&lt; typename Geometry::IdxRangeFdistribu &gt; allfdistribu, double dt) const = 0<br>_operates a transport of the distribution function._  |
| virtual  | [**~IAdvectionSpatial**](classIAdvectionSpatial.md#function-iadvectionspatial) () = default<br> |






















































## Public Functions Documentation




### function BslAdvectionSpatial 

_Constructor_  __
```C++
inline explicit BslAdvectionSpatial::BslAdvectionSpatial (
    PreallocatableInterpolatorType const & interpolator_x
) 
```





**Parameters:**


* `interpolator_x` interpolator along the [**GridX**](structGridX.md) direction which refers to the spatial space. 
 




        

<hr>



### function operator() 

_Advects fdistribu along_ [_**GridX**_](structGridX.md) _for a duration dt._
```C++
inline Field< double, IdxRangeFdistrib > BslAdvectionSpatial::operator() (
    Field< double, IdxRangeFdistrib > const allfdistribu,
    double const dt
) override const
```





**Parameters:**


* `allfdistribu` Reference to the whole distribution function for one species, allocated on the device (ie it lets the choice of the location depend on the build configuration). 
* `dt` Time step 



**Returns:**

A reference to the allfdistribu array containing the value of the function at the coordinates. 





        

<hr>



### function ~BslAdvectionSpatial 

```C++
BslAdvectionSpatial::~BslAdvectionSpatial () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/bsl_advection_x.hpp`

