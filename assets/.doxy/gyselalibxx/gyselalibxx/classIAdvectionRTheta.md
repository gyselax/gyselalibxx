

# Class IAdvectionRTheta



[**ClassList**](annotated.md) **>** [**IAdvectionRTheta**](classIAdvectionRTheta.md)



_Define the base class of 2D advection operators in polar index range._ 

* `#include <iadvection_rtheta.hpp>`





Inherited by the following classes: [BslAdvectionRTheta](classBslAdvectionRTheta.md),  [BslAdvectionRTheta](classBslAdvectionRTheta.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DFieldRTheta | [**operator()**](#function-operator) (DFieldRTheta allfdistribu, [**DConstVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; advection\_field, double const dt) const = 0<br>_Advect a function along the advection field given on dt with a given advection field along XY._  |
| virtual DFieldRTheta | [**operator()**](#function-operator_1) (DFieldRTheta allfdistribu, [**DConstVectorFieldRTheta**](classVectorField.md)&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; advection\_field, [**DVector**](classTensor.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & advection\_field\_xy\_centre, double const dt) const = 0<br>_Advect a function along the advection field given on dt with a given advection field along RTheta._  |
| virtual  | [**~IAdvectionRTheta**](#function-iadvectionrtheta) () = default<br> |




























## Public Functions Documentation




### function operator() 

_Advect a function along the advection field given on dt with a given advection field along XY._ 
```C++
virtual DFieldRTheta IAdvectionRTheta::operator() (
    DFieldRTheta allfdistribu,
    DConstVectorFieldRTheta < X , Y > advection_field,
    double const dt
) const = 0
```





**Parameters:**


* `allfdistribu` The function to be advected. 
* `advection_field` The advection field along the physical index range axes, XY. 
* `dt` The time step.



**Returns:**

A Field of the advected function (allfdistribu). 





        

<hr>



### function operator() 

_Advect a function along the advection field given on dt with a given advection field along RTheta._ 
```C++
virtual DFieldRTheta IAdvectionRTheta::operator() (
    DFieldRTheta allfdistribu,
    DConstVectorFieldRTheta < R , Theta > advection_field,
    DVector < X , Y > const & advection_field_xy_centre,
    double const dt
) const = 0
```





**Parameters:**


* `allfdistribu` The function to be advected. 
* `advection_field` The advection field on the contravariant basis of the logical domain. 
* `advection_field_xy_centre` A vector in the Cartesian basis, containing the value of the advection field at the O-point. 
* `dt` The time step.



**Returns:**

A Field of the advected function (allfdistribu). 





        

<hr>



### function ~IAdvectionRTheta 

```C++
virtual IAdvectionRTheta::~IAdvectionRTheta () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/advection/iadvection_rtheta.hpp`

