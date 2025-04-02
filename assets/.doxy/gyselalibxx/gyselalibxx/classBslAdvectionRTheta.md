

# Class BslAdvectionRTheta

**template &lt;class FootFinder, class Mapping&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvectionRTheta**](classBslAdvectionRTheta.md)



_Define an advection operator on 2D_  _index range._[More...](#detailed-description)

* `#include <bsl_advection_rtheta.hpp>`



Inherits the following classes: [IAdvectionRTheta](classIAdvectionRTheta.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvectionRTheta**](#function-bsladvectionrtheta) ([**PreallocatableSplineInterpolatorType**](classPreallocatableSplineInterpolator2D.md) const & function\_interpolator, FootFinder const & foot\_finder, Mapping const & mapping) <br>_Instantiate an advection operator._  |
| virtual DFieldRTheta | [**operator()**](#function-operator) (DFieldRTheta allfdistribu, [**DConstVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; advection\_field\_xy, double dt) override const<br>_Allocate a Field of the advected function._  |
| virtual DFieldRTheta | [**operator()**](#function-operator_1) (DFieldRTheta allfdistribu, [**DConstVectorFieldRTheta**](classVectorField.md)&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; advection\_field\_rtheta, [**DVector**](classTensor.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & advection\_field\_xy\_centre, double dt) override const<br>_Allocate a Field to the advected function._  |
|   | [**~BslAdvectionRTheta**](#function-bsladvectionrtheta) () override<br> |


## Public Functions inherited from IAdvectionRTheta

See [IAdvectionRTheta](classIAdvectionRTheta.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldRTheta | [**operator()**](classIAdvectionRTheta.md#function-operator) (DFieldRTheta allfdistribu, [**DConstVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; advection\_field, double const dt) const = 0<br>_Advect a function along the advection field given on dt with a given advection field along XY._  |
| virtual DFieldRTheta | [**operator()**](classIAdvectionRTheta.md#function-operator_1) (DFieldRTheta allfdistribu, [**DConstVectorFieldRTheta**](classVectorField.md)&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; advection\_field, [**DVector**](classTensor.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; const & advection\_field\_xy\_centre, double const dt) const = 0<br>_Advect a function along the advection field given on dt with a given advection field along RTheta._  |
| virtual  | [**~IAdvectionRTheta**](classIAdvectionRTheta.md#function-iadvectionrtheta) () = default<br> |






















































## Detailed Description


The advection operator uses a backward semi-Lagrangian method. The method is based on the property that the solution is constant along the characteristics.


For the following equation: 


we write the characteristics: 


Then the property gives us: 


So the first step of the advection operator is to compute the feet of the characteristics  for each mesh point .


For the second step, we interpolate the function at the computed feet of the characteristics, and obtain the function at the next time step: .


Different time integration methods are implemented to solve the equation of the characteristics. They are defined in the [**IPolarFootFinder**](classIPolarFootFinder.md) class.


The feet can be advected on different domains (physical domain or pseudo-physical domain) which are determined in the [**SplinePolarFootFinder**](classSplinePolarFootFinder.md) operator.


The interpolation of the function is always done in the logical domain, where the B-splines are defined.




**See also:** [**IPolarFootFinder**](classIPolarFootFinder.md) 



    
## Public Functions Documentation




### function BslAdvectionRTheta 

_Instantiate an advection operator._ 
```C++
inline BslAdvectionRTheta::BslAdvectionRTheta (
    PreallocatableSplineInterpolatorType const & function_interpolator,
    FootFinder const & foot_finder,
    Mapping const & mapping
) 
```





**Parameters:**


* `function_interpolator` The polar interpolator to interpolate the function once the characteristics have been computed. 
* `foot_finder` An IFootFinder which computes the feet of the characteristics. 
* `mapping` The mapping function from the logical domain to the physical domain.



**Template parameters:**


* `IFootFinder` A child class of IFootFinder. 




        

<hr>



### function operator() 

_Allocate a Field of the advected function._ 
```C++
inline virtual DFieldRTheta BslAdvectionRTheta::operator() (
    DFieldRTheta allfdistribu,
    DConstVectorFieldRTheta < X , Y > advection_field_xy,
    double dt
) override const
```





**Parameters:**


* `allfdistribu` A Field containing the values of the function we want to advect. 
* `advection_field_xy` A DConstVectorFieldRTheta containing the values of the advection field on the physical domain axes. 
* `dt` A time step used.



**Returns:**

A Field to allfdistribu advected on the time step given. 





        
Implements [*IAdvectionRTheta::operator()*](classIAdvectionRTheta.md#function-operator)


<hr>



### function operator() 

_Allocate a Field to the advected function._ 
```C++
inline virtual DFieldRTheta BslAdvectionRTheta::operator() (
    DFieldRTheta allfdistribu,
    DConstVectorFieldRTheta < R , Theta > advection_field_rtheta,
    DVector < X , Y > const & advection_field_xy_centre,
    double dt
) override const
```





**Parameters:**


* `allfdistribu` A Field containing the values of the function we want to advect. 
* `advection_field_rtheta` A DConstVectorFieldRTheta containing the values of the advection field on the logical index range axis. It is expressed on the contravariant basis. 
* `advection_field_xy_centre` A vector in the Cartesian basis, containing the value of the advection field at the O-point. 
* `dt` A time step used.



**Returns:**

A Field to allfdistribu advected on the time step given. 





        
Implements [*IAdvectionRTheta::operator()*](classIAdvectionRTheta.md#function-operator_1)


<hr>



### function ~BslAdvectionRTheta 

```C++
BslAdvectionRTheta::~BslAdvectionRTheta () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/advection/bsl_advection_rtheta.hpp`

