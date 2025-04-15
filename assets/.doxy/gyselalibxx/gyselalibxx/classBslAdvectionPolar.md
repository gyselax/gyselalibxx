

# Class BslAdvectionPolar

**template &lt;class FootFinder, class LogicalToPhysicalMapping, class InterpolatorPolar&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvectionPolar**](classBslAdvectionPolar.md)



_Define an advection operator on 2D_  _domain._[More...](#detailed-description)

* `#include <bsl_advection_polar.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvectionPolar**](#function-bsladvectionpolar) (InterpolatorPolar const & function\_interpolator, FootFinder const & foot\_finder, LogicalToPhysicalMapping const & logical\_to\_physical\_mapping) <br>_Instantiate an advection operator._  |
|  DFieldFDistribu | [**operator()**](#function-operator) (DFieldFDistribu allfdistribu, [**DVectorConstFieldAdvectionXY**](classVectorField.md) advection\_field\_xy, double dt) const<br>_Allocate a Field of the advected function._  |
|  DFieldFDistribu | [**operator()**](#function-operator_1) (DFieldFDistribu allfdistribu, [**DVectorConstFieldAdvectionRTheta**](classVectorField.md) advection\_field\_rtheta, [**DTensor**](classTensor.md)&lt; CartesianBasis &gt; const & advection\_field\_xy\_centre, double dt) const<br>_Allocate a Field to the advected function._  |
|   | [**~BslAdvectionPolar**](#function-bsladvectionpolar) () = default<br> |




























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




### function BslAdvectionPolar 

_Instantiate an advection operator._ 
```C++
inline BslAdvectionPolar::BslAdvectionPolar (
    InterpolatorPolar const & function_interpolator,
    FootFinder const & foot_finder,
    LogicalToPhysicalMapping const & logical_to_physical_mapping
) 
```





**Parameters:**


* `function_interpolator` The polar interpolator to interpolate the function once the characteristics have been computed. 
* `foot_finder` An IFootFinder which computes the feet of the characteristics. 
* `logical_to_physical_mapping` The mapping function from the logical domain to the physical domain.



**Template parameters:**


* `IFootFinder` A child class of IFootFinder. 




        

<hr>



### function operator() 

_Allocate a Field of the advected function._ 
```C++
inline DFieldFDistribu BslAdvectionPolar::operator() (
    DFieldFDistribu allfdistribu,
    DVectorConstFieldAdvectionXY advection_field_xy,
    double dt
) const
```





**Parameters:**


* `allfdistribu` A Field containing the values of the function we want to advect. 
* `advection_field_xy` A field of vectors defined on the Cartesian basis containing the values of the advection field at each point on the logical grid. 
* `dt` A time step used.



**Returns:**

A Field to allfdistribu advected on the time step given. 





        

<hr>



### function operator() 

_Allocate a Field to the advected function._ 
```C++
inline DFieldFDistribu BslAdvectionPolar::operator() (
    DFieldFDistribu allfdistribu,
    DVectorConstFieldAdvectionRTheta advection_field_rtheta,
    DTensor < CartesianBasis > const & advection_field_xy_centre,
    double dt
) const
```





**Parameters:**


* `allfdistribu` A Field containing the values of the function we want to advect. 
* `advection_field_rtheta` A field of vectors defined on the Curvilinear basis containing the values of the advection field at each point on the logical grid. It is expressed on the contravariant basis. 
* `advection_field_xy_centre` A vector in the Cartesian basis, containing the value of the advection field at the O-point. 
* `dt` A time step used.



**Returns:**

A Field to allfdistribu advected on the time step given. 





        

<hr>



### function ~BslAdvectionPolar 

```C++
BslAdvectionPolar::~BslAdvectionPolar () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/bsl_advection_polar.hpp`

