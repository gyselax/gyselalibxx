

# Class BslAdvectionVelocity

**template &lt;class Geometry, concepts::Interpolation1D FunctionInterpolator, class DataType&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvectionVelocity**](classBslAdvectionVelocity.md)



_A class which computes the velocity advection along the dimension of interest GridV. Working for every Cartesian geometry._ 

* `#include <bsl_advection_vx.hpp>`



Inherits the following classes: [IAdvectionVelocity](classIAdvectionVelocity.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvectionVelocity**](#function-bsladvectionvelocity-12) (FunctionBuilder const & function\_builder, FunctionEvaluator const & function\_evaluator) <br>_Constructor._  |
|   | [**BslAdvectionVelocity**](#function-bsladvectionvelocity-22) (FunctionInterpolator const & function\_interpolator) <br>_Constructor._  |
|  Field&lt; DataType, IdxRangeFdistribu &gt; | [**operator()**](#function-operator) (Field&lt; DataType, IdxRangeFdistribu &gt; const allfdistribu, ConstField&lt; DataType, IdxRangeSpatial &gt; const electric\_field, DataType const dt) override const<br>_Advects fdistribu along GridV for a duration dt._  |
|   | [**~BslAdvectionVelocity**](#function-bsladvectionvelocity) () override<br> |


## Public Functions inherited from IAdvectionVelocity

See [IAdvectionVelocity](classIAdvectionVelocity.md)

| Type | Name |
| ---: | :--- |
| virtual Field&lt; DataType, typename Geometry::IdxRangeFdistribu &gt; | [**operator()**](classIAdvectionVelocity.md#function-operator) (Field&lt; DataType, typename Geometry::IdxRangeFdistribu &gt; allfdistribu, ConstField&lt; DataType, typename Geometry::IdxRangeSpatial &gt; electric\_field, DataType dt) const = 0<br>_operates a transport of the distribution function._  |
| virtual  | [**~IAdvectionVelocity**](classIAdvectionVelocity.md#function-iadvectionvelocity) () = default<br> |






















































## Public Functions Documentation




### function BslAdvectionVelocity [1/2]

_Constructor._ 
```C++
inline explicit BslAdvectionVelocity::BslAdvectionVelocity (
    FunctionBuilder const & function_builder,
    FunctionEvaluator const & function_evaluator
) 
```





**Parameters:**


* `function_builder` Builder along the GridV direction used to build the interpolation representation of the advected function. 
* `function_evaluator` Evaluator along the GridV direction used to evaluate the advected function at the characteristic feet. 




        

<hr>



### function BslAdvectionVelocity [2/2]

_Constructor._ 
```C++
inline explicit BslAdvectionVelocity::BslAdvectionVelocity (
    FunctionInterpolator const & function_interpolator
) 
```





**Parameters:**


* `function_interpolator` Interpolator along the [**GridVx**](structGridVx.md) direction used to build and evaluate the interpolation representation of the advected function. 




        

<hr>



### function operator() 

_Advects fdistribu along GridV for a duration dt._ 
```C++
inline Field< DataType, IdxRangeFdistribu > BslAdvectionVelocity::operator() (
    Field< DataType, IdxRangeFdistribu > const allfdistribu,
    ConstField< DataType, IdxRangeSpatial > const electric_field,
    DataType const dt
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

