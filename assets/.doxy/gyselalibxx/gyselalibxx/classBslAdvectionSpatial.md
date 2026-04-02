

# Class BslAdvectionSpatial

**template &lt;class Geometry, concepts::Interpolation1D FunctionInterpolator, class DataType&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvectionSpatial**](classBslAdvectionSpatial.md)



_A class which computes the spatial advection along the dimension of interest_ [_**GridX**_](structGridX.md) _. Working for every Cartesian geometry._

* `#include <bsl_advection_x.hpp>`



Inherits the following classes: [IAdvectionSpatial](classIAdvectionSpatial.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvectionSpatial**](#function-bsladvectionspatial-12) (FunctionBuilder const & function\_builder, FunctionEvaluator const & function\_evaluator) <br>_Constructor._  |
|   | [**BslAdvectionSpatial**](#function-bsladvectionspatial-22) (FunctionInterpolator const & function\_interpolator) <br>_Constructor._  |
|  Field&lt; DataType, IdxRangeFdistrib &gt; | [**operator()**](#function-operator) (Field&lt; DataType, IdxRangeFdistrib &gt; const allfdistribu, DataType const dt) override const<br>_Advects fdistribu along_ [_**GridX**_](structGridX.md) _for a duration dt._ |
|   | [**~BslAdvectionSpatial**](#function-bsladvectionspatial) () override<br> |


## Public Functions inherited from IAdvectionSpatial

See [IAdvectionSpatial](classIAdvectionSpatial.md)

| Type | Name |
| ---: | :--- |
| virtual Field&lt; DataType, typename Geometry::IdxRangeFdistribu &gt; | [**operator()**](classIAdvectionSpatial.md#function-operator) (Field&lt; DataType, typename Geometry::IdxRangeFdistribu &gt; allfdistribu, DataType dt) const = 0<br>_operates a transport of the distribution function._  |
| virtual  | [**~IAdvectionSpatial**](classIAdvectionSpatial.md#function-iadvectionspatial) () = default<br> |






















































## Public Functions Documentation




### function BslAdvectionSpatial [1/2]

_Constructor._ 
```C++
inline explicit BslAdvectionSpatial::BslAdvectionSpatial (
    FunctionBuilder const & function_builder,
    FunctionEvaluator const & function_evaluator
) 
```





**Parameters:**


* `function_builder` Builder along the [**GridX**](structGridX.md) direction used to build the interpolation representation of the advected function. 
* `function_evaluator` Evaluator along the [**GridX**](structGridX.md) direction used to evaluate the advected function at the characteristic feet. 




        

<hr>



### function BslAdvectionSpatial [2/2]

_Constructor._ 
```C++
inline explicit BslAdvectionSpatial::BslAdvectionSpatial (
    FunctionInterpolator const & function_interpolator
) 
```





**Parameters:**


* `function_interpolator` Interpolator along the [**GridX**](structGridX.md) direction used to build and evaluate the interpolation representation of the advected function. 




        

<hr>



### function operator() 

_Advects fdistribu along_ [_**GridX**_](structGridX.md) _for a duration dt._
```C++
inline Field< DataType, IdxRangeFdistrib > BslAdvectionSpatial::operator() (
    Field< DataType, IdxRangeFdistrib > const allfdistribu,
    DataType const dt
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

