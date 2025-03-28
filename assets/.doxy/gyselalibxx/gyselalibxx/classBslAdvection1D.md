

# Class BslAdvection1D

**template &lt;class GridInterest, class IdxRangeAdvection, class IdxRangeFunction, class AdvectionFieldBuilder, class AdvectionFieldEvaluator, class TimeStepper&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvection1D**](classBslAdvection1D.md)



_A class which computes the advection along the dimension of interest GridInterest._ [More...](#detailed-description)

* `#include <bsl_advection_1d.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvection1D**](#function-bsladvection1d) (FunctionPreallocatableInterpolatorType const & function\_interpolator, AdvectionFieldBuilder const & adv\_field\_builder, AdvectionFieldEvaluator const & adv\_field\_evaluator, TimeStepper const & time\_stepper) <br>_Constructor when the advection index range and the function index range are different._  |
|  FunctionField | [**operator()**](#function-operator) (FunctionField const allfdistribu, AdvecField const advection\_field, double const dt, std::optional&lt; AdvecFieldDerivConstField &gt; const advection\_field\_derivatives\_min=std::nullopt, std::optional&lt; AdvecFieldDerivConstField &gt; const advection\_field\_derivatives\_max=std::nullopt) const<br>_Advects allfdistribu along the advection dimension GridInterest for a duration dt._  |
|   | [**~BslAdvection1D**](#function-bsladvection1d) () = default<br> |




























## Detailed Description


This operator solves the following equation type





with
* , a function defined on an index range ;
* , an advection field defined on subindex range ;
* , an advection dimension.




The characteristic equation is solved on the advection index range . Then the feet on  are computed from the characteristic feet on  and the function  is interpolated at the feet in .


The characteristic equation is solved using a time integration method ([**ITimeStepper**](classITimeStepper.md)).




**Template parameters:**


* `GridInterest` The dimension along which the advection is computed. It refers to the dimension of  in the equation. 
* `IdxRangeAdvection` The index range  where the characteristic equation is solved. It also refers to the index range of the advection field. It had to also be defined on the GridInterest for the time integration method. 
* `IdxRangeFunction` The index range  where allfdistribu is defined. 
* `AdvectionFieldBuilder` The type of the spline builder for the advection field (see SplineBuilder). 
* `AdvectionFieldEvaluator` The type of the spline evaluator for the advection field (see SplineEvaluator). 
 
* `TimeStepper` The time integration method applied to solve the characteristic equation. The method is picked among the child classes of [**ITimeStepper**](classITimeStepper.md). 




    
## Public Functions Documentation




### function BslAdvection1D 

_Constructor when the advection index range and the function index range are different._ 
```C++
inline explicit BslAdvection1D::BslAdvection1D (
    FunctionPreallocatableInterpolatorType const & function_interpolator,
    AdvectionFieldBuilder const & adv_field_builder,
    AdvectionFieldEvaluator const & adv_field_evaluator,
    TimeStepper const & time_stepper
) 
```



When IdxRangeAdvection and IdxRangeFunction are different, we need one interpolator for each index range.


We can also use it when we want two different interpolators but defined on the same index range (e.g. different boundary conditions for the evaluators).




**Parameters:**


* `function_interpolator` interpolator along the GridInterest direction to interpolate the advected function (allfdistribu) on the index range of the function. 
* `adv_field_builder` builder along the GridInterest direction to build a spline representation of the advection field on the function index range. 
* `adv_field_evaluator` evaluator along the GridInterest direction to evaluate the advection field spline representation on the function index range. 
 
* `time_stepper` time integration method for the characteristic equation. 




        

<hr>



### function operator() 

_Advects allfdistribu along the advection dimension GridInterest for a duration dt._ 
```C++
inline FunctionField BslAdvection1D::operator() (
    FunctionField const allfdistribu,
    AdvecField const advection_field,
    double const dt,
    std::optional< AdvecFieldDerivConstField > const advection_field_derivatives_min=std::nullopt,
    std::optional< AdvecFieldDerivConstField > const advection_field_derivatives_max=std::nullopt
) const
```





**Parameters:**


* `allfdistribu` Reference to the advected function, allocated on the device 
* `advection_field` Reference to the advection field, allocated on the device. 
* `advection_field_derivatives_min` Reference to the advection field derivatives at the left side of the interest dimension, allocated on the device. 
* `advection_field_derivatives_max` Reference to the advection field derivatives at the right side of the interest dimension, allocated on the device. 
* `dt` Time step.



**Returns:**

A reference to the allfdistribu array after advection on dt. 





        

<hr>



### function ~BslAdvection1D 

```C++
BslAdvection1D::~BslAdvection1D () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/bsl_advection_1d.hpp`

