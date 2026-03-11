

# Class BslAdvection1D

**template &lt;class GridInterest, class IdxRangeAdvection, class IdxRangeFunction, concepts::Interpolation FunctionInterpolator, concepts::Interpolation AdvectionFieldInterpolator, class TimeStepperBuilder, class DataType&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvection1D**](classBslAdvection1D.md)



_A class which computes the advection along the dimension of interest GridInterest._ [More...](#detailed-description)

* `#include <bsl_advection_1d.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename AdvectionFieldInterpolator::BuilderType | [**AdvectionFieldBuilder**](#typedef-advectionfieldbuilder)  <br>_The type of the spline builder for the advection field (see SplineBuilder)._  |
| typedef typename AdvectionFieldInterpolator::EvaluatorType | [**AdvectionFieldEvaluator**](#typedef-advectionfieldevaluator)  <br>_The type of the spline evaluator for the advection field (see SplineEvaluator)._  |
| typedef typename FunctionInterpolator::BuilderType | [**FunctionBuilder**](#typedef-functionbuilder)  <br>_The type of the interpolation builder for the advected function along GridInterest._  |
| typedef typename FunctionInterpolator::EvaluatorType | [**FunctionEvaluator**](#typedef-functionevaluator)  <br>_The type of the interpolation evaluator for the advected function along GridInterest._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvection1D**](#function-bsladvection1d-12) ([**FunctionBuilder**](classBslAdvection1D.md#typedef-functionbuilder) const & function\_builder, [**FunctionEvaluator**](classBslAdvection1D.md#typedef-functionevaluator) const & function\_evaluator, [**AdvectionFieldBuilder**](classBslAdvection1D.md#typedef-advectionfieldbuilder) const & adv\_field\_builder, [**AdvectionFieldEvaluator**](classBslAdvection1D.md#typedef-advectionfieldevaluator) const & adv\_field\_evaluator, TimeStepperBuilder const & time\_stepper\_builder) <br>_Constructor when the advection domain and the function domain are different._  |
|   | [**BslAdvection1D**](#function-bsladvection1d-22) (FunctionInterpolator const & function\_interpolator, AdvectionFieldInterpolator const & adv\_field\_interpolator, TimeStepperBuilder const & time\_stepper\_builder) <br>_Constructor when the advection domain and the function domain are different._  |
|  FunctionField | [**operator()**](#function-operator) (FunctionField const allfdistribu, AdvecField const advection\_field, DataType const dt, std::optional&lt; AdvecFieldDerivConstField &gt; const advection\_field\_derivatives\_min=std::nullopt, std::optional&lt; AdvecFieldDerivConstField &gt; const advection\_field\_derivatives\_max=std::nullopt) const<br>_Advects allfdistribu along the advection dimension GridInterest for a duration dt._  |
|   | [**~BslAdvection1D**](#function-bsladvection1d) () = default<br> |




























## Detailed Description


This operator solves the following equation type



\[\partial_t f_s(t,x) + A_{s, x_i} (x') \cdot \partial_{x_i} f_s (t, x) = 0, \qquad x\in \Omega, x'\in\Omega'\]



with



* \(f\), a function defined on an domain \(\Omega\);
* \(A\), an advection field defined on subdomain \(\Omega'\subset \Omega\);
* \(x_i\), an advection dimension.




The characteristic equation is solved on the advection domain \(\Omega'\). Then the feet on \(\Omega\) are computed from the characteristic feet on \(\Omega'\) and the function \(f\) is interpolated at the feet in \(\Omega\).


The characteristic equation is solved using a time integration method ([**ITimeStepper**](classITimeStepper.md)).




**Template parameters:**


* `GridInterest` The dimension along which the advection is computed. It refers to the dimension of \(x_i\) in the equation. 
* `IdxRangeAdvection` The index range for the interpolation points defined on \(\Omega'\) where the characteristic equation is solved. \(\Omega'\) also refers to the domain of the advection field. It had to also be defined on the GridInterest for the time integration method. 
* `IdxRangeFunction` The index range of \(\Omega\) where allfdistribu is defined. 
* `FunctionInterpolator` The type of the interpolation for the advected function along GridInterest. 
* `AdvectionFieldInterpolator` The type of the interpolation for the advection field (see SplineBuilder). 
* `TimeStepperBuilder` A time stepper builder indicating which time integration method should be applied to solve the characteristic equation. 




    
## Public Types Documentation




### typedef AdvectionFieldBuilder 

_The type of the spline builder for the advection field (see SplineBuilder)._ 
```C++
using BslAdvection1D< GridInterest, IdxRangeAdvection, IdxRangeFunction, FunctionInterpolator, AdvectionFieldInterpolator, TimeStepperBuilder, DataType >::AdvectionFieldBuilder =  typename AdvectionFieldInterpolator::BuilderType;
```




<hr>



### typedef AdvectionFieldEvaluator 

_The type of the spline evaluator for the advection field (see SplineEvaluator)._ 
```C++
using BslAdvection1D< GridInterest, IdxRangeAdvection, IdxRangeFunction, FunctionInterpolator, AdvectionFieldInterpolator, TimeStepperBuilder, DataType >::AdvectionFieldEvaluator =  typename AdvectionFieldInterpolator::EvaluatorType;
```




<hr>



### typedef FunctionBuilder 

_The type of the interpolation builder for the advected function along GridInterest._ 
```C++
using BslAdvection1D< GridInterest, IdxRangeAdvection, IdxRangeFunction, FunctionInterpolator, AdvectionFieldInterpolator, TimeStepperBuilder, DataType >::FunctionBuilder =  typename FunctionInterpolator::BuilderType;
```




<hr>



### typedef FunctionEvaluator 

_The type of the interpolation evaluator for the advected function along GridInterest._ 
```C++
using BslAdvection1D< GridInterest, IdxRangeAdvection, IdxRangeFunction, FunctionInterpolator, AdvectionFieldInterpolator, TimeStepperBuilder, DataType >::FunctionEvaluator =  typename FunctionInterpolator::EvaluatorType;
```




<hr>
## Public Functions Documentation




### function BslAdvection1D [1/2]

_Constructor when the advection domain and the function domain are different._ 
```C++
inline explicit BslAdvection1D::BslAdvection1D (
    FunctionBuilder const & function_builder,
    FunctionEvaluator const & function_evaluator,
    AdvectionFieldBuilder const & adv_field_builder,
    AdvectionFieldEvaluator const & adv_field_evaluator,
    TimeStepperBuilder const & time_stepper_builder
) 
```



When IdxRangeAdvection and IdxRangeFunction are different, we need one builder and evaluator for each index range.


We can also use it when we want two different builders/evaluators but defined on the same domain (e.g. different boundary conditions for the evaluators).




**Parameters:**


* `function_builder` Builder along the GridInterest direction used to build the interpolation representation of the advected function (allfdistribu) on the domain of the function. 
* `function_evaluator` Evaluator along the GridInterest direction used to evaluate the advected function at the characteristic feet. 
* `adv_field_builder` Builder along the GridInterest direction to build a spline representation of the advection field on the function domain. 
* `adv_field_evaluator` Evaluator along the GridInterest direction to evaluate the advection field spline representation on the function domain. 
* `time_stepper_builder` A builder for the time integration method used for the characteristic equation. 




        

<hr>



### function BslAdvection1D [2/2]

_Constructor when the advection domain and the function domain are different._ 
```C++
inline explicit BslAdvection1D::BslAdvection1D (
    FunctionInterpolator const & function_interpolator,
    AdvectionFieldInterpolator const & adv_field_interpolator,
    TimeStepperBuilder const & time_stepper_builder
) 
```



When IdxRangeAdvection and IdxRangeFunction are different, we need one builder and evaluator for each index range.


We can also use it when we want two different builders/evaluators but defined on the same domain (e.g. different boundary conditions for the evaluators).




**Parameters:**


* `function_interpolator` Interpolator along the GridInterest direction used to build a continuous representation and evaluate the advected function at the characteristic feet. 
* `adv_field_interpolator` Interpolator along the GridInterest direction used to build a continuous representation and evaluate the advection field at the characteristic feet. 
* `time_stepper_builder` A builder for the time integration method used for the characteristic equation. 




        

<hr>



### function operator() 

_Advects allfdistribu along the advection dimension GridInterest for a duration dt._ 
```C++
inline FunctionField BslAdvection1D::operator() (
    FunctionField const allfdistribu,
    AdvecField const advection_field,
    DataType const dt,
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

