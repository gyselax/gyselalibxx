

# Class PolarSplineEvaluator

**template &lt;class PolarBSplinesType, class OuterExtrapolationRule&gt;**



[**ClassList**](annotated.md) **>** [**PolarSplineEvaluator**](classPolarSplineEvaluator.md)



_Define an evaluator on polar B-splines._ [More...](#detailed-description)

* `#include <polar_spline_evaluator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename PolarBSplinesType::BSplinesR\_tag | [**BSplinesR**](#typedef-bsplinesr)  <br>_Tag the type of first dimension B-splines._  |
| typedef typename PolarBSplinesType::BSplinesTheta\_tag | [**BSplinesTheta**](#typedef-bsplinestheta)  <br>_Tag the type of second dimension B-splines._  |
| typedef typename BSplinesR::continuous\_dimension\_type | [**DimR**](#typedef-dimr)  <br>_Tag the first dimension of the space._  |
| typedef typename BSplinesTheta::continuous\_dimension\_type | [**DimTheta**](#typedef-dimtheta)  <br>_Tag the second dimension of the space._  |
| typedef PolarBSplinesType | [**bsplines\_type**](#typedef-bsplines_type)  <br>_Tag the type of B-splines._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  int constexpr | [**continuity**](#variable-continuity)   = `PolarBSplinesType::continuity`<br>_Tag the order of continuity of the B-splines basis at the O-point._  |














## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarSplineEvaluator**](#function-polarsplineevaluator-14) () = delete<br> |
|   | [**PolarSplineEvaluator**](#function-polarsplineevaluator-24) (OuterExtrapolationRule const & outer\_bc) <br>_Instantiate a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _with boundary evaluation conditions._ |
|   | [**PolarSplineEvaluator**](#function-polarsplineevaluator-34) ([**PolarSplineEvaluator**](classPolarSplineEvaluator.md) const & x) = default<br>_Instantiate a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _from another._ |
|   | [**PolarSplineEvaluator**](#function-polarsplineevaluator-44) ([**PolarSplineEvaluator**](classPolarSplineEvaluator.md) && x) = default<br>_Instantiate a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _from another temporary._ |
|  double | [**deriv\_1\_and\_2**](#function-deriv_1_and_2) (ddc::Coordinate&lt; [**DimR**](classPolarSplineEvaluator.md#typedef-dimr), [**DimTheta**](classPolarSplineEvaluator.md#typedef-dimtheta) &gt; coord\_eval, host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef) const<br>_Get the value of the cross derivative of the spline function._  |
|  double | [**deriv\_dim\_1**](#function-deriv_dim_1-12) (ddc::Coordinate&lt; [**DimR**](classPolarSplineEvaluator.md#typedef-dimr), [**DimTheta**](classPolarSplineEvaluator.md#typedef-dimtheta) &gt; coord\_eval, host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef) const<br>_Get the value of the derivative of the spline function on the first dimension._  |
|  void | [**deriv\_dim\_1**](#function-deriv_dim_1-22) (DField&lt; Domain, Kokkos::HostSpace &gt; const spline\_eval, ConstField&lt; ddc::Coordinate&lt; [**DimR**](classPolarSplineEvaluator.md#typedef-dimr), [**DimTheta**](classPolarSplineEvaluator.md#typedef-dimtheta) &gt;, Domain, Kokkos::HostSpace &gt; const coords\_eval, host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef) const<br>_Get the values of the derivative of the spline function on the first dimension._  |
|  void | [**deriv\_dim\_1\_and\_2**](#function-deriv_dim_1_and_2) (DField&lt; Domain, Kokkos::HostSpace &gt; const spline\_eval, ConstField&lt; ddc::Coordinate&lt; [**DimR**](classPolarSplineEvaluator.md#typedef-dimr), [**DimTheta**](classPolarSplineEvaluator.md#typedef-dimtheta) &gt;, Domain, Kokkos::HostSpace &gt; const coords\_eval, host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef) const<br>_Get the values of the cross derivative of the spline function._  |
|  double | [**deriv\_dim\_2**](#function-deriv_dim_2-12) (ddc::Coordinate&lt; [**DimR**](classPolarSplineEvaluator.md#typedef-dimr), [**DimTheta**](classPolarSplineEvaluator.md#typedef-dimtheta) &gt; coord\_eval, host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef) const<br>_Get the value of the derivative of the spline function on the second dimension._  |
|  void | [**deriv\_dim\_2**](#function-deriv_dim_2-22) (DField&lt; Domain, Kokkos::HostSpace &gt; const spline\_eval, ConstField&lt; ddc::Coordinate&lt; [**DimR**](classPolarSplineEvaluator.md#typedef-dimr), [**DimTheta**](classPolarSplineEvaluator.md#typedef-dimtheta) &gt;, Domain, Kokkos::HostSpace &gt; const coords\_eval, host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef) const<br>_Get the values of the derivative of the spline function on the second dimension._  |
|  double | [**integrate**](#function-integrate) (host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef, Mapping const mapping) const<br>_Get the in integral of a spline function over the domain._  |
|  double | [**operator()**](#function-operator) (ddc::Coordinate&lt; [**DimR**](classPolarSplineEvaluator.md#typedef-dimr), [**DimTheta**](classPolarSplineEvaluator.md#typedef-dimtheta) &gt; coord\_eval, host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef) const<br>_Get the value of the spline function at a given coordinate._  |
|  void | [**operator()**](#function-operator_1) (DField&lt; Domain, Kokkos::HostSpace &gt; const spline\_eval, ConstField&lt; ddc::Coordinate&lt; [**DimR**](classPolarSplineEvaluator.md#typedef-dimr), [**DimTheta**](classPolarSplineEvaluator.md#typedef-dimtheta) &gt;, Domain, Kokkos::HostSpace &gt; const coords\_eval, host\_t&lt; [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType &gt; &gt; const spline\_coef) const<br>_Get the values of the spline function on a domain._  |
|  [**PolarSplineEvaluator**](classPolarSplineEvaluator.md) & | [**operator=**](#function-operator_2) ([**PolarSplineEvaluator**](classPolarSplineEvaluator.md) const & x) = default<br>_Assign a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _from another._ |
|  [**PolarSplineEvaluator**](classPolarSplineEvaluator.md) & | [**operator=**](#function-operator_3) ([**PolarSplineEvaluator**](classPolarSplineEvaluator.md) && x) = default<br>_Assign a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _from another temporary._ |
|   | [**~PolarSplineEvaluator**](#function-polarsplineevaluator) () = default<br> |




























## Detailed Description




**See also:** [**PolarBSplines**](classPolarBSplines.md) 



    
## Public Types Documentation




### typedef BSplinesR 

_Tag the type of first dimension B-splines._ 
```C++
using PolarSplineEvaluator< PolarBSplinesType, OuterExtrapolationRule >::BSplinesR =  typename PolarBSplinesType::BSplinesR_tag;
```




<hr>



### typedef BSplinesTheta 

_Tag the type of second dimension B-splines._ 
```C++
using PolarSplineEvaluator< PolarBSplinesType, OuterExtrapolationRule >::BSplinesTheta =  typename PolarBSplinesType::BSplinesTheta_tag;
```




<hr>



### typedef DimR 

_Tag the first dimension of the space._ 
```C++
using PolarSplineEvaluator< PolarBSplinesType, OuterExtrapolationRule >::DimR =  typename BSplinesR::continuous_dimension_type;
```




<hr>



### typedef DimTheta 

_Tag the second dimension of the space._ 
```C++
using PolarSplineEvaluator< PolarBSplinesType, OuterExtrapolationRule >::DimTheta =  typename BSplinesTheta::continuous_dimension_type;
```




<hr>



### typedef bsplines\_type 

_Tag the type of B-splines._ 
```C++
using PolarSplineEvaluator< PolarBSplinesType, OuterExtrapolationRule >::bsplines_type =  PolarBSplinesType;
```




<hr>
## Public Static Attributes Documentation




### variable continuity 

_Tag the order of continuity of the B-splines basis at the O-point._ 
```C++
int constexpr PolarSplineEvaluator< PolarBSplinesType, OuterExtrapolationRule >::continuity;
```




<hr>
## Public Functions Documentation




### function PolarSplineEvaluator [1/4]

```C++
PolarSplineEvaluator::PolarSplineEvaluator () = delete
```




<hr>



### function PolarSplineEvaluator [2/4]

_Instantiate a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _with boundary evaluation conditions._
```C++
inline explicit PolarSplineEvaluator::PolarSplineEvaluator (
    OuterExtrapolationRule const & outer_bc
) 
```



Instantiate a [**PolarSplineEvaluator**](classPolarSplineEvaluator.md) by specifying how points lying outside the domain should be evaluated. The domain is the domain on which the polar splines are defined.




**Parameters:**


* `outer_bc` A class containing an operator which can be called to provide a boundary value to evaluate a point lying outside the domain. 




        

<hr>



### function PolarSplineEvaluator [3/4]

_Instantiate a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _from another._
```C++
PolarSplineEvaluator::PolarSplineEvaluator (
    PolarSplineEvaluator const & x
) = default
```





**Parameters:**


* `x` Another [**PolarSplineEvaluator**](classPolarSplineEvaluator.md) (lvalue). 




        

<hr>



### function PolarSplineEvaluator [4/4]

_Instantiate a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _from another temporary._
```C++
PolarSplineEvaluator::PolarSplineEvaluator (
    PolarSplineEvaluator && x
) = default
```





**Parameters:**


* `x` Another temporary [**PolarSplineEvaluator**](classPolarSplineEvaluator.md) (rvalue). 




        

<hr>



### function deriv\_1\_and\_2 

_Get the value of the cross derivative of the spline function._ 
```C++
inline double PolarSplineEvaluator::deriv_1_and_2 (
    ddc::Coordinate< DimR , DimTheta > coord_eval,
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef
) const
```





**Parameters:**


* `coord_eval` The coordinate where we want to evaluate. 
* `spline_coef` The B-splines coefficients of the function we want to evaluate.



**Returns:**

The value of the cross derivative of the spline function 





        

<hr>



### function deriv\_dim\_1 [1/2]

_Get the value of the derivative of the spline function on the first dimension._ 
```C++
inline double PolarSplineEvaluator::deriv_dim_1 (
    ddc::Coordinate< DimR , DimTheta > coord_eval,
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef
) const
```





**Parameters:**


* `coord_eval` The coordinate where we want to evaluate. 
* `spline_coef` The B-splines coefficients of the function we want to evaluate.



**Returns:**

The value of the derivative of the spline function on the first dimension. 





        

<hr>



### function deriv\_dim\_1 [2/2]

_Get the values of the derivative of the spline function on the first dimension._ 
```C++
template<class Domain>
inline void PolarSplineEvaluator::deriv_dim_1 (
    DField< Domain, Kokkos::HostSpace > const spline_eval,
    ConstField< ddc::Coordinate< DimR , DimTheta >, Domain, Kokkos::HostSpace > const coords_eval,
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef
) const
```





**Parameters:**


* `spline_eval` The values of the function evaluated on the domain. 
* `coords_eval` The coordinates where we want to evaluate. 
* `spline_coef` The B-splines coefficients of the spline function we want to evaluate. 




        

<hr>



### function deriv\_dim\_1\_and\_2 

_Get the values of the cross derivative of the spline function._ 
```C++
template<class Domain>
inline void PolarSplineEvaluator::deriv_dim_1_and_2 (
    DField< Domain, Kokkos::HostSpace > const spline_eval,
    ConstField< ddc::Coordinate< DimR , DimTheta >, Domain, Kokkos::HostSpace > const coords_eval,
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef
) const
```





**Parameters:**


* `spline_eval` The values of the function evaluated on the domain. 
* `coords_eval` The coordinates where we want to evaluate. 
* `spline_coef` The B-splines coefficients of the splinefunction we want to evaluate. 




        

<hr>



### function deriv\_dim\_2 [1/2]

_Get the value of the derivative of the spline function on the second dimension._ 
```C++
inline double PolarSplineEvaluator::deriv_dim_2 (
    ddc::Coordinate< DimR , DimTheta > coord_eval,
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef
) const
```





**Parameters:**


* `coord_eval` The coordinate where we want to evaluate. 
* `spline_coef` The B-splines coefficients of the function we want to evaluate.



**Returns:**

The value of the derivative of the spline function on the second dimension. 





        

<hr>



### function deriv\_dim\_2 [2/2]

_Get the values of the derivative of the spline function on the second dimension._ 
```C++
template<class Domain>
inline void PolarSplineEvaluator::deriv_dim_2 (
    DField< Domain, Kokkos::HostSpace > const spline_eval,
    ConstField< ddc::Coordinate< DimR , DimTheta >, Domain, Kokkos::HostSpace > const coords_eval,
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef
) const
```





**Parameters:**


* `spline_eval` The values of the function evaluated on the domain. 
* `coords_eval` The coordinates where we want to evaluate. 
* `spline_coef` The B-splines coefficients of the spline function we want to evaluate.. 




        

<hr>



### function integrate 

_Get the in integral of a spline function over the domain._ 
```C++
template<class Mapping>
inline double PolarSplineEvaluator::integrate (
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef,
    Mapping const mapping
) const
```





**Parameters:**


* `spline_coef` The B-splines coefficients of the function we want to evaluate. 
* `mapping` The mapping function.



**Returns:**

The integral of the spline function over the domain. 





        

<hr>



### function operator() 

_Get the value of the spline function at a given coordinate._ 
```C++
inline double PolarSplineEvaluator::operator() (
    ddc::Coordinate< DimR , DimTheta > coord_eval,
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef
) const
```





**Parameters:**


* `coord_eval` The coordinate where we want to evaluate. 
* `spline_coef` The B-splines coefficients of the function we want to evaluate.



**Returns:**

A double with value of the spline function at the given coordinate. 





        

<hr>



### function operator() 

_Get the values of the spline function on a domain._ 
```C++
template<class Domain>
inline void PolarSplineEvaluator::operator() (
    DField< Domain, Kokkos::HostSpace > const spline_eval,
    ConstField< ddc::Coordinate< DimR , DimTheta >, Domain, Kokkos::HostSpace > const coords_eval,
    host_t< ConstPolarSpline < PolarBSplinesType > > const spline_coef
) const
```





**Parameters:**


* `spline_eval` The values of the function evaluated on the domain. 
* `coords_eval` The coordinates where we want to evaluate. 
* `spline_coef` The B-splines coefficients of the spline function we want to evaluate. 




        

<hr>



### function operator= 

_Assign a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _from another._
```C++
PolarSplineEvaluator & PolarSplineEvaluator::operator= (
    PolarSplineEvaluator const & x
) = default
```





**Parameters:**


* `x` Another [**PolarSplineEvaluator**](classPolarSplineEvaluator.md) (lvalue).



**Returns:**

[**PolarSplineEvaluator**](classPolarSplineEvaluator.md) assigned. 





        

<hr>



### function operator= 

_Assign a_ [_**PolarSplineEvaluator**_](classPolarSplineEvaluator.md) _from another temporary._
```C++
PolarSplineEvaluator & PolarSplineEvaluator::operator= (
    PolarSplineEvaluator && x
) = default
```





**Parameters:**


* `x` Another temporary [**PolarSplineEvaluator**](classPolarSplineEvaluator.md) (rvalue).



**Returns:**

[**PolarSplineEvaluator**](classPolarSplineEvaluator.md) assigned. 





        

<hr>



### function ~PolarSplineEvaluator 

```C++
PolarSplineEvaluator::~PolarSplineEvaluator () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_spline_evaluator.hpp`

