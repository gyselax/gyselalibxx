

# Struct PolarSpline

**template &lt;class PolarBSplinesType, class MemSpace&gt;**



[**ClassList**](annotated.md) **>** [**PolarSpline**](structPolarSpline.md)



_A structure containing the two Fields necessary to define a reference to a spline on a set of polar basis splines._ [More...](#detailed-description)

* `#include <polar_spline.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename PolarBSplinesType::BSplinesR\_tag | [**BSplinesR**](#typedef-bsplinesr)  <br>_The radial bspline from which the polar B-splines are constructed._  |
| typedef typename PolarBSplinesType::BSplinesTheta\_tag | [**BSplinesTheta**](#typedef-bsplinestheta)  <br>_The poloidal bspline from which the polar B-splines are constructed._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  DField&lt; IdxRange&lt; PolarBSplinesType &gt;, MemSpace &gt; | [**singular\_spline\_coef**](#variable-singular_spline_coef)  <br> |
|  DField&lt; IdxRange&lt; [**BSplinesR**](structPolarSpline.md#typedef-bsplinesr), [**BSplinesTheta**](structPolarSpline.md#typedef-bsplinestheta) &gt;, MemSpace &gt; | [**spline\_coef**](#variable-spline_coef)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarSpline**](#function-polarspline) ([**PolarSplineMem**](structPolarSplineMem.md)&lt; PolarBSplinesType, MemSpace &gt; & spl) <br> |
|  [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; | [**span\_cview**](#function-span_cview) () const<br> |
|  [**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; | [**span\_view**](#function-span_view) () <br> |




























## Detailed Description




**Template parameters:**


* `PolarBSplinesType` The type of the polar B-splines on which this spline is defined. 




    
## Public Types Documentation




### typedef BSplinesR 

_The radial bspline from which the polar B-splines are constructed._ 
```C++
using PolarSpline< PolarBSplinesType, MemSpace >::BSplinesR =  typename PolarBSplinesType::BSplinesR_tag;
```




<hr>



### typedef BSplinesTheta 

_The poloidal bspline from which the polar B-splines are constructed._ 
```C++
using PolarSpline< PolarBSplinesType, MemSpace >::BSplinesTheta =  typename PolarBSplinesType::BSplinesTheta_tag;
```




<hr>
## Public Attributes Documentation




### variable singular\_spline\_coef 

```C++
DField<IdxRange<PolarBSplinesType>, MemSpace> PolarSpline< PolarBSplinesType, MemSpace >::singular_spline_coef;
```



A Field containing the coefficients in front of the b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines. 


        

<hr>



### variable spline\_coef 

```C++
DField<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> PolarSpline< PolarBSplinesType, MemSpace >::spline_coef;
```



A Field containing the coefficients in front of the b-spline elements which can be expressed as a tensor product of 1D B-splines. 


        

<hr>
## Public Functions Documentation




### function PolarSpline 

```C++
inline explicit PolarSpline::PolarSpline (
    PolarSplineMem < PolarBSplinesType, MemSpace > & spl
) 
```



Construct a reference to a [**PolarSplineMem**](structPolarSplineMem.md).




**Parameters:**


* `spl` The [**PolarSplineMem**](structPolarSplineMem.md) being referenced. 




        

<hr>



### function span\_cview 

```C++
inline ConstPolarSpline < PolarBSplinesType, MemSpace > PolarSpline::span_cview () const
```



Get a constant reference to the polar spline referenced by this polar spline view.




**Returns:**

A constant reference to a polar spline. 





        

<hr>



### function span\_view 

```C++
inline PolarSpline < PolarBSplinesType, MemSpace > PolarSpline::span_view () 
```



Get a modifiable reference to the polar spline referenced by this polar spline view.




**Returns:**

A modifiable reference to a polar spline. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_spline.hpp`

