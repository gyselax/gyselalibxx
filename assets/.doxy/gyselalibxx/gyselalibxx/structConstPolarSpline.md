

# Struct ConstPolarSpline

**template &lt;class PolarBSplinesType, class MemSpace&gt;**



[**ClassList**](annotated.md) **>** [**ConstPolarSpline**](structConstPolarSpline.md)



_A structure containing the two ConstFields necessary to define a constant reference to a spline on a set of polar basis splines._ [More...](#detailed-description)

* `#include <polar_spline.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename PolarBSplinesType::BSplinesR\_tag | [**BSplinesR**](#typedef-bsplinesr)  <br>_The radial bspline from which the polar B-splines are constructed._  |
| typedef typename PolarBSplinesType::BSplinesTheta\_tag | [**BSplinesTheta**](#typedef-bsplinestheta)  <br>_The poloidal bspline from which the polar B-splines are constructed._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  DConstField&lt; IdxRange&lt; PolarBSplinesType &gt;, MemSpace &gt; const | [**singular\_spline\_coef**](#variable-singular_spline_coef)  <br> |
|  DConstField&lt; IdxRange&lt; [**BSplinesR**](structConstPolarSpline.md#typedef-bsplinesr), [**BSplinesTheta**](structConstPolarSpline.md#typedef-bsplinestheta) &gt;, MemSpace &gt; const | [**spline\_coef**](#variable-spline_coef)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ConstPolarSpline**](#function-constpolarspline-12) ([**PolarSplineMem**](structPolarSplineMem.md)&lt; PolarBSplinesType, MemSpace &gt; const & spl) <br> |
|   | [**ConstPolarSpline**](#function-constpolarspline-22) ([**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; const & spl) <br> |
|  [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; | [**span\_cview**](#function-span_cview) () const<br> |
|  [**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; | [**span\_view**](#function-span_view) () const<br> |




























## Detailed Description




**Template parameters:**


* `PolarBSplinesType` The type of the polar B-splines on which this spline is defined. 




    
## Public Types Documentation




### typedef BSplinesR 

_The radial bspline from which the polar B-splines are constructed._ 
```C++
using ConstPolarSpline< PolarBSplinesType, MemSpace >::BSplinesR =  typename PolarBSplinesType::BSplinesR_tag;
```




<hr>



### typedef BSplinesTheta 

_The poloidal bspline from which the polar B-splines are constructed._ 
```C++
using ConstPolarSpline< PolarBSplinesType, MemSpace >::BSplinesTheta =  typename PolarBSplinesType::BSplinesTheta_tag;
```




<hr>
## Public Attributes Documentation




### variable singular\_spline\_coef 

```C++
DConstField<IdxRange<PolarBSplinesType>, MemSpace> const ConstPolarSpline< PolarBSplinesType, MemSpace >::singular_spline_coef;
```



A ConstField containing the coefficients in front of the b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines. 


        

<hr>



### variable spline\_coef 

```C++
DConstField<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> const ConstPolarSpline< PolarBSplinesType, MemSpace >::spline_coef;
```



A ConstField containing the coefficients in front of the b-spline elements which can be expressed as a tensor product of 1D B-splines. 


        

<hr>
## Public Functions Documentation




### function ConstPolarSpline [1/2]

```C++
inline explicit ConstPolarSpline::ConstPolarSpline (
    PolarSplineMem < PolarBSplinesType, MemSpace > const & spl
) 
```



Construct a constant reference to a [**PolarSplineMem**](structPolarSplineMem.md).




**Parameters:**


* `spl` The [**PolarSplineMem**](structPolarSplineMem.md) being referenced. 




        

<hr>



### function ConstPolarSpline [2/2]

```C++
inline explicit ConstPolarSpline::ConstPolarSpline (
    PolarSpline < PolarBSplinesType, MemSpace > const & spl
) 
```



Construct a constant reference to a [**PolarSplineMem**](structPolarSplineMem.md) from a [**PolarSpline**](structPolarSpline.md)




**Parameters:**


* `spl` The [**PolarSplineMem**](structPolarSplineMem.md) being referenced. 




        

<hr>



### function span\_cview 

```C++
inline ConstPolarSpline < PolarBSplinesType, MemSpace > ConstPolarSpline::span_cview () const
```



Get a constant reference to the polar spline referenced by this polar spline view.




**Returns:**

A constant reference to a polar spline. 





        

<hr>



### function span\_view 

```C++
inline PolarSpline < PolarBSplinesType, MemSpace > ConstPolarSpline::span_view () const
```



Get a reference to the polar spline referenced by this polar spline view.




**Returns:**

A reference to a polar spline. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_spline.hpp`

