

# Struct PolarSplineMem

**template &lt;class PolarBSplinesType, class MemSpace&gt;**



[**ClassList**](annotated.md) **>** [**PolarSplineMem**](structPolarSplineMem.md)



_A structure containing the two FieldMems necessary to define a spline on a set of polar basis splines._ [More...](#detailed-description)

* `#include <polar_spline.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename PolarBSplinesType::BSplinesR\_tag | [**BSplinesR**](#typedef-bsplinesr)  <br>_The radial bspline from which the polar B-splines are constructed._  |
| typedef typename PolarBSplinesType::BSplinesTheta\_tag | [**BSplinesTheta**](#typedef-bsplinestheta)  <br>_The poloidal bspline from which the polar B-splines are constructed._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  DFieldMem&lt; IdxRange&lt; PolarBSplinesType &gt;, MemSpace &gt; | [**singular\_spline\_coef**](#variable-singular_spline_coef)  <br> |
|  DFieldMem&lt; IdxRange&lt; [**BSplinesR**](structPolarSplineMem.md#typedef-bsplinesr), [**BSplinesTheta**](structPolarSplineMem.md#typedef-bsplinestheta) &gt;, MemSpace &gt; | [**spline\_coef**](#variable-spline_coef)  <br> |
















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PolarSplineMem**](#function-polarsplinemem-12) (IdxRange&lt; [**BSplinesR**](structPolarSplineMem.md#typedef-bsplinesr), [**BSplinesTheta**](structPolarSplineMem.md#typedef-bsplinestheta) &gt; domain) <br>_A constructor for the_ [_**PolarSplineMem**_](structPolarSplineMem.md) _._ |
|   | [**PolarSplineMem**](#function-polarsplinemem-22) (IdxRange&lt; PolarBSplinesType &gt; singular\_domain, IdxRange&lt; [**BSplinesR**](structPolarSplineMem.md#typedef-bsplinesr), [**BSplinesTheta**](structPolarSplineMem.md#typedef-bsplinestheta) &gt; domain) <br>_Initialise the FieldMems containing the coefficients of the spline representation._  |
|  [**ConstPolarSpline**](structConstPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; | [**span\_cview**](#function-span_cview) () const<br> |
|  [**PolarSpline**](structPolarSpline.md)&lt; PolarBSplinesType, MemSpace &gt; | [**span\_view**](#function-span_view) () <br> |




























## Detailed Description




**Template parameters:**


* `PolarBSplinesType` The type of the polar B-splines on which this spline is defined. 




    
## Public Types Documentation




### typedef BSplinesR 

_The radial bspline from which the polar B-splines are constructed._ 
```C++
using PolarSplineMem< PolarBSplinesType, MemSpace >::BSplinesR =  typename PolarBSplinesType::BSplinesR_tag;
```




<hr>



### typedef BSplinesTheta 

_The poloidal bspline from which the polar B-splines are constructed._ 
```C++
using PolarSplineMem< PolarBSplinesType, MemSpace >::BSplinesTheta =  typename PolarBSplinesType::BSplinesTheta_tag;
```




<hr>
## Public Attributes Documentation




### variable singular\_spline\_coef 

```C++
DFieldMem<IdxRange<PolarBSplinesType>, MemSpace> PolarSplineMem< PolarBSplinesType, MemSpace >::singular_spline_coef;
```



A FieldMem containing the coefficients in front of the b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines. 


        

<hr>



### variable spline\_coef 

```C++
DFieldMem<IdxRange<BSplinesR, BSplinesTheta>, MemSpace> PolarSplineMem< PolarBSplinesType, MemSpace >::spline_coef;
```



A FieldMem containing the coefficients in front of the b-spline elements which can be expressed as a tensor product of 1D B-splines. 


        

<hr>
## Public Functions Documentation




### function PolarSplineMem [1/2]

_A constructor for the_ [_**PolarSplineMem**_](structPolarSplineMem.md) _._
```C++
inline explicit PolarSplineMem::PolarSplineMem (
    IdxRange< BSplinesR , BSplinesTheta > domain
) 
```



A constructor for the [**PolarSplineMem**](structPolarSplineMem.md) which takes a 2D domain of B-splines. This domain is then split into the domains which are relevant for a polar bspline constructed from these 2D B-splines. These new domains are then used to create the chunks which define the spline.




**Parameters:**


* `domain` A 2D domain of B-splines. 




        

<hr>



### function PolarSplineMem [2/2]

_Initialise the FieldMems containing the coefficients of the spline representation._ 
```C++
inline PolarSplineMem::PolarSplineMem (
    IdxRange< PolarBSplinesType > singular_domain,
    IdxRange< BSplinesR , BSplinesTheta > domain
) 
```





**Parameters:**


* `singular_domain` The domain for the coefficients in front of the b-spline elements near the singular point. These are the elements which cannot be expressed as a tensor product of 1D B-splines. 
* `domain` The domain for the coefficients in front of the b-spline elements which can be expressed as a tensor product of 1D B-splines. 




        

<hr>



### function span\_cview 

```C++
inline ConstPolarSpline < PolarBSplinesType, MemSpace > PolarSplineMem::span_cview () const
```



Get a constant reference to this polar spline view.




**Returns:**

A constant reference to this polar spline. 





        

<hr>



### function span\_view 

```C++
inline PolarSpline < PolarBSplinesType, MemSpace > PolarSplineMem::span_view () 
```



Get a modifiable reference to this polar spline.




**Returns:**

A modifiable reference to this polar spline. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_spline.hpp`

