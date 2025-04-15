

# Class PolarBSplines

**template &lt;class [**BSplinesR**](structBSplinesR.md), class [**BSplinesTheta**](structBSplinesTheta.md), int C&gt;**



[**ClassList**](annotated.md) **>** [**PolarBSplines**](classPolarBSplines.md)



[More...](#detailed-description)

* `#include <polar_bsplines.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**Impl**](classPolarBSplines_1_1Impl.md) &lt;class DDim, class MemorySpace&gt;<br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**BSplinesR**](structBSplinesR.md) | [**BSplinesR\_tag**](#typedef-bsplinesr_tag)  <br>_The radial bspline from which the polar B-splines are constructed._  |
| typedef [**BSplinesTheta**](structBSplinesTheta.md) | [**BSplinesTheta\_tag**](#typedef-bsplinestheta_tag)  <br>_The poloidal bspline from which the polar B-splines are constructed._  |
| typedef typename BSplinesR::continuous\_dimension\_type | [**DimR**](#typedef-dimr)  <br>_The tag for the radial direction of the B-splines._  |
| typedef typename BSplinesTheta::continuous\_dimension\_type | [**DimTheta**](#typedef-dimtheta)  <br>_The tag for the poloidal direction of the B-splines._  |
| typedef [**PolarBSplines**](classPolarBSplines.md) | [**discrete\_dimension\_type**](#typedef-discrete_dimension_type)  <br> |
| typedef IdxRange&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**tensor\_product\_idx\_range\_type**](#typedef-tensor_product_idx_range_type)  <br> |
| typedef IdxStep&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**tensor\_product\_idx\_step\_type**](#typedef-tensor_product_idx_step_type)  <br> |
| typedef Idx&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**tensor\_product\_index\_type**](#typedef-tensor_product_index_type)  <br> |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  int constexpr | [**continuity**](#variable-continuity)   = `C`<br>_The continuity enforced by the B-splines at the singular point._  |
















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION [**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type) | [**get\_2d\_index**](#function-get_2d_index) (Idx&lt; DDim &gt; const & idx) <br> |
|  KOKKOS\_FUNCTION Idx&lt; DDim &gt; | [**get\_polar\_index**](#function-get_polar_index) ([**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type) const & idx) <br> |
|  constexpr std::size\_t | [**n\_singular\_basis**](#function-n_singular_basis) () <br> |
|  constexpr IdxRange&lt; DDim &gt; | [**singular\_idx\_range**](#function-singular_idx_range) () <br>_Get the IdxRange containing the indices of the b-splines which traverse the singular point._  |


























## Detailed Description


A class containing all information describing polar B-splines.


Polar B-splines are 2D B-splines with a special treatment for the central singular point of a polar domain. At this singular point new B-splines are created which traverse the singular point and ensure the desired continuity condition.




**Template parameters:**


* [**BSplinesR**](structBSplinesR.md) The basis of radial B-splines from which the polar B-splines are constructed. 
* [**BSplinesTheta**](structBSplinesTheta.md) The poloidal bspline from which the polar B-splines are constructed. 
* `C` The continuity condition. The resulting splines will be continuously differentiable C times. If C == -1 then the resulting spline representation will be discontinuous at the singular point. 




    
## Public Types Documentation




### typedef BSplinesR\_tag 

_The radial bspline from which the polar B-splines are constructed._ 
```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::BSplinesR_tag =  BSplinesR;
```




<hr>



### typedef BSplinesTheta\_tag 

_The poloidal bspline from which the polar B-splines are constructed._ 
```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::BSplinesTheta_tag =  BSplinesTheta;
```




<hr>



### typedef DimR 

_The tag for the radial direction of the B-splines._ 
```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::DimR =  typename BSplinesR::continuous_dimension_type;
```




<hr>



### typedef DimTheta 

_The tag for the poloidal direction of the B-splines._ 
```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::DimTheta =  typename BSplinesTheta::continuous_dimension_type;
```




<hr>



### typedef discrete\_dimension\_type 

```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::discrete_dimension_type =  PolarBSplines;
```



The tag denoting the discrete dimension described by this class.


This is the tag which should be used to create a Field whose contents are each associated with a PolarBSpline. In other words a spline defined on this basis would have the type: DField&lt;IdxRange&lt;PolarBSplines&gt;; 


        

<hr>



### typedef tensor\_product\_idx\_range\_type 

```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::tensor_product_idx_range_type =  IdxRange<BSplinesR, BSplinesTheta>;
```



The type of the 2D idx\_range for the subset of the polar B-splines which can be expressed as a tensor product of 1D B-splines. 


        

<hr>



### typedef tensor\_product\_idx\_step\_type 

```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::tensor_product_idx_step_type =  IdxStep<BSplinesR, BSplinesTheta>;
```



The type of a 2D vector for the subset of the polar B-splines which can be expressed as a tensor product of 1D B-splines. 


        

<hr>



### typedef tensor\_product\_index\_type 

```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::tensor_product_index_type =  Idx<BSplinesR, BSplinesTheta>;
```



The type of a 2D index for the subset of the polar B-splines which can be expressed as a tensor product of 1D B-splines. 


        

<hr>
## Public Static Attributes Documentation




### variable continuity 

_The continuity enforced by the B-splines at the singular point._ 
```C++
int constexpr PolarBSplines< BSplinesR, BSplinesTheta, C >::continuity;
```




<hr>
## Public Static Functions Documentation




### function get\_2d\_index 

```C++
template<class DDim>
static inline KOKKOS_FUNCTION tensor_product_index_type PolarBSplines::get_2d_index (
    Idx< DDim > const & idx
) 
```



Get the 2D index of the tensor product bspline which, when evaluated at the same point, returns the same values as the polar bspline indicated by the index passed as an argument.




**Parameters:**


* `idx` The index of the basis spline in the PolarBSpline index range.



**Returns:**

The index of the equivalent 2D BSpline expressed as a 2D tensor product of 1D BSplines. 





        

<hr>



### function get\_polar\_index 

```C++
template<class DDim>
static inline KOKKOS_FUNCTION Idx< DDim > PolarBSplines::get_polar_index (
    tensor_product_index_type const & idx
) 
```



Get the index of the polar bspline which, when evaluated at the same point, returns the same values as the 2D tensor product bspline indicated by the index passed as an argument.




**Parameters:**


* `idx` The index of a 2D BSpline which is expressed as a tensor product of 1D BSplines.



**Returns:**

The index of the basis spline in the PolarBSpline index range. 





        

<hr>



### function n\_singular\_basis 

```C++
static inline constexpr std::size_t PolarBSplines::n_singular_basis () 
```



Get the number of singular B-splines i.e. B-splines which traverse the singular point.




**Returns:**

The number of B-splines which traverse the singular point. 





        

<hr>



### function singular\_idx\_range 

_Get the IdxRange containing the indices of the b-splines which traverse the singular point._ 
```C++
template<class DDim>
static inline constexpr IdxRange< DDim > PolarBSplines::singular_idx_range () 
```





**Returns:**

The IdxRange containing the indices of the b-splines which traverse the singular point. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_bsplines.hpp`

