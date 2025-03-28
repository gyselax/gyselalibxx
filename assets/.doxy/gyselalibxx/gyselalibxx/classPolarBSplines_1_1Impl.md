

# Class PolarBSplines::Impl

**template &lt;class DDim, class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**PolarBSplines**](classPolarBSplines.md) **>** [**Impl**](classPolarBSplines_1_1Impl.md)



[More...](#detailed-description)

* `#include <polar_bsplines.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Corner1Tag**](structPolarBSplines_1_1Impl_1_1Corner1Tag.md) <br>_The tag for the first corner of the Barycentric coordinates._  |
| struct | [**Corner2Tag**](structPolarBSplines_1_1Impl_1_1Corner2Tag.md) <br>_The tag for the second corner of the Barycentric coordinates._  |
| struct | [**Corner3Tag**](structPolarBSplines_1_1Impl_1_1Corner3Tag.md) <br>_The tag for the third corner of the Barycentric coordinates._  |
| struct | [**IntermediateBernsteinBasis**](structPolarBSplines_1_1Impl_1_1IntermediateBernsteinBasis.md) &lt;class DiscreteMapping&gt;<br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**PolarBSplines**](classPolarBSplines.md) | [**discrete\_dimension\_type**](#typedef-discrete_dimension_type)  <br>_The tag which should be used to create a Field whose contents are each associated with a PolarBSpline._  |
| typedef IdxRange&lt; DDim &gt; | [**discrete\_domain\_type**](#typedef-discrete_domain_type)  <br>_The type of a index range of_ [_**PolarBSplines**_](classPolarBSplines.md) _._ |
| typedef Idx&lt; DDim &gt; | [**discrete\_element\_type**](#typedef-discrete_element_type)  <br>_The type of an index associated with a PolarBSpline._  |
| typedef IdxStep&lt; DDim &gt; | [**discrete\_vector\_type**](#typedef-discrete_vector_type)  <br>_The type of a vector associated with a PolarBSpline._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Impl**](#function-impl-25) (const DiscreteMapping & curvilinear\_to\_cartesian) <br> |
|   | [**Impl**](#function-impl-35) ([**Impl**](classPolarBSplines_1_1Impl.md)&lt; DDim, OriginMemorySpace &gt; const & impl) <br> |
|   | [**Impl**](#function-impl-45) ([**Impl**](classPolarBSplines_1_1Impl.md) const & x) = default<br> |
|   | [**Impl**](#function-impl-55) ([**Impl**](classPolarBSplines_1_1Impl.md) && x) = default<br> |
|  KOKKOS\_FUNCTION [**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type) | [**eval\_basis**](#function-eval_basis) (DSpan1D singular\_values, DSpan2D values, ddc::Coordinate&lt; [**DimR**](classPolarBSplines.md#typedef-dimr), [**DimTheta**](classPolarBSplines.md#typedef-dimtheta) &gt; p) const<br>_Evaluate the polar basis splines at the coordinate p._  |
|  KOKKOS\_FUNCTION [**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type) | [**eval\_deriv\_r**](#function-eval_deriv_r) (DSpan1D singular\_derivs, DSpan2D derivs, ddc::Coordinate&lt; [**DimR**](classPolarBSplines.md#typedef-dimr), [**DimTheta**](classPolarBSplines.md#typedef-dimtheta) &gt; p) const<br>_Evaluate the radial derivative of the polar basis splines at the coordinate p._  |
|  KOKKOS\_FUNCTION [**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type) | [**eval\_deriv\_r\_and\_theta**](#function-eval_deriv_r_and_theta) (DSpan1D singular\_derivs, DSpan2D derivs, ddc::Coordinate&lt; [**DimR**](classPolarBSplines.md#typedef-dimr), [**DimTheta**](classPolarBSplines.md#typedef-dimtheta) &gt; p) const<br>_Evaluate the second order derivative of the polar basis splines in the radial and poloidal directions, at the coordinate p._  |
|  KOKKOS\_FUNCTION [**tensor\_product\_index\_type**](classPolarBSplines.md#typedef-tensor_product_index_type) | [**eval\_deriv\_theta**](#function-eval_deriv_theta) (DSpan1D singular\_derivs, DSpan2D derivs, ddc::Coordinate&lt; [**DimR**](classPolarBSplines.md#typedef-dimr), [**DimTheta**](classPolarBSplines.md#typedef-dimtheta) &gt; p) const<br>_Evaluate the poloidal derivative of the polar basis splines at the coordinate p._  |
|  [**discrete\_domain\_type**](classPolarBSplines_1_1Impl.md#typedef-discrete_domain_type) | [**full\_domain**](#function-full_domain) () noexcept const<br> |
|  void | [**integrals**](#function-integrals-22) ([**PolarSpline**](structPolarSpline.md)&lt; DDim, MemorySpace2 &gt; int\_vals) const<br> |
|  std::size\_t | [**nbasis**](#function-nbasis) () noexcept const<br> |
|  [**Impl**](classPolarBSplines_1_1Impl.md) & | [**operator=**](#function-operator) ([**Impl**](classPolarBSplines_1_1Impl.md) const & x) = default<br> |
|  [**Impl**](classPolarBSplines_1_1Impl.md) & | [**operator=**](#function-operator_1) ([**Impl**](classPolarBSplines_1_1Impl.md) && x) = default<br> |
|  [**discrete\_domain\_type**](classPolarBSplines_1_1Impl.md#typedef-discrete_domain_type) | [**tensor\_bspline\_idx\_range**](#function-tensor_bspline_idx_range) () noexcept const<br>_Returns the IdxRange containing the indices of the b-splines which don't traverse the singular point and can be expressed as a tensor-product of 1D b-splines._  |
|   | [**~Impl**](#function-impl) () = default<br> |




























## Detailed Description


The [**Impl**](classPolarBSplines_1_1Impl.md) class holds the implementation of the [**PolarBSplines**](classPolarBSplines.md). The implementation is specific to the memory space so that the Fields can be defined with index ranges related to instances of this class.




**Template parameters:**


* `MemorySpace` Indicates where the object is saved. This is either on the host or the device. 




    
## Public Types Documentation




### typedef discrete\_dimension\_type 

_The tag which should be used to create a Field whose contents are each associated with a PolarBSpline._ 
```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::Impl< DDim, MemorySpace >::discrete_dimension_type =  PolarBSplines;
```




<hr>



### typedef discrete\_domain\_type 

_The type of a index range of_ [_**PolarBSplines**_](classPolarBSplines.md) _._
```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::Impl< DDim, MemorySpace >::discrete_domain_type =  IdxRange<DDim>;
```




<hr>



### typedef discrete\_element\_type 

_The type of an index associated with a PolarBSpline._ 
```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::Impl< DDim, MemorySpace >::discrete_element_type =  Idx<DDim>;
```




<hr>



### typedef discrete\_vector\_type 

_The type of a vector associated with a PolarBSpline._ 
```C++
using PolarBSplines< BSplinesR, BSplinesTheta, C >::Impl< DDim, MemorySpace >::discrete_vector_type =  IdxStep<DDim>;
```




<hr>
## Public Functions Documentation




### function Impl [2/5]

```C++
template<class DiscreteMapping>
inline explicit PolarBSplines::Impl::Impl (
    const DiscreteMapping & curvilinear_to_cartesian
) 
```



A constructor for the [**PolarBSplines**](classPolarBSplines.md).




**Parameters:**


* `curvilinear_to_cartesian` A mapping from curvilinear to Cartesian coordinates. This is used to find the singular point and determine the Barycentric coordinates which are used to define the new basis splines which cross the singular point. 




        

<hr>



### function Impl [3/5]

```C++
template<class OriginMemorySpace>
inline explicit PolarBSplines::Impl::Impl (
    Impl < DDim, OriginMemorySpace > const & impl
) 
```



A copy constructor for the [**PolarBSplines**](classPolarBSplines.md).




**Parameters:**


* `impl` The [**PolarBSplines**](classPolarBSplines.md) being copied. 




        

<hr>



### function Impl [4/5]

```C++
PolarBSplines::Impl::Impl (
    Impl const & x
) = default
```



A copy constructor for the [**PolarBSplines**](classPolarBSplines.md).




**Parameters:**


* `x` The [**PolarBSplines**](classPolarBSplines.md) being copied. 




        

<hr>



### function Impl [5/5]

```C++
PolarBSplines::Impl::Impl (
    Impl && x
) = default
```



A copy constructor for the [**PolarBSplines**](classPolarBSplines.md) taking a temporary r-value.




**Parameters:**


* `x` The [**PolarBSplines**](classPolarBSplines.md) being copied. 




        

<hr>



### function eval\_basis 

_Evaluate the polar basis splines at the coordinate p._ 
```C++
KOKKOS_FUNCTION tensor_product_index_type PolarBSplines::Impl::eval_basis (
    DSpan1D singular_values,
    DSpan2D values,
    ddc::Coordinate< DimR , DimTheta > p
) const
```



Evaluate all the b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines, as well as the non-zero b-spline elements which can be expressed as a tensor product of 1D B-splines.




**Parameters:**


* `singular_values` The value of the b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines, evaluated at the coordinate p. 
* `values` The value of the non-zero b-spline elements which can be expressed as a tensor product of 1D B-splines. 
* `p` The coordinate where the basis functions are evaluated.



**Returns:**

The 2D tensor product index of the first b-spline element in the values array. 





        

<hr>



### function eval\_deriv\_r 

_Evaluate the radial derivative of the polar basis splines at the coordinate p._ 
```C++
KOKKOS_FUNCTION tensor_product_index_type PolarBSplines::Impl::eval_deriv_r (
    DSpan1D singular_derivs,
    DSpan2D derivs,
    ddc::Coordinate< DimR , DimTheta > p
) const
```



Evaluate the radial derivative of all the b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines, as well as the non-zero b-spline elements which can be expressed as a tensor product of 1D B-splines.




**Parameters:**


* `singular_derivs` The value of the radial derivative b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines, evaluated at the coordinate p. 
* `derivs` The value of the radial derivative of the non-zero b-spline elements which can be expressed as a tensor product of 1D B-splines. 
* `p` The coordinate where the basis functions are evaluated.



**Returns:**

The 2D tensor product index of the first b-spline element in the values array. 





        

<hr>



### function eval\_deriv\_r\_and\_theta 

_Evaluate the second order derivative of the polar basis splines in the radial and poloidal directions, at the coordinate p._ 
```C++
KOKKOS_FUNCTION tensor_product_index_type PolarBSplines::Impl::eval_deriv_r_and_theta (
    DSpan1D singular_derivs,
    DSpan2D derivs,
    ddc::Coordinate< DimR , DimTheta > p
) const
```



Evaluate the 2nd order derivative of all the b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines, as well as the non-zero b-spline elements which can be expressed as a tensor product of 1D B-splines.




**Parameters:**


* `singular_derivs` The value of the 2nd order derivative b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines, evaluated at the coordinate p. 
* `derivs` The value of the 2nd order derivative of the non-zero b-spline elements which can be expressed as a tensor product of 1D B-splines. 
* `p` The coordinate where the basis functions are evaluated.



**Returns:**

The 2D tensor product index of the first b-spline element in the values array. 





        

<hr>



### function eval\_deriv\_theta 

_Evaluate the poloidal derivative of the polar basis splines at the coordinate p._ 
```C++
KOKKOS_FUNCTION tensor_product_index_type PolarBSplines::Impl::eval_deriv_theta (
    DSpan1D singular_derivs,
    DSpan2D derivs,
    ddc::Coordinate< DimR , DimTheta > p
) const
```



Evaluate the poloidal derivative of all the b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines, as well as the non-zero b-spline elements which can be expressed as a tensor product of 1D B-splines.




**Parameters:**


* `singular_derivs` The value of the poloidal derivative b-spline elements near the singular point which cannot be expressed as a tensor product of 1D B-splines, evaluated at the coordinate p. 
* `derivs` The value of the poloidal derivative of the non-zero b-spline elements which can be expressed as a tensor product of 1D B-splines. 
* `p` The coordinate where the basis functions are evaluated.



**Returns:**

The 2D tensor product index of the first b-spline element in the values array. 





        

<hr>



### function full\_domain 

```C++
inline discrete_domain_type PolarBSplines::Impl::full_domain () noexcept const
```



Returns the index range containing the indices of all the polar b-splines.




**Returns:**

The index range containing the indices of all the polar b-splines. 





        

<hr>



### function integrals [2/2]

```C++
template<class MemorySpace2>
void PolarBSplines::Impl::integrals (
    PolarSpline < DDim, MemorySpace2 > int_vals
) const
```



Calculate the integrals of each of the basis splines.




**Parameters:**


* `int_vals` The integrals of the basis splines. 




        

<hr>



### function nbasis 

```C++
inline std::size_t PolarBSplines::Impl::nbasis () noexcept const
```



Get the total number of basis functions.




**Returns:**

The number of basis functions. 





        

<hr>



### function operator= 

```C++
Impl & PolarBSplines::Impl::operator= (
    Impl const & x
) = default
```



A copy operator for the [**PolarBSplines**](classPolarBSplines.md).




**Parameters:**


* `x` The [**PolarBSplines**](classPolarBSplines.md) being copied.



**Returns:**

A reference to this PolarBSpline. 





        

<hr>



### function operator= 

```C++
Impl & PolarBSplines::Impl::operator= (
    Impl && x
) = default
```



A copy operator for the [**PolarBSplines**](classPolarBSplines.md) taking a temporary r-value.




**Parameters:**


* `x` The [**PolarBSplines**](classPolarBSplines.md) being copied.



**Returns:**

A reference to this PolarBSpline. 





        

<hr>



### function tensor\_bspline\_idx\_range 

_Returns the IdxRange containing the indices of the b-splines which don't traverse the singular point and can be expressed as a tensor-product of 1D b-splines._ 
```C++
inline discrete_domain_type PolarBSplines::Impl::tensor_bspline_idx_range () noexcept const
```





**Returns:**

The IdxRange containing the indices of the b-splines which don't traverse the singular point. 





        

<hr>



### function ~Impl 

```C++
PolarBSplines::Impl::~Impl () = default
```



The destructor for the [**PolarBSplines**](classPolarBSplines.md). 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_bsplines.hpp`

