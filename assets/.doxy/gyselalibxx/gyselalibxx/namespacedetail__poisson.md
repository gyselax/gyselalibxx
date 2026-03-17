

# Namespace detail\_poisson



[**Namespace List**](namespaces.md) **>** [**detail\_poisson**](namespacedetail__poisson.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION std::conditional\_t&lt; calculate\_derivs, [**DVector**](classTensor.md)&lt; typename PolarBSplinesRTheta::R::Dual, typename PolarBSplinesRTheta::Theta::Dual &gt;, void &gt; | [**get\_polar\_bspline\_vals\_and\_derivs**](#function-get_polar_bspline_vals_and_derivs) (double & val, Coord&lt; typename [**PolarBSplinesRTheta::R**](classPolarBSplines.md#typedef-r), typename [**PolarBSplinesRTheta::Theta**](classPolarBSplines.md#typedef-theta) &gt; coord, Idx&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; idx) <br>_Get the value and derivative of the specified polar bspline at the specified quadrature point._  |
|  KOKKOS\_FUNCTION IdxRange&lt; QDimRMesh, QDimThetaMesh &gt; | [**get\_quadrature\_between\_knots**](#function-get_quadrature_between_knots) (Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; [**BSplinesR**](structBSplinesR.md) &gt; &gt; start\_knot\_r, Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; [**BSplinesR**](structBSplinesR.md) &gt; &gt; end\_knot\_r, Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; [**BSplinesTheta**](structBSplinesTheta.md) &gt; &gt; start\_knot\_theta, Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; [**BSplinesTheta**](structBSplinesTheta.md) &gt; &gt; end\_knot\_theta, Idx&lt; QDimRMesh, QDimThetaMesh &gt; idx\_quad\_front) <br>_Compute the quadrature range between a provided set of knots._  |
|  KOKKOS\_FUNCTION Idx&lt; DDom &gt; | [**mod\_add**](#function-mod_add) (Idx&lt; DDom &gt; idx, IdxStep&lt; DDom &gt; idx\_step, IdxRange&lt; DDom &gt; idx\_range) <br> |
|  KOKKOS\_FUNCTION IdxStep&lt; [**BSplinesTheta**](structBSplinesTheta.md) &gt; | [**theta\_mod**](#function-theta_mod) (IdxStep&lt; [**BSplinesTheta**](structBSplinesTheta.md) &gt; idx\_theta) <br>_Calculates the modulo idx\_theta in relation to cells number along_ \(\theta\) _direction._ |
|  KOKKOS\_INLINE\_FUNCTION IdxType | [**theta\_mod**](#function-theta_mod) (IdxType idx) <br>_Calculates the index which is inside the poloidal domain using the periodicity properties._  |




























## Public Functions Documentation




### function get\_polar\_bspline\_vals\_and\_derivs 

_Get the value and derivative of the specified polar bspline at the specified quadrature point._ 
```C++
template<typename PolarBSplinesRTheta, bool calculate_derivs>
KOKKOS_FUNCTION std::conditional_t< calculate_derivs, DVector < typename PolarBSplinesRTheta::R::Dual, typename PolarBSplinesRTheta::Theta::Dual >, void > detail_poisson::get_polar_bspline_vals_and_derivs (
    double & val,
    Coord< typename PolarBSplinesRTheta::R , typename PolarBSplinesRTheta::Theta > coord,
    Idx< PolarBSplinesRTheta > idx
) 
```



This method calculates the value and the derivatives of polar bsplines. It is templated by calculate\_derivs to avoid code duplication between get\_polar\_bspline\_vals\_and\_derivs and get\_polar\_bspline\_vals. The calling method should not need to use the template parameter.




**Template parameters:**


* [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) The Polar B-Splines to be used. 
* `calculate_derivs` Returns the derivatives at the coordinate if true. If false, only the value is calculated.



**Parameters:**


* `val` The value of the specified polar bspline at the specified point. 
* `coord` The coordinate where the value of the polar bspline should be calculated. 
* `idx` The polar bspline of interest. 



**Returns:**

The derivative of the polar bspline (only returned if calculate\_derivs is true). 





        

<hr>



### function get\_quadrature\_between\_knots 

_Compute the quadrature range between a provided set of knots._ 
```C++
template<typename QDimRMesh, typename QDimThetaMesh, typename BSplinesR, typename BSplinesTheta>
KOKKOS_FUNCTION IdxRange< QDimRMesh, QDimThetaMesh > detail_poisson::get_quadrature_between_knots (
    Idx< ddc::knot_discrete_dimension_t< BSplinesR > > start_knot_r,
    Idx< ddc::knot_discrete_dimension_t< BSplinesR > > end_knot_r,
    Idx< ddc::knot_discrete_dimension_t< BSplinesTheta > > start_knot_theta,
    Idx< ddc::knot_discrete_dimension_t< BSplinesTheta > > end_knot_theta,
    Idx< QDimRMesh, QDimThetaMesh > idx_quad_front
) 
```



Compute the range of quadrature points which are found between a set of knots in both the radial and poloidal directions. In order to return a contiguous range the result may include indices which are outside the domain. A modulo operator should be applied before using the indices.




**Template parameters:**


* `QDimRMesh` [**Quadrature**](classQuadrature.md) mesh in radial direction. 
* `QDimThetaMesh` [**Quadrature**](classQuadrature.md) mesh in poloidal direction. 
* [**BSplinesR**](structBSplinesR.md) Splines in radial direction. 
* [**BSplinesTheta**](structBSplinesTheta.md) Splines in poloidal direction.



**Parameters:**


* `start_knot_r` The index of the knot describing the lower bound of the domain of interest in the radial direction. 
* `end_knot_r` The index of the knot describing the upper bound of the domain of interest in the radial direction. 
* `start_knot_theta` The index of the knot describing the lower bound of the domain of interest in the poloidal direction. 
* `end_knot_theta` The index of the knot describing the upper bound of the domain of interest in the poloidal direction. 
* `idx_quad_front` The first index of the index range of the quadrature points. 



**Returns:**

The range of quadrature points in the specified domain. 





        

<hr>



### function mod\_add 

```C++
template<typename DDom, std::enable_if_t< DDom::continuous_dimension_type::PERIODIC, bool >>
KOKKOS_FUNCTION Idx< DDom > detail_poisson::mod_add (
    Idx< DDom > idx,
    IdxStep< DDom > idx_step,
    IdxRange< DDom > idx_range
) 
```




<hr>



### function theta\_mod 

_Calculates the modulo idx\_theta in relation to cells number along_ \(\theta\) _direction._
```C++
template<typename BSplinesTheta>
KOKKOS_FUNCTION IdxStep< BSplinesTheta > detail_poisson::theta_mod (
    IdxStep< BSplinesTheta > idx_theta
) 
```





**Template parameters:**


* [**BSplinesTheta**](structBSplinesTheta.md) Periodic B-Splines in poloidal direction.



**Parameters:**


* `idx_theta` \(\theta\) index.



**Returns:**

The corresponding indice modulo \(\theta\) direction cells number 





        

<hr>



### function theta\_mod 

_Calculates the index which is inside the poloidal domain using the periodicity properties._ 
```C++
template<typename BSplinesTheta, class IdxType>
KOKKOS_INLINE_FUNCTION IdxType detail_poisson::theta_mod (
    IdxType idx
) 
```





**Template parameters:**


* [**BSplinesTheta**](structBSplinesTheta.md) Periodic B-Splines in poloidal direction. 
* `IdxType` The Type of the input index.



**Parameters:**


* `idx` A multi-dimensional index including the polar bspline index.



**Returns:**

The corresponding index inside the domain. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polarpoissonlikeassembler.hpp`

