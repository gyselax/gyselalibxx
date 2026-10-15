

# Struct detail\_poisson::TensorBSplineStash

**template &lt;typename QDimRMesh, typename QDimThetaMesh, typename [**BSplinesR**](structBSplinesR.md), typename [**BSplinesTheta**](structBSplinesTheta.md)&gt;**



[**ClassList**](annotated.md) **>** [**detail\_poisson**](namespacedetail__poisson.md) **>** [**TensorBSplineStash**](structdetail__poisson_1_1TensorBSplineStash.md)



_The values and first derivatives of the non-zero 1D B-splines at the 1D quadrature points._ [More...](#detailed-description)

* `#include <polar_spline_fem_poisson_like_assembler.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  ConstField&lt; Idx&lt; [**BSplinesR**](structBSplinesR.md) &gt;, IdxRange&lt; QDimRMesh &gt; &gt; | [**first\_basis\_r**](#variable-first_basis_r)  <br>_The first non-zero radial B-spline at each radial quadrature point._  |
|  ConstField&lt; Idx&lt; [**BSplinesTheta**](structBSplinesTheta.md) &gt;, IdxRange&lt; QDimThetaMesh &gt; &gt; | [**first\_basis\_theta**](#variable-first_basis_theta)  <br>_The first non-zero poloidal B-spline at each poloidal quadrature point._  |
|  DConstField&lt; IdxRange&lt; QDimRMesh, [**BSplinesR**](structBSplinesR.md), ddc::Deriv&lt; [**R**](structR.md) &gt; &gt; &gt; | [**vals\_r**](#variable-vals_r)  <br>_The values (derivative order 0) and derivatives (derivative order 1) of the radial B-splines._  |
|  DConstField&lt; IdxRange&lt; QDimThetaMesh, [**BSplinesTheta**](structBSplinesTheta.md), ddc::Deriv&lt; [**Theta**](structTheta.md) &gt; &gt; &gt; | [**vals\_theta**](#variable-vals_theta)  <br>_The values (derivative order 0) and derivatives (derivative order 1) of the poloidal B-splines._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION double | [**vals**](#function-vals) (Idx&lt; QDimRMesh, QDimThetaMesh &gt; idx\_quad, Idx&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; idx) const<br>_Get the value of a tensor product B-spline at a quadrature point._  |
|  KOKKOS\_FUNCTION std::tuple&lt; double, [**DVector**](classTensor.md)&lt; [**R\_cov**](structR__cov.md), [**Theta\_cov**](structTheta__cov.md) &gt; &gt; | [**vals\_and\_derivs**](#function-vals_and_derivs) (Idx&lt; QDimRMesh, QDimThetaMesh &gt; idx\_quad, Idx&lt; [**BSplinesR**](structBSplinesR.md), [**BSplinesTheta**](structBSplinesTheta.md) &gt; idx) const<br>_Get the value and the derivatives of a tensor product B-spline at a quadrature point._  |




























## Detailed Description


The tensor product B-splines and the quadrature grid both have a tensor product structure. The value of a tensor product B-spline at a quadrature point can therefore be computed from the values of the 1D B-splines at the 1D quadrature points. The memory required to store these values is O(N\_r + N\_theta).


The B-spline index stored in the stash is the offset from the first non-zero B-spline at the quadrature point.




**Template parameters:**


* `QDimRMesh` [**Quadrature**](classQuadrature.md) mesh in radial direction. 
* `QDimThetaMesh` [**Quadrature**](classQuadrature.md) mesh in poloidal direction. 
* [**BSplinesR**](structBSplinesR.md) Splines in radial direction. 
* [**BSplinesTheta**](structBSplinesTheta.md) Splines in poloidal direction. 




    
## Public Attributes Documentation




### variable first\_basis\_r 

_The first non-zero radial B-spline at each radial quadrature point._ 
```C++
ConstField<Idx<BSplinesR>, IdxRange<QDimRMesh> > detail_poisson::TensorBSplineStash< QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta >::first_basis_r;
```




<hr>



### variable first\_basis\_theta 

_The first non-zero poloidal B-spline at each poloidal quadrature point._ 
```C++
ConstField<Idx<BSplinesTheta>, IdxRange<QDimThetaMesh> > detail_poisson::TensorBSplineStash< QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta >::first_basis_theta;
```




<hr>



### variable vals\_r 

_The values (derivative order 0) and derivatives (derivative order 1) of the radial B-splines._ 
```C++
DConstField<IdxRange<QDimRMesh, BSplinesR, ddc::Deriv<R> > > detail_poisson::TensorBSplineStash< QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta >::vals_r;
```




<hr>



### variable vals\_theta 

_The values (derivative order 0) and derivatives (derivative order 1) of the poloidal B-splines._ 
```C++
DConstField<IdxRange<QDimThetaMesh, BSplinesTheta, ddc::Deriv<Theta> > > detail_poisson::TensorBSplineStash< QDimRMesh, QDimThetaMesh, BSplinesR, BSplinesTheta >::vals_theta;
```




<hr>
## Public Functions Documentation




### function vals 

_Get the value of a tensor product B-spline at a quadrature point._ 
```C++
inline KOKKOS_FUNCTION double detail_poisson::TensorBSplineStash::vals (
    Idx< QDimRMesh, QDimThetaMesh > idx_quad,
    Idx< BSplinesR , BSplinesTheta > idx
) const
```





**Parameters:**


* `idx` The 2D index of the tensor product B-spline. 
* `idx_quad` The index of the quadrature point. 



**Returns:**

The value of the B-spline at the quadrature point. 





        

<hr>



### function vals\_and\_derivs 

_Get the value and the derivatives of a tensor product B-spline at a quadrature point._ 
```C++
inline KOKKOS_FUNCTION std::tuple< double, DVector < R_cov , Theta_cov > > detail_poisson::TensorBSplineStash::vals_and_derivs (
    Idx< QDimRMesh, QDimThetaMesh > idx_quad,
    Idx< BSplinesR , BSplinesTheta > idx
) const
```





**Parameters:**


* `idx_quad` The index of the quadrature point. 
* `idx` The 2D index of the tensor product B-spline. 



**Returns:**

The value and derivatives of the B-spline at the quadrature point. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polar_spline_fem_poisson_like_assembler.hpp`

