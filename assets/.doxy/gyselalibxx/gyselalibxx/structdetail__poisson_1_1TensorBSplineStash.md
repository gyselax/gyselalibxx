

# Struct detail\_poisson::TensorBSplineStash

**template &lt;typename QDim1Mesh, typename QDim2Mesh, typename BSplines1, typename BSplines2&gt;**



[**ClassList**](annotated.md) **>** [**detail\_poisson**](namespacedetail__poisson.md) **>** [**TensorBSplineStash**](structdetail__poisson_1_1TensorBSplineStash.md)



_The values and first derivatives of the non-zero 1D B-splines at the 1D quadrature points._ [More...](#detailed-description)

* `#include <fem_poisson_assembler_utils.hpp>`





















## Public Attributes

| Type | Name |
| ---: | :--- |
|  ConstField&lt; Idx&lt; BSplines1 &gt;, IdxRange&lt; QDim1Mesh &gt; &gt; | [**first\_basis\_1**](#variable-first_basis_1)  <br>_The first non-zero B-spline at each quadrature point in the first dimension._  |
|  ConstField&lt; Idx&lt; BSplines2 &gt;, IdxRange&lt; QDim2Mesh &gt; &gt; | [**first\_basis\_2**](#variable-first_basis_2)  <br>_The first non-zero B-spline at each quadrature point in the second dimension._  |
|  DConstField&lt; IdxRange&lt; QDim1Mesh, BSplines1, ddc::Deriv&lt; Dim1 &gt; &gt; &gt; | [**vals\_1**](#variable-vals_1)  <br>_The values (derivative order 0) and derivatives (derivative order 1) of the B-splines in the first dimension._  |
|  DConstField&lt; IdxRange&lt; QDim2Mesh, BSplines2, ddc::Deriv&lt; Dim2 &gt; &gt; &gt; | [**vals\_2**](#variable-vals_2)  <br>_The values (derivative order 0) and derivatives (derivative order 1) of the B-splines in the second dimension._  |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION double | [**vals**](#function-vals) (Idx&lt; QDim1Mesh, QDim2Mesh &gt; idx\_quad, Idx&lt; BSplines1, BSplines2 &gt; idx) const<br>_Get the value of a tensor product B-spline at a quadrature point._  |
|  KOKKOS\_FUNCTION std::tuple&lt; double, DVector&lt; Dim1\_cov, Dim2\_cov &gt; &gt; | [**vals\_and\_derivs**](#function-vals_and_derivs) (Idx&lt; QDim1Mesh, QDim2Mesh &gt; idx\_quad, Idx&lt; BSplines1, BSplines2 &gt; idx) const<br>_Get the value and the derivatives of a tensor product B-spline at a quadrature point._  |




























## Detailed Description


The tensor product B-splines and the quadrature grid both have a tensor product structure. The value of a tensor product B-spline at a quadrature point can therefore be computed from the values of the 1D B-splines at the 1D quadrature points. The memory required to store these values is O(N\_1 + N\_2).


The B-spline index stored in the stash is the offset from the first non-zero B-spline at the quadrature point.




**Template parameters:**


* `QDim1Mesh` [**Quadrature**](classQuadrature.md) mesh in the first dimension. 
* `QDim2Mesh` [**Quadrature**](classQuadrature.md) mesh in the second dimension. 
* `BSplines1` Splines in the first dimension. 
* `BSplines2` Splines in the second dimension. 




    
## Public Attributes Documentation




### variable first\_basis\_1 

_The first non-zero B-spline at each quadrature point in the first dimension._ 
```C++
ConstField<Idx<BSplines1>, IdxRange<QDim1Mesh> > detail_poisson::TensorBSplineStash< QDim1Mesh, QDim2Mesh, BSplines1, BSplines2 >::first_basis_1;
```




<hr>



### variable first\_basis\_2 

_The first non-zero B-spline at each quadrature point in the second dimension._ 
```C++
ConstField<Idx<BSplines2>, IdxRange<QDim2Mesh> > detail_poisson::TensorBSplineStash< QDim1Mesh, QDim2Mesh, BSplines1, BSplines2 >::first_basis_2;
```




<hr>



### variable vals\_1 

_The values (derivative order 0) and derivatives (derivative order 1) of the B-splines in the first dimension._ 
```C++
DConstField<IdxRange<QDim1Mesh, BSplines1, ddc::Deriv<Dim1> > > detail_poisson::TensorBSplineStash< QDim1Mesh, QDim2Mesh, BSplines1, BSplines2 >::vals_1;
```




<hr>



### variable vals\_2 

_The values (derivative order 0) and derivatives (derivative order 1) of the B-splines in the second dimension._ 
```C++
DConstField<IdxRange<QDim2Mesh, BSplines2, ddc::Deriv<Dim2> > > detail_poisson::TensorBSplineStash< QDim1Mesh, QDim2Mesh, BSplines1, BSplines2 >::vals_2;
```




<hr>
## Public Functions Documentation




### function vals 

_Get the value of a tensor product B-spline at a quadrature point._ 
```C++
inline KOKKOS_FUNCTION double detail_poisson::TensorBSplineStash::vals (
    Idx< QDim1Mesh, QDim2Mesh > idx_quad,
    Idx< BSplines1, BSplines2 > idx
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
inline KOKKOS_FUNCTION std::tuple< double, DVector< Dim1_cov, Dim2_cov > > detail_poisson::TensorBSplineStash::vals_and_derivs (
    Idx< QDim1Mesh, QDim2Mesh > idx_quad,
    Idx< BSplines1, BSplines2 > idx
) const
```





**Parameters:**


* `idx_quad` The index of the quadrature point. 
* `idx` The 2D index of the tensor product B-spline. 



**Returns:**

The value and derivatives of the B-spline at the quadrature point. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/fem_poisson_assembler_utils.hpp`

