

# Namespace detail\_poisson



[**Namespace List**](namespaces.md) **>** [**detail\_poisson**](namespacedetail__poisson.md)




















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**TensorBSplineStash**](structdetail__poisson_1_1TensorBSplineStash.md) &lt;typename QDim1Mesh, typename QDim2Mesh, typename BSplines1, typename BSplines2&gt;<br>_The values and first derivatives of the non-zero 1D B-splines at the 1D quadrature points._  |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  Kokkos::View&lt; double \*\*, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace &gt; | [**as\_batched\_vector\_view**](#function-as_batched_vector_view) (DField&lt; IdxRange&lt; Dims... &gt; &gt; field) <br>_Get a 2D Kokkos view of a 1D field with N batch dimensions._  |
|  KOKKOS\_FUNCTION std::array&lt; IdxRange&lt; BSplines &gt;, 2 &gt; | [**get\_overlapping\_bsplines\_1d**](#function-get_overlapping_bsplines_1d) (Idx&lt; BSplines &gt; idx\_test, IdxRange&lt; BSplines &gt; idx\_range\_basis) <br>_Get the 1D B-splines in a basis whose supports overlap the support of a given B-spline._  |
|  KOKKOS\_INLINE\_FUNCTION IdxRange&lt; ddc::knot\_discrete\_dimension\_t&lt; BSplines &gt; &gt; | [**get\_overlapping\_knots\_1d**](#function-get_overlapping_knots_1d) (Idx&lt; BSplines &gt; idx\_test, Idx&lt; BSplines &gt; idx\_trial) <br>_Compute the range of knots between which two 1D B-splines are both non-zero._  |
|  KOKKOS\_FUNCTION std::conditional\_t&lt; calculate\_derivs, std::tuple&lt; double, DVector&lt; typename PolarBSplinesRTheta::R::Dual, typename PolarBSplinesRTheta::Theta::Dual &gt; &gt;, double &gt; | [**get\_polar\_bspline\_vals\_and\_derivs**](#function-get_polar_bspline_vals_and_derivs) (IdxQuad idx\_quad, Idx&lt; [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md) &gt; idx) <br>_Get the value and derivatives of a polar bspline which covers the singular point at the specified quadrature point._  |
|  KOKKOS\_FUNCTION IdxRange&lt; QDim1Mesh, QDim2Mesh &gt; | [**get\_quadrature\_between\_knots**](#function-get_quadrature_between_knots) (Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; BSplines1 &gt; &gt; start\_knot\_1, Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; BSplines1 &gt; &gt; end\_knot\_1, Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; BSplines2 &gt; &gt; start\_knot\_2, Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; BSplines2 &gt; &gt; end\_knot\_2, Idx&lt; QDim1Mesh, QDim2Mesh &gt; idx\_quad\_front) <br>_Compute the quadrature range between a provided set of knots._  |
|  KOKKOS\_FUNCTION IdxRange&lt; QDimMesh &gt; | [**get\_quadrature\_between\_knots\_1d**](#function-get_quadrature_between_knots_1d) (Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; BSplines &gt; &gt; start\_knot, Idx&lt; ddc::knot\_discrete\_dimension\_t&lt; BSplines &gt; &gt; end\_knot, Idx&lt; QDimMesh &gt; idx\_quad\_front) <br>_Compute the 1D quadrature range between a provided pair of knots._  |
|  KOKKOS\_FUNCTION IdxRange&lt; QDim1Mesh, QDim2Mesh &gt; | [**get\_support\_quadrature**](#function-get_support_quadrature) (Idx&lt; BSplines1, BSplines2 &gt; idx, Idx&lt; QDim1Mesh, QDim2Mesh &gt; idx\_quad\_front) <br>_Compute the range of quadrature points in the support of a tensor product B-spline._  |
|  KOKKOS\_FUNCTION IdxRange&lt; QDimMesh &gt; | [**get\_support\_quadrature\_1d**](#function-get_support_quadrature_1d) (Idx&lt; BSplines &gt; idx, Idx&lt; QDimMesh &gt; idx\_quad\_front) <br>_Compute the range of quadrature points in the support of a 1D B-spline._  |
|  KOKKOS\_FUNCTION std::tuple&lt; double, DVector&lt; typename QDim1Mesh::continuous\_dimension\_type, typename QDim2Mesh::continuous\_dimension\_type &gt; &gt; | [**get\_test\_weights**](#function-get_test_weights) (double basis\_val\_test, DVector&lt; typename QDim1Mesh::continuous\_dimension\_type::Dual, typename QDim2Mesh::continuous\_dimension\_type::Dual &gt; basis\_derivs\_test, Idx&lt; QDim1Mesh, QDim2Mesh &gt; idx\_quad, IdxBatch idx\_batch, CoeffAlpha const & coeff\_alpha, CoeffBeta const & coeff\_beta, MetricTensorType const & metric\_tensor, DField&lt; IdxRange&lt; QDim1Mesh, QDim2Mesh &gt; &gt; int\_volume) <br>_Computes the test\_stiffness\_weight and test\_mass\_weight for a given quadrature point._  |
|  KOKKOS\_FUNCTION Idx&lt; DDom &gt; | [**mod\_add**](#function-mod_add) (Idx&lt; DDom &gt; idx, IdxStep&lt; DDom &gt; idx\_step, IdxRange&lt; DDom &gt; idx\_range) <br>_Add an idx\_step to an idx in such a way that the result remains in the idx\_range._  |
|  KOKKOS\_INLINE\_FUNCTION int | [**n\_nonzeros\_in\_rows**](#function-n_nonzeros_in_rows) (int i, int nbasis) <br>_Get the number of non-zero elements in the first i rows of a 1D stiffness matrix._  |
|  KOKKOS\_FUNCTION IdxStep&lt; BSplines &gt; | [**periodic\_mod**](#function-periodic_mod) (IdxStep&lt; BSplines &gt; idx\_step) <br>_Calculates the modulo idx\_step in relation to the number of basis functions of periodic B-splines._  |
|  KOKKOS\_INLINE\_FUNCTION IdxType | [**periodic\_mod**](#function-periodic_mod) (IdxType idx) <br>_Calculates the index which is inside the domain of periodic B-splines using the periodicity properties._  |
|  KOKKOS\_INLINE\_FUNCTION Idx&lt; QDim1Mesh, QDim2Mesh &gt; | [**periodic\_wrap**](#function-periodic_wrap) (Idx&lt; QDim1Mesh, QDim2Mesh &gt; idx\_quad, IdxRange&lt; QDim1Mesh, QDim2Mesh &gt; const & target\_range, IdxRange&lt; QDim1Mesh, QDim2Mesh &gt; const & full\_quad\_idx\_range) <br>_Get the periodic image of a quadrature index which is found in a target range._  |




























## Public Functions Documentation




### function as\_batched\_vector\_view 

_Get a 2D Kokkos view of a 1D field with N batch dimensions._ 
```C++
template<class... Dims>
Kokkos::View< double **, Kokkos::LayoutRight, Kokkos::DefaultExecutionSpace > detail_poisson::as_batched_vector_view (
    DField< IdxRange< Dims... > > field
) 
```



The returned view shares the memory of the field. The first dimension of the view iterates over all the batch dimensions (in layout\_right order), the second dimension iterates over the last dimension of the field.




**Parameters:**


* `field` A contiguous layout\_right field whose last dimension is the dimension of the matrix equation and whose other dimensions are batch dimensions.



**Returns:**

A 2D Kokkos view of shape (number of batches, matrix size). 





        

<hr>



### function get\_overlapping\_bsplines\_1d 

_Get the 1D B-splines in a basis whose supports overlap the support of a given B-spline._ 
```C++
template<typename BSplines>
KOKKOS_FUNCTION std::array< IdxRange< BSplines >, 2 > detail_poisson::get_overlapping_bsplines_1d (
    Idx< BSplines > idx_test,
    IdxRange< BSplines > idx_range_basis
) 
```



The result is returned as 2 contiguous index ranges. The indices in the first range are all smaller than the indices in the second range. This allows the B-splines to be iterated over in order, even when the overlap wraps around the boundary of a periodic basis. The second range is empty if the overlap does not wrap around the boundary.




**Template parameters:**


* `BSplines` The 1D B-splines.



**Parameters:**


* `idx_test` The index of the B-spline (the test function). 
* `idx_range_basis` The index range of the B-splines in the basis. In a periodic dimension this must contain all the B-splines. 



**Returns:**

The 2 sorted index ranges containing all the B-splines in the basis which overlap the given B-spline. 





        

<hr>



### function get\_overlapping\_knots\_1d 

_Compute the range of knots between which two 1D B-splines are both non-zero._ 
```C++
template<typename BSplines>
KOKKOS_INLINE_FUNCTION IdxRange< ddc::knot_discrete_dimension_t< BSplines > > detail_poisson::get_overlapping_knots_1d (
    Idx< BSplines > idx_test,
    Idx< BSplines > idx_trial
) 
```



If the B-splines are periodic then, in order to return a contiguous range, the result may include knots beyond the end of the domain.




**Template parameters:**


* `BSplines` The 1D B-splines.



**Parameters:**


* `idx_test` The index of the first B-spline (the test function). 
* `idx_trial` The index of the second B-spline (the trial function). 



**Returns:**

The range of knots describing the cells where both B-splines are non-zero. 





        

<hr>



### function get\_polar\_bspline\_vals\_and\_derivs 

_Get the value and derivatives of a polar bspline which covers the singular point at the specified quadrature point._ 
```C++
template<bool calculate_derivs, class PolarBSplinesRTheta, class IdxQuad>
KOKKOS_FUNCTION std::conditional_t< calculate_derivs, std::tuple< double, DVector< typename PolarBSplinesRTheta::R::Dual, typename PolarBSplinesRTheta::Theta::Dual > >, double > detail_poisson::get_polar_bspline_vals_and_derivs (
    IdxQuad idx_quad,
    Idx< PolarBSplinesRTheta > idx
) 
```



The polar bsplines which cover the singular point are linear combinations of tensor product bsplines so they are evaluated using PolarBSplines::eval\_basis\_and\_n\_derivs. The polar bsplines which do not cover the singular point should be passed as a 2D index (see the overload below).




**Template parameters:**


* `calculate_derivs` If true, the value and the derivatives are calculated. If false, only the value is calculated.



**Parameters:**


* `idx_quad` The quadrature point where the polar bspline should be evaluated. 
* `idx` The index of a polar bspline which covers the singular point. 



**Returns:**

The value of the polar bspline if calculate\_derivs is false. Otherwise a tuple containing the value and the derivatives of the polar bspline. 





        

<hr>



### function get\_quadrature\_between\_knots 

_Compute the quadrature range between a provided set of knots._ 
```C++
template<typename QDim1Mesh, typename QDim2Mesh, typename BSplines1, typename BSplines2>
KOKKOS_FUNCTION IdxRange< QDim1Mesh, QDim2Mesh > detail_poisson::get_quadrature_between_knots (
    Idx< ddc::knot_discrete_dimension_t< BSplines1 > > start_knot_1,
    Idx< ddc::knot_discrete_dimension_t< BSplines1 > > end_knot_1,
    Idx< ddc::knot_discrete_dimension_t< BSplines2 > > start_knot_2,
    Idx< ddc::knot_discrete_dimension_t< BSplines2 > > end_knot_2,
    Idx< QDim1Mesh, QDim2Mesh > idx_quad_front
) 
```



Compute the range of quadrature points which are found between a set of knots in both dimensions. The index of the front of the returned range in the second dimension is always inside the domain. However, in order to return a contiguous range, the result may include indices beyond the end of the domain in the second dimension. A modulo operator should be applied before using the indices. See get\_quadrature\_between\_knots\_1d for implementation details.




**Template parameters:**


* `QDim1Mesh` [**Quadrature**](classQuadrature.md) mesh in the first dimension. 
* `QDim2Mesh` [**Quadrature**](classQuadrature.md) mesh in the second dimension. 
* `BSplines1` Splines in the first dimension. 
* `BSplines2` Splines in the second dimension.



**Parameters:**


* `start_knot_1` The index of the knot describing the lower bound of the domain of interest in the first dimension. 
* `end_knot_1` The index of the knot describing the upper bound of the domain of interest in the first dimension. 
* `start_knot_2` The index of the knot describing the lower bound of the domain of interest in the second dimension. 
* `end_knot_2` The index of the knot describing the upper bound of the domain of interest in the second dimension. 
* `idx_quad_front` The first index of the index range of the quadrature points. 



**Returns:**

The range of quadrature points in the specified domain. 





        

<hr>



### function get\_quadrature\_between\_knots\_1d 

_Compute the 1D quadrature range between a provided pair of knots._ 
```C++
template<typename QDimMesh, typename BSplines>
KOKKOS_FUNCTION IdxRange< QDimMesh > detail_poisson::get_quadrature_between_knots_1d (
    Idx< ddc::knot_discrete_dimension_t< BSplines > > start_knot,
    Idx< ddc::knot_discrete_dimension_t< BSplines > > end_knot,
    Idx< QDimMesh > idx_quad_front
) 
```



Compute the range of quadrature points which are found between two knots. If the B-splines are periodic then the front of the returned range is always inside the domain. However, in order to return a contiguous range, the result may include indices beyond the end of the domain. A modulo operator should be applied before using the indices.




**Template parameters:**


* `QDimMesh` [**Quadrature**](classQuadrature.md) mesh. 
* `BSplines` The 1D B-splines.



**Parameters:**


* `start_knot` The index of the knot describing the lower bound of the domain of interest. 
* `end_knot` The index of the knot describing the upper bound of the domain of interest. 
* `idx_quad_front` The first index of the index range of the quadrature points. 



**Returns:**

The range of quadrature points in the specified domain. 





        

<hr>



### function get\_support\_quadrature 

_Compute the range of quadrature points in the support of a tensor product B-spline._ 
```C++
template<typename QDim1Mesh, typename QDim2Mesh, typename BSplines1, typename BSplines2>
KOKKOS_FUNCTION IdxRange< QDim1Mesh, QDim2Mesh > detail_poisson::get_support_quadrature (
    Idx< BSplines1, BSplines2 > idx,
    Idx< QDim1Mesh, QDim2Mesh > idx_quad_front
) 
```



In order to return a contiguous range the result may include indices which are outside the domain in periodic dimensions. A modulo operator should be applied before using the indices.




**Template parameters:**


* `QDim1Mesh` [**Quadrature**](classQuadrature.md) mesh in the first dimension. 
* `QDim2Mesh` [**Quadrature**](classQuadrature.md) mesh in the second dimension. 
* `BSplines1` Splines in the first dimension. 
* `BSplines2` Splines in the second dimension.



**Parameters:**


* `idx` The 2D index of the tensor product B-spline. 
* `idx_quad_front` The first index of the index range of the quadrature points. 



**Returns:**

The range of quadrature points in the support of the B-spline. 





        

<hr>



### function get\_support\_quadrature\_1d 

_Compute the range of quadrature points in the support of a 1D B-spline._ 
```C++
template<typename QDimMesh, typename BSplines>
KOKKOS_FUNCTION IdxRange< QDimMesh > detail_poisson::get_support_quadrature_1d (
    Idx< BSplines > idx,
    Idx< QDimMesh > idx_quad_front
) 
```



If the B-splines are periodic then, in order to return a contiguous range, the result may include indices which are outside the domain. A modulo operator should be applied before using the indices.




**Template parameters:**


* `QDimMesh` [**Quadrature**](classQuadrature.md) mesh. 
* `BSplines` The 1D B-splines.



**Parameters:**


* `idx` The index of the B-spline. 
* `idx_quad_front` The first index of the index range of the quadrature points. 



**Returns:**

The range of quadrature points in the support of the B-spline. 





        

<hr>



### function get\_test\_weights 

_Computes the test\_stiffness\_weight and test\_mass\_weight for a given quadrature point._ 
```C++
template<class QDim1Mesh, class QDim2Mesh, class IdxBatch, class CoeffAlpha, class CoeffBeta, concepts::MetricTensor MetricTensorType>
KOKKOS_FUNCTION std::tuple< double, DVector< typename QDim1Mesh::continuous_dimension_type, typename QDim2Mesh::continuous_dimension_type > > detail_poisson::get_test_weights (
    double basis_val_test,
    DVector< typename QDim1Mesh::continuous_dimension_type::Dual, typename QDim2Mesh::continuous_dimension_type::Dual > basis_derivs_test,
    Idx< QDim1Mesh, QDim2Mesh > idx_quad,
    IdxBatch idx_batch,
    CoeffAlpha const & coeff_alpha,
    CoeffBeta const & coeff_beta,
    MetricTensorType const & metric_tensor,
    DField< IdxRange< QDim1Mesh, QDim2Mesh > > int_volume
) 
```



The weak integral element is \(w (\alpha \nabla B_i \cdot G^{-1} \nabla B_j + \beta B_i B_j)\) where \(w\) is the integral volume, \(B_i\) is the test function and \(B_j\) is the trial function. This function computes \(w \alpha G^{-1} \nabla B_i\) and \(w \beta B_i\) so that the weak integral element can be computed from the values and derivatives of the trial function alone.




**Template parameters:**


* `QDim1Mesh` [**Quadrature**](classQuadrature.md) mesh in the first dimension. 
* `QDim2Mesh` [**Quadrature**](classQuadrature.md) mesh in the second dimension.



**Parameters:**


* `basis_val_test` The value of the test basis function at the quadrature point. 
* `basis_derivs_test` The gradient of the test basis function at the quadrature point. 
* `idx_quad` The index for the point in the quadrature scheme. 
* `idx_batch` The index of the batch for which the element is being calculated. 
* `coeff_alpha` A callable with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)` returning \(\alpha\) at the given batch and coordinate. 
* `coeff_beta` A callable with signature `double operator()(IdxBatch, Coord<Dim1, Dim2>)` returning \(\beta\) at the given batch and coordinate. 
* `metric_tensor` The metric tensor associated with the mapping from the logical domain to the physical domain where the equation is defined. 
* `int_volume` The integral volume associated with each point used in the quadrature scheme. 



**Returns:**

A tuple containing
* test\_mass\_weight, the quantity \(w \beta B_i\) which multiplies the value of the trial function,
* test\_stiffness\_weight, the quantity \(w \alpha G^{-1} \nabla B_i\) which multiplies the gradient of the trial function. 







        

<hr>



### function mod\_add 

_Add an idx\_step to an idx in such a way that the result remains in the idx\_range._ 
```C++
template<typename DDom, std::enable_if_t< DDom::continuous_dimension_type::PERIODIC, bool >>
KOKKOS_FUNCTION Idx< DDom > detail_poisson::mod_add (
    Idx< DDom > idx,
    IdxStep< DDom > idx_step,
    IdxRange< DDom > idx_range
) 
```




<hr>



### function n\_nonzeros\_in\_rows 

_Get the number of non-zero elements in the first i rows of a 1D stiffness matrix._ 
```C++
template<typename BSplines>
KOKKOS_INLINE_FUNCTION int detail_poisson::n_nonzeros_in_rows (
    int i,
    int nbasis
) 
```



Row j of the 1D stiffness matrix contains one non-zero element for each B-spline whose support overlaps the support of the j-th B-spline (including the j-th B-spline itself). In a periodic basis all rows contain the same number of non-zero elements. In a non-periodic basis this number decreases near the boundaries. The result is therefore the cumulative quantity required to build the nnz data structure of a CSR matrix. When i == nbasis the result is the total number of non-zero elements in the 1D stiffness matrix.




**Template parameters:**


* `BSplines` The B-splines in the dimension of interest.



**Parameters:**


* `i` The number of rows (counted from the start of the basis) whose non-zero elements are counted. 
* `nbasis` The number of B-splines in the basis. 



**Returns:**

The number of non-zero elements in rows 0 to i-1 of the 1D stiffness matrix. 





        

<hr>



### function periodic\_mod 

_Calculates the modulo idx\_step in relation to the number of basis functions of periodic B-splines._ 
```C++
template<typename BSplines>
KOKKOS_FUNCTION IdxStep< BSplines > detail_poisson::periodic_mod (
    IdxStep< BSplines > idx_step
) 
```





**Template parameters:**


* `BSplines` Periodic B-Splines.



**Parameters:**


* `idx_step` A step between two B-spline indices.



**Returns:**

The corresponding step modulo the number of basis functions. 





        

<hr>



### function periodic\_mod 

_Calculates the index which is inside the domain of periodic B-splines using the periodicity properties._ 
```C++
template<typename BSplines, class IdxType>
KOKKOS_INLINE_FUNCTION IdxType detail_poisson::periodic_mod (
    IdxType idx
) 
```





**Template parameters:**


* `BSplines` Periodic B-Splines. 
* `IdxType` The Type of the input index.



**Parameters:**


* `idx` A multi-dimensional index including the B-spline index.



**Returns:**

The corresponding index inside the domain. 





        

<hr>



### function periodic\_wrap 

_Get the periodic image of a quadrature index which is found in a target range._ 
```C++
template<typename BSplines1, typename BSplines2, typename QDim1Mesh, typename QDim2Mesh>
KOKKOS_INLINE_FUNCTION Idx< QDim1Mesh, QDim2Mesh > detail_poisson::periodic_wrap (
    Idx< QDim1Mesh, QDim2Mesh > idx_quad,
    IdxRange< QDim1Mesh, QDim2Mesh > const & target_range,
    IdxRange< QDim1Mesh, QDim2Mesh > const & full_quad_idx_range
) 
```



The quadrature ranges returned by get\_support\_quadrature and get\_quadrature\_between\_knots are contiguous so they may extend past the end of the domain in periodic dimensions. This function shifts an index by a period in each periodic dimension where it is outside the target range. It can be used to move an index from such a range into the domain (target\_range = full\_quad\_idx\_range) or from the domain into such a range.




**Template parameters:**


* `BSplines1` Splines in the first dimension. 
* `BSplines2` Splines in the second dimension. 
* `QDim1Mesh` [**Quadrature**](classQuadrature.md) mesh in the first dimension. 
* `QDim2Mesh` [**Quadrature**](classQuadrature.md) mesh in the second dimension.



**Parameters:**


* `idx_quad` The index of a quadrature point. 
* `target_range` The range in which the returned index should be found. 
* `full_quad_idx_range` The index range of all the quadrature points (one period). 



**Returns:**

The equivalent index inside target\_range. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/fem_poisson_assembler_utils.hpp`

