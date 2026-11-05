

# Class DerivReconstructionCoeffs

**template &lt;class EdgePerpGrid1Type, class EdgePerpGrid2Type, class BSplinesPerp1Type, class BSplinesPerp2Type, Extremity extremity\_1, Extremity extremity\_2&gt;**



[**ClassList**](annotated.md) **>** [**DerivReconstructionCoeffs**](classDerivReconstructionCoeffs.md)



_Compute the coefficients a, b and c of the interface derivative reconstruction method._ [More...](#detailed-description)

* `#include <deriv_reconstruction_coeffs.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**DerivReconstructionCoeffs**](#function-derivreconstructioncoeffs-14) ([**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) const & idx\_range\_1d\_1, [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) const & idx\_range\_1d\_2, ddc::SplineBuilderClosure const & Closure1=ddc::SplineBuilderClosure::HERMITE, ddc::SplineBuilderClosure const & Closure2=ddc::SplineBuilderClosure::HERMITE) <br>_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _._ |
|   | [**DerivReconstructionCoeffs**](#function-derivreconstructioncoeffs-24) (IdxRangeA const & idx\_range\_a, IdxRangeB const & idx\_range\_b, ddc::SplineBuilderClosure const & Closure1=ddc::SplineBuilderClosure::HERMITE, ddc::SplineBuilderClosure const & Closure2=ddc::SplineBuilderClosure::HERMITE) <br>_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _._ |
|   | [**DerivReconstructionCoeffs**](#function-derivreconstructioncoeffs-34) ([**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) const & idx\_range\_1d\_1, [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) const & idx\_range\_1d\_2, std::size\_t const number\_chosen\_cells) <br>_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _. This constructor calculates an approximation of the formula._ |
|   | [**DerivReconstructionCoeffs**](#function-derivreconstructioncoeffs-44) (IdxRangeA const & idx\_range\_a, IdxRangeB const & idx\_range\_b, std::size\_t const number\_chosen\_cells) <br>_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _. This constructor calculates an approximation of the formula._ |
|  double | [**get\_approx\_cross\_deriv**](#function-get_approx_cross_deriv) (DConstField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1), Kokkos::HostSpace, Layout1 &gt; const & derivs\_1, DConstField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2), Kokkos::HostSpace, Layout2 &gt; const & derivs\_2, int const sign\_deriv\_2) const<br>_Get the linear combination of the derivatives values (c)._  |
|  double | [**get\_approx\_deriv**](#function-get_approx_deriv-12) (DConstField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1), Kokkos::HostSpace, Layout1 &gt; const & function\_1, DConstField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2), Kokkos::HostSpace, Layout2 &gt; const & function\_2) const<br>_Get the linear combination of the function values (c)._  |
|  double | [**get\_approx\_deriv**](#function-get_approx_deriv-22) (DConstField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2), Kokkos::HostSpace, Layout2 &gt; const & function\_2, DConstField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1), Kokkos::HostSpace, Layout1 &gt; const & function\_1) const<br>_Get the linear combination of the function values (c). See_ get\_approx\_deriv _._ |
|  double | [**get\_coeff\_deriv\_patch\_1**](#function-get_coeff_deriv_patch_1) () const<br>_Get the coefficient (a) in front of the derivative on the patch 1 of the given_ [_**Interface**_](structInterface.md) _._ |
|  double | [**get\_coeff\_deriv\_patch\_2**](#function-get_coeff_deriv_patch_2) () const<br>_Get the coefficient (b) in front of the derivative on the patch 2 of the given_ [_**Interface**_](structInterface.md) _._ |




## Protected Types

| Type | Name |
| ---: | :--- |
| typedef BSplinesPerp1Type | [**BSplinesPerp1**](#typedef-bsplinesperp1)  <br>_The B-splines in the direction perpendicular to the interface on patch 1._  |
| typedef BSplinesPerp2Type | [**BSplinesPerp2**](#typedef-bsplinesperp2)  <br>_The B-splines in the direction perpendicular to the interface on patch 2._  |
| typedef ddc::knot\_discrete\_dimension\_t&lt; [**BSplinesPerp1**](classDerivReconstructionCoeffs.md#typedef-bsplinesperp1) &gt; | [**GridBreakPt1**](#typedef-gridbreakpt1)  <br>_The grid of break points of the B-splines on patch 1._  |
| typedef ddc::knot\_discrete\_dimension\_t&lt; [**BSplinesPerp2**](classDerivReconstructionCoeffs.md#typedef-bsplinesperp2) &gt; | [**GridBreakPt2**](#typedef-gridbreakpt2)  <br>_The grid of break points of the B-splines on patch 2._  |
| typedef EdgePerpGrid1Type | [**GridEdgePerp1**](#typedef-gridedgeperp1)  <br>_The grid of interpolation points perpendicular to the interface on patch 1._  |
| typedef EdgePerpGrid2Type | [**GridEdgePerp2**](#typedef-gridedgeperp2)  <br>_The grid of interpolation points perpendicular to the interface on patch 2._  |
| typedef Idx&lt; [**GridEdgePerp1**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp1) &gt; | [**IdxEdgePerp1**](#typedef-idxedgeperp1)  <br>_An index of_ [_**GridEdgePerp1**_](classDerivReconstructionCoeffs.md#typedef-gridedgeperp1) _._ |
| typedef Idx&lt; [**GridEdgePerp2**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp2) &gt; | [**IdxEdgePerp2**](#typedef-idxedgeperp2)  <br>_An index of_ [_**GridEdgePerp2**_](classDerivReconstructionCoeffs.md#typedef-gridedgeperp2) _._ |
| typedef IdxRange&lt; [**GridEdgePerp1**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp1) &gt; | [**IdxRange1DPerp\_1**](#typedef-idxrange1dperp_1)  <br>_The index range over the interpolation points perpendicular to the interface on patch 1._  |
| typedef IdxRange&lt; [**GridEdgePerp2**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp2) &gt; | [**IdxRange1DPerp\_2**](#typedef-idxrange1dperp_2)  <br>_The index range over the interpolation points perpendicular to the interface on patch 2._  |




## Protected Attributes

| Type | Name |
| ---: | :--- |
|  double | [**m\_coeff\_deriv\_patch\_1**](#variable-m_coeff_deriv_patch_1)  <br>_The coefficient (a) in front of the derivative on patch 1._  |
|  double | [**m\_coeff\_deriv\_patch\_2**](#variable-m_coeff_deriv_patch_2)  <br>_The coefficient (b) in front of the derivative on patch 2._  |
|  [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) const | [**m\_idx\_range\_perp\_1**](#variable-m_idx_range_perp_1)  <br>_The index range of the interpolation points used on patch 1._  |
|  [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) const | [**m\_idx\_range\_perp\_2**](#variable-m_idx_range_perp_2)  <br>_The index range of the interpolation points used on patch 2._  |
|  const bool | [**m\_is\_cell\_bound\_1\_with\_extra\_interpol\_pt**](#variable-m_is_cell_bound_1_with_extra_interpol_pt)  <br>_True if patch 1 has an additional interpolation point in its boundary cell (Greville closure)._  |
|  const bool | [**m\_is\_cell\_bound\_2\_with\_extra\_interpol\_pt**](#variable-m_is_cell_bound_2_with_extra_interpol_pt)  <br>_True if patch 2 has an additional interpolation point in its boundary cell (Greville closure)._  |
|  host\_t&lt; DField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) &gt; &gt; | [**m\_weights\_patch\_1**](#variable-m_weights_patch_1)  <br>_The weights_ \(\omega_k\) _of the function values on patch 1 in the coefficient (c)._ |
|  host\_t&lt; DFieldMem&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) &gt; &gt; | [**m\_weights\_patch\_1\_alloc**](#variable-m_weights_patch_1_alloc)  <br>_The memory allocated for the weights of the function values on patch 1._  |
|  host\_t&lt; DField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) &gt; &gt; | [**m\_weights\_patch\_2**](#variable-m_weights_patch_2)  <br>_The weights_ \(\omega_k\) _of the function values on patch 2 in the coefficient (c)._ |
|  host\_t&lt; DFieldMem&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) &gt; &gt; | [**m\_weights\_patch\_2\_alloc**](#variable-m_weights_patch_2_alloc)  <br>_The memory allocated for the weights of the function values on patch 2._  |


## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  std::size\_t | [**n**](#variable-n)   = `2`<br>_The position of the n-th index in the arrays used for the recursive calculation._  |
|  std::size\_t | [**n\_minus\_1**](#variable-n_minus_1)   = `1`<br>_The position of the (n-1)-th index in the arrays used for the recursive calculation._  |
|  std::size\_t | [**n\_minus\_2**](#variable-n_minus_2)   = `0`<br>_The position of the (n-2)-th index in the arrays used for the recursive calculation._  |
|  Extremity constexpr | [**s\_extremity\_1**](#variable-s_extremity_1)   = `extremity\_1`<br>_The extremity of patch 1 where the interface is located._  |
|  Extremity constexpr | [**s\_extremity\_2**](#variable-s_extremity_2)   = `extremity\_2`<br>_The extremity of patch 2 where the interface is located._  |














## Protected Functions

| Type | Name |
| ---: | :--- |
|  void | [**check\_additional\_interpolation\_points\_location**](#function-check_additional_interpolation_points_location) (IdxRange&lt; EdgePerpGrid &gt; const & idx\_range\_perp, Extremity const extremity) <br>_Check that the additional interpolation point is located in the boundary cell of the patch opposite the interface. An exception is raised if this is not the case._  |
|  void | [**check\_break\_points\_are\_interpolation\_points**](#function-check_break_points_are_interpolation_points) (bool const is\_cell\_bound\_with\_extra\_interpol\_pt, IdxRange&lt; EdgePerpGrid &gt; const & idx\_range\_perp, Extremity const extremity) const<br>_Check that the break points of the B-splines are interpolation points. An exception is raised if this is not the case._  |
|  void | [**correction\_boundary**](#function-correction_boundary) (int const n\_points, Idx1D const interface\_idx, std::array&lt; double, 3 &gt; & coeff\_deriv\_patch, std::array&lt; double, 3 &gt; & coeff\_deriv\_other\_patch, std::array&lt; host\_t&lt; DField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) &gt; &gt;, 3 &gt; & weights\_patch1, std::array&lt; host\_t&lt; DField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) &gt; &gt;, 3 &gt; & weights\_patch2) <br>_Correction of the local coefficients on the boundary for additional interpolation points in the middle of a cell._  |
|  double | [**get\_approx\_with\_change\_of\_sign**](#function-get_approx_with_change_of_sign) (DConstField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1), Kokkos::HostSpace, Layout1 &gt; const & function\_1, DConstField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2), Kokkos::HostSpace, Layout2 &gt; const & function\_2, int const sign\_function\_2) const<br>_Get the linear combination of the function values (c), changing the sign of the values on patch 2 if required._  |
|  std::tuple&lt; double, double &gt; | [**get\_cell\_lengths**](#function-get_cell_lengths) (Idx1D const & idx) const<br>_Get the lengths of the cell on the right side and the left side of the given index._  |
|  KOKKOS\_INLINE\_FUNCTION Idx&lt; Grid1D &gt; | [**get\_extremity\_idx**](#function-get_extremity_idx) (Extremity const extremity, IdxRange&lt; Grid1D &gt; const & idx\_range) const<br>_Get the index at the given extremity of an index range._  |
|  std::array&lt; double, 3 &gt; | [**get\_gammas**](#function-get_gammas) (double const cell\_length\_left, double const cell\_length\_right) const<br>_Compute the coefficients gammas\_i._  |
|  void | [**recursion**](#function-recursion) (int const n\_cells, Idx1D const interface\_idx, std::array&lt; double, 3 &gt; & coeff\_deriv\_patch, std::array&lt; double, 3 &gt; & coeff\_deriv\_other\_patch, std::array&lt; host\_t&lt; DField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) &gt; &gt;, 3 &gt; & weights\_patch1, std::array&lt; host\_t&lt; DField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) &gt; &gt;, 3 &gt; & weights\_patch2) <br>_Recursion to compute the coefficients a, b and the weights omega._  |
|  void | [**set\_coefficients\_non\_uniform\_case**](#function-set_coefficients_non_uniform_case) () <br>_Compute the coefficients a, b and the weights omega applying the recursive formula._  |
|  void | [**set\_coefficients\_uniform\_per\_patch\_case**](#function-set_coefficients_uniform_per_patch_case) () <br>_Compute the coefficients a, b and the weights omega applying the explicit formula for the uniform case._  |


## Protected Static Functions

| Type | Name |
| ---: | :--- |
|  std::array&lt; DField&lt; IdxRangeType, Kokkos::HostSpace &gt;, 3 &gt; | [**build\_array**](#function-build_array) (int input\_position, DField&lt; IdxRangeType, Kokkos::HostSpace &gt; persistent\_memory, DField&lt; IdxRangeType, Kokkos::HostSpace &gt; buffer\_1, DField&lt; IdxRangeType, Kokkos::HostSpace &gt; buffer\_2) <br>_Build an array with the persistent\_memory at the requested input\_position._  |
|  void | [**cycle\_buffers**](#function-cycle_buffers-12) (std::array&lt; double, 3 &gt; & buffer) <br>_Cycle the roles of the doubles in the array of values._  |
|  void | [**cycle\_buffers**](#function-cycle_buffers-22) (std::array&lt; DField&lt; IdxRangeType, MemorySpace &gt;, 3 &gt; & buffer) <br>_Cycle the memory pointers in the array._  |
|  double | [**get\_alpha**](#function-get_alpha) (double const cell\_length\_left, double const cell\_length\_right) <br>_Compute the coefficient alpha\_i._  |
|  double | [**get\_beta**](#function-get_beta) (double const cell\_length\_left, double const cell\_length\_right) <br>_Compute the coefficient beta\_i._  |
|  void | [**set\_local\_uniform\_weights**](#function-set_local_uniform_weights) (host\_t&lt; DField&lt; IdxRange&lt; EdgePerpGrid &gt; &gt; &gt; weights, double factor, double un\_other\_patch, Idx&lt; EdgePerpGrid &gt; idx\_interface) <br>_Apply the explicit formula for the uniform case to weights on one side of the interface._  |


## Detailed Description


For a given [**Interface**](structInterface.md), this operator computes the coefficients a, b and c of the following relation between the interface derivative and the derivatives at the other boundaries of the two patches: \(s'(X_I) = c + a s'(X_{I+1}) + b s'(X_{I+1})\),


with
* s'(X\_I) the derivative at the [**Interface**](structInterface.md) of the two local splines on the patches.
* s'(X\_{I+1}) the derivative at the next (or right) [**Interface**](structInterface.md) of the local splines on patch 1.
* s'(X\_{I-1}) the derivative at the previous (or left) [**Interface**](structInterface.md) of the local splines on patch 2.
* a,b scalars depending on the meshes of the two patches.
* c a linear combination of the function values on the two patches. The weights in the linear combination depend on the meshes of the two patches. It can be written as a sum: \(c = \sum_{k = -N2}^{N1} \omega_k f_k\), with \(\omega_k\) the weights, \(f_k\) the function values,




By default, \(N1\) the number of cells in patch 1, and \(N2\) the number of cells in patch 2.


But a constructor is also provided to choose $N\_{trunc}=N1=N2$ in order to calculate an approximation of the derivative: \(s'(X_I) \simeq c_{trunc} = \sum_{k = -N_{trunc}}^{N_{trunc}} \omega_k f_k\).


Scheme of the two patches: X\_{I-1} X\_I X\_{I+1} 
 \| \| \| \| 1 \| 2 \| \|\_\_\_\_\_\_\_\|\_\_\_\_\_\_\_\|
* ← → +
* → ← - Orientation of an equivalent global mesh : →




All the formulae and more details are given in the README.md.




**Template parameters:**


* `InterfaceType` The interface between two patches where we want to compute the derivatives. 




**Warning:**

The applied method only works for interpolation points located on the break points. In the case where "ddc::SplineBuilderClosure::GREVILLE" is specified, additional interpolation points are also placed in the first or last cell of the patch. Please be sure to initialise the discrete space of your Grid on the break points (especially in the non-uniform case). 





    
## Public Functions Documentation




### function DerivReconstructionCoeffs [1/4]

_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _._
```C++
inline DerivReconstructionCoeffs::DerivReconstructionCoeffs (
    IdxRange1DPerp_1 const & idx_range_1d_1,
    IdxRange1DPerp_2 const & idx_range_1d_2,
    ddc::SplineBuilderClosure const & Closure1=ddc::SplineBuilderClosure::HERMITE,
    ddc::SplineBuilderClosure const & Closure2=ddc::SplineBuilderClosure::HERMITE
) 
```



 It computes the coefficients a and b, and the weights \(\omega_k\) and stores the values in the class. If you want to compute the derivatives with an exact formula, please provide the index ranges on the full domain of the patches. If you want to use an approximation (e.g. use only 5 cells), please provide the index ranges on a reduced number of cells (e.g. 5 cells).


If the interpolation points are uniform, it computes the coefficients with an explicit formula. Otherwise, it uses the recursive formula (see README.md).




**Parameters:**


* `idx_range_1d_1` 1D index range perpendicular to the [**Interface**](structInterface.md), on the patch 1. 
* `idx_range_1d_2` 1D index range perpendicular to the [**Interface**](structInterface.md), on the patch 2. 
* `Closure1` The spline closure type on the opposite edge of the interface on the patch 1. By default, the value is set to ddc::SplineBuilderClosure::HERMITE. If ddc::SplineBuilderClosure::GREVILLE is given, a treatment will be applied to consider the additional interpolation point. Giving ddc::SplineBuilderClosure::PERIODIC does not make sense. 
* `Closure2` The spline closure type on the opposite edge of the interface on the patch 2. By default, the value is set to ddc::SplineBuilderClosure::HERMITE. If ddc::SplineBuilderClosure::GREVILLE is given, a treatment will be applied to consider the additional interpolation point. Giving ddc::SplineBuilderClosure::PERIODIC does not make sense. 




        

<hr>



### function DerivReconstructionCoeffs [2/4]

_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _._
```C++
template<class IdxRangeA, class IdxRangeB>
inline DerivReconstructionCoeffs::DerivReconstructionCoeffs (
    IdxRangeA const & idx_range_a,
    IdxRangeB const & idx_range_b,
    ddc::SplineBuilderClosure const & Closure1=ddc::SplineBuilderClosure::HERMITE,
    ddc::SplineBuilderClosure const & Closure2=ddc::SplineBuilderClosure::HERMITE
) 
```





**Parameters:**


* `idx_range_a` Index range on one patch. 
* `idx_range_b` Index range on the other patch. 
* `Closure1` The spline closure type on the opposite edge of the interface on the patch 1. By default, the value is set to ddc::SplineBuilderClosure::HERMITE. 
* `Closure2` The spline closure type on the opposite edge of the interface on the patch 2. By default, the value is set to ddc::SplineBuilderClosure::HERMITE. 




        

<hr>



### function DerivReconstructionCoeffs [3/4]

_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _. This constructor calculates an approximation of the formula._
```C++
inline DerivReconstructionCoeffs::DerivReconstructionCoeffs (
    IdxRange1DPerp_1 const & idx_range_1d_1,
    IdxRange1DPerp_2 const & idx_range_1d_2,
    std::size_t const number_chosen_cells
) 
```





**Parameters:**


* `idx_range_1d_1` 1D index range perpendicular to the [**Interface**](structInterface.md), on the patch 1. 
* `idx_range_1d_2` 1D index range perpendicular to the [**Interface**](structInterface.md), on the patch 2. 
* `number_chosen_cells` The number of cells we choose for the computation of the interface derivatives. Same number for Patch1 and Patch2. 



**Warning:**

If the mesh has additional interpolation points, please be sure to not take these points. 





        

<hr>



### function DerivReconstructionCoeffs [4/4]

_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _. This constructor calculates an approximation of the formula._
```C++
template<class IdxRangeA, class IdxRangeB>
inline DerivReconstructionCoeffs::DerivReconstructionCoeffs (
    IdxRangeA const & idx_range_a,
    IdxRangeB const & idx_range_b,
    std::size_t const number_chosen_cells
) 
```





**Parameters:**


* `idx_range_a` Index range on one patch. 
* `idx_range_b` Index range on the other patch. 
* `number_chosen_cells` The number of cells we choose for the computation of the interface derivatives. Same number for Patch1 and Patch2. 



**Warning:**

If the mesh has additional interpolation points, please be sure to not take these points. 





        

<hr>



### function get\_approx\_cross\_deriv 

_Get the linear combination of the derivatives values (c)._ 
```C++
template<class Layout1, class Layout2>
inline double DerivReconstructionCoeffs::get_approx_cross_deriv (
    DConstField< IdxRange1DPerp_1 , Kokkos::HostSpace, Layout1 > const & derivs_1,
    DConstField< IdxRange1DPerp_2 , Kokkos::HostSpace, Layout2 > const & derivs_2,
    int const sign_deriv_2
) const
```



This linear combination corresponds to an approximation of the interface cross-derivative.


It is more appropriate than get\_approx\_deriv to compute the cross-derivatives. The computation is the same. However, in addition, the sign of the derivatives on the patch 2 is changed in the case where the orientations of the two patches do not match. E.g. dx(f(-x)) = - dxf(-x)




**Parameters:**


* `derivs_1` Derivatives at the interpolation points on patch 1. 
* `derivs_2` Derivatives at the interpolation points on patch 2. 
* `sign_deriv_2` The sign (+1 or -1) applied to the derivatives on patch 2 to account for the relative orientation of the two patches. 



**Returns:**

the linear combination of the function values (c). 





        

<hr>



### function get\_approx\_deriv [1/2]

_Get the linear combination of the function values (c)._ 
```C++
template<class Layout1, class Layout2>
inline double DerivReconstructionCoeffs::get_approx_deriv (
    DConstField< IdxRange1DPerp_1 , Kokkos::HostSpace, Layout1 > const & function_1,
    DConstField< IdxRange1DPerp_2 , Kokkos::HostSpace, Layout2 > const & function_2
) const
```



This linear combination corresponds to an approximation of the interface derivative.




**Parameters:**


* `function_1` Function values at the interpolation points on patch 1. 
* `function_2` Function values at the interpolation points on patch 2. 



**Returns:**

the linear combination of the function values (c).


Remark: for the approximation case, the cells were selected by the given index ranges in the instantiation. Here, function\_1 and function\_2 can be given on the whole domain and the operator will select the correct values. But, it could be more optimised to directly give a slice on the selected cells. E.g. get\_approx\_deriv(function\_1[idx\_interface\_1][selected\_cells\_idx\_range\_1], 
                 function\_2[idx\_interface\_2][selected\_cells\_idx\_range\_2]); and get\_approx\_deriv(function\_1[idx\_interface\_1], 
                 function\_2[idx\_interface\_2]); will return the same value. 


        

<hr>



### function get\_approx\_deriv [2/2]

_Get the linear combination of the function values (c). See_ get\_approx\_deriv _._
```C++
template<class Layout1, class Layout2>
inline double DerivReconstructionCoeffs::get_approx_deriv (
    DConstField< IdxRange1DPerp_2 , Kokkos::HostSpace, Layout2 > const & function_2,
    DConstField< IdxRange1DPerp_1 , Kokkos::HostSpace, Layout1 > const & function_1
) const
```





**Parameters:**


* `function_1` Function values at the interpolation points on patch 1. 
* `function_2` Function values at the interpolation points on patch 2. 



**Returns:**

the linear combination of the function values (c). 





        

<hr>



### function get\_coeff\_deriv\_patch\_1 

_Get the coefficient (a) in front of the derivative on the patch 1 of the given_ [_**Interface**_](structInterface.md) _._
```C++
inline double DerivReconstructionCoeffs::get_coeff_deriv_patch_1 () const
```





**Returns:**

The value of the coefficient a. 





        

<hr>



### function get\_coeff\_deriv\_patch\_2 

_Get the coefficient (b) in front of the derivative on the patch 2 of the given_ [_**Interface**_](structInterface.md) _._
```C++
inline double DerivReconstructionCoeffs::get_coeff_deriv_patch_2 () const
```





**Returns:**

The value of the coefficient b. 





        

<hr>
## Protected Types Documentation




### typedef BSplinesPerp1 

_The B-splines in the direction perpendicular to the interface on patch 1._ 
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::BSplinesPerp1 =  BSplinesPerp1Type;
```




<hr>



### typedef BSplinesPerp2 

_The B-splines in the direction perpendicular to the interface on patch 2._ 
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::BSplinesPerp2 =  BSplinesPerp2Type;
```




<hr>



### typedef GridBreakPt1 

_The grid of break points of the B-splines on patch 1._ 
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::GridBreakPt1 =  ddc::knot_discrete_dimension_t<BSplinesPerp1>;
```




<hr>



### typedef GridBreakPt2 

_The grid of break points of the B-splines on patch 2._ 
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::GridBreakPt2 =  ddc::knot_discrete_dimension_t<BSplinesPerp2>;
```




<hr>



### typedef GridEdgePerp1 

_The grid of interpolation points perpendicular to the interface on patch 1._ 
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::GridEdgePerp1 =  EdgePerpGrid1Type;
```




<hr>



### typedef GridEdgePerp2 

_The grid of interpolation points perpendicular to the interface on patch 2._ 
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::GridEdgePerp2 =  EdgePerpGrid2Type;
```




<hr>



### typedef IdxEdgePerp1 

_An index of_ [_**GridEdgePerp1**_](classDerivReconstructionCoeffs.md#typedef-gridedgeperp1) _._
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::IdxEdgePerp1 =  Idx<GridEdgePerp1>;
```




<hr>



### typedef IdxEdgePerp2 

_An index of_ [_**GridEdgePerp2**_](classDerivReconstructionCoeffs.md#typedef-gridedgeperp2) _._
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::IdxEdgePerp2 =  Idx<GridEdgePerp2>;
```




<hr>



### typedef IdxRange1DPerp\_1 

_The index range over the interpolation points perpendicular to the interface on patch 1._ 
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::IdxRange1DPerp_1 =  IdxRange<GridEdgePerp1>;
```




<hr>



### typedef IdxRange1DPerp\_2 

_The index range over the interpolation points perpendicular to the interface on patch 2._ 
```C++
using DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::IdxRange1DPerp_2 =  IdxRange<GridEdgePerp2>;
```




<hr>
## Protected Attributes Documentation




### variable m\_coeff\_deriv\_patch\_1 

_The coefficient (a) in front of the derivative on patch 1._ 
```C++
double DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_coeff_deriv_patch_1;
```




<hr>



### variable m\_coeff\_deriv\_patch\_2 

_The coefficient (b) in front of the derivative on patch 2._ 
```C++
double DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_coeff_deriv_patch_2;
```




<hr>



### variable m\_idx\_range\_perp\_1 

_The index range of the interpolation points used on patch 1._ 
```C++
IdxRange1DPerp_1 const DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_idx_range_perp_1;
```




<hr>



### variable m\_idx\_range\_perp\_2 

_The index range of the interpolation points used on patch 2._ 
```C++
IdxRange1DPerp_2 const DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_idx_range_perp_2;
```




<hr>



### variable m\_is\_cell\_bound\_1\_with\_extra\_interpol\_pt 

_True if patch 1 has an additional interpolation point in its boundary cell (Greville closure)._ 
```C++
const bool DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_is_cell_bound_1_with_extra_interpol_pt;
```




<hr>



### variable m\_is\_cell\_bound\_2\_with\_extra\_interpol\_pt 

_True if patch 2 has an additional interpolation point in its boundary cell (Greville closure)._ 
```C++
const bool DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_is_cell_bound_2_with_extra_interpol_pt;
```




<hr>



### variable m\_weights\_patch\_1 

_The weights_ \(\omega_k\) _of the function values on patch 1 in the coefficient (c)._
```C++
host_t<DField<IdxRange1DPerp_1> > DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_weights_patch_1;
```




<hr>



### variable m\_weights\_patch\_1\_alloc 

_The memory allocated for the weights of the function values on patch 1._ 
```C++
host_t<DFieldMem<IdxRange1DPerp_1> > DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_weights_patch_1_alloc;
```




<hr>



### variable m\_weights\_patch\_2 

_The weights_ \(\omega_k\) _of the function values on patch 2 in the coefficient (c)._
```C++
host_t<DField<IdxRange1DPerp_2> > DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_weights_patch_2;
```




<hr>



### variable m\_weights\_patch\_2\_alloc 

_The memory allocated for the weights of the function values on patch 2._ 
```C++
host_t<DFieldMem<IdxRange1DPerp_2> > DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::m_weights_patch_2_alloc;
```




<hr>
## Protected Static Attributes Documentation




### variable n 

_The position of the n-th index in the arrays used for the recursive calculation._ 
```C++
std::size_t DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::n;
```




<hr>



### variable n\_minus\_1 

_The position of the (n-1)-th index in the arrays used for the recursive calculation._ 
```C++
std::size_t DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::n_minus_1;
```




<hr>



### variable n\_minus\_2 

_The position of the (n-2)-th index in the arrays used for the recursive calculation._ 
```C++
std::size_t DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::n_minus_2;
```




<hr>



### variable s\_extremity\_1 

_The extremity of patch 1 where the interface is located._ 
```C++
Extremity constexpr DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::s_extremity_1;
```




<hr>



### variable s\_extremity\_2 

_The extremity of patch 2 where the interface is located._ 
```C++
Extremity constexpr DerivReconstructionCoeffs< EdgePerpGrid1Type, EdgePerpGrid2Type, BSplinesPerp1Type, BSplinesPerp2Type, extremity_1, extremity_2 >::s_extremity_2;
```




<hr>
## Protected Functions Documentation




### function check\_additional\_interpolation\_points\_location 

_Check that the additional interpolation point is located in the boundary cell of the patch opposite the interface. An exception is raised if this is not the case._ 
```C++
template<typename BSplinesPerp, typename GridBreakPt, typename EdgePerpGrid>
inline void DerivReconstructionCoeffs::check_additional_interpolation_points_location (
    IdxRange< EdgePerpGrid > const & idx_range_perp,
    Extremity const extremity
) 
```





**Parameters:**


* `idx_range_perp` The index range of the interpolation points used on the patch. 
* `extremity` The extremity of the patch on which the interface is located. 




        

<hr>



### function check\_break\_points\_are\_interpolation\_points 

_Check that the break points of the B-splines are interpolation points. An exception is raised if this is not the case._ 
```C++
template<typename BSplinesPerp, typename GridBreakPt, typename EdgePerpGrid>
inline void DerivReconstructionCoeffs::check_break_points_are_interpolation_points (
    bool const is_cell_bound_with_extra_interpol_pt,
    IdxRange< EdgePerpGrid > const & idx_range_perp,
    Extremity const extremity
) const
```





**Parameters:**


* `is_cell_bound_with_extra_interpol_pt` True if the patch has an additional interpolation point in its boundary cell. 
* `idx_range_perp` The index range of the interpolation points used on the patch. 
* `extremity` The extremity of the patch on which the interface is located. 




        

<hr>



### function correction\_boundary 

_Correction of the local coefficients on the boundary for additional interpolation points in the middle of a cell._ 
```C++
template<bool is_on_patch1, class Idx1D>
inline void DerivReconstructionCoeffs::correction_boundary (
    int const n_points,
    Idx1D const interface_idx,
    std::array< double, 3 > & coeff_deriv_patch,
    std::array< double, 3 > & coeff_deriv_other_patch,
    std::array< host_t< DField< IdxRange1DPerp_1 > >, 3 > & weights_patch1,
    std::array< host_t< DField< IdxRange1DPerp_2 > >, 3 > & weights_patch2
) 
```





**Parameters:**


* `n_points` Number of interpolation points on the patch where the correction is made. 
* `interface_idx` Index at the interface.

The following parameters are array with three elements corresponding to the values at the step n-1, n and n+1 in the recursion. 

**Parameters:**


* `coeff_deriv_patch` Coefficient in front of the derivative of the patch. E.g. Patch1: coeff\_deriv\_patch = b \| Patch2: coeff\_deriv\_patch = a 
* `coeff_deriv_other_patch` Coefficient in front of the derivative of the other patch. E.g. Patch1: coeff\_deriv\_other\_patch = a \| Patch2: coeff\_deriv\_other\_patch = b 
* `weights_patch1` Weights in front of the function values on patch 1. 
* `weights_patch2` Weights in front of the function values on patch 2. 




        

<hr>



### function get\_approx\_with\_change\_of\_sign 

_Get the linear combination of the function values (c), changing the sign of the values on patch 2 if required._ 
```C++
template<class Layout1, class Layout2>
inline double DerivReconstructionCoeffs::get_approx_with_change_of_sign (
    DConstField< IdxRange1DPerp_1 , Kokkos::HostSpace, Layout1 > const & function_1,
    DConstField< IdxRange1DPerp_2 , Kokkos::HostSpace, Layout2 > const & function_2,
    int const sign_function_2
) const
```





**Parameters:**


* `function_1` Function values at the interpolation points on patch 1. 
* `function_2` Function values at the interpolation points on patch 2. 
* `sign_function_2` The sign (+1 or -1) applied to the function values on patch 2. 



**Returns:**

the linear combination of the function values (c). 





        

<hr>



### function get\_cell\_lengths 

_Get the lengths of the cell on the right side and the left side of the given index._ 
```C++
template<bool is_on_patch1, class Idx1D>
inline std::tuple< double, double > DerivReconstructionCoeffs::get_cell_lengths (
    Idx1D const & idx
) const
```




<hr>



### function get\_extremity\_idx 

_Get the index at the given extremity of an index range._ 
```C++
template<class Grid1D>
inline KOKKOS_INLINE_FUNCTION Idx< Grid1D > DerivReconstructionCoeffs::get_extremity_idx (
    Extremity const extremity,
    IdxRange< Grid1D > const & idx_range
) const
```





**Parameters:**


* `extremity` The extremity (FRONT or BACK). 
* `idx_range` The index range. 



**Returns:**

The first index if extremity is FRONT, the last index otherwise. 





        

<hr>



### function get\_gammas 

_Compute the coefficients gammas\_i._ 
```C++
inline std::array< double, 3 > DerivReconstructionCoeffs::get_gammas (
    double const cell_length_left,
    double const cell_length_right
) const
```




<hr>



### function recursion 

_Recursion to compute the coefficients a, b and the weights omega._ 
```C++
template<bool is_on_patch1, class Idx1D>
inline void DerivReconstructionCoeffs::recursion (
    int const n_cells,
    Idx1D const interface_idx,
    std::array< double, 3 > & coeff_deriv_patch,
    std::array< double, 3 > & coeff_deriv_other_patch,
    std::array< host_t< DField< IdxRange1DPerp_1 > >, 3 > & weights_patch1,
    std::array< host_t< DField< IdxRange1DPerp_2 > >, 3 > & weights_patch2
) 
```





**Parameters:**


* `n_cells` Number of the cells of the patch where the recursion is made. 
* `interface_idx` Index at the interface.

The following parameters are array with three elements corresponding to the values at the step n-1, n and n+1 in the recursion. 

**Parameters:**


* `coeff_deriv_patch` Coefficient in front of the derivative of the patch. E.g. Patch1: coeff\_deriv\_patch = b \| Patch2: coeff\_deriv\_patch = a 
* `coeff_deriv_other_patch` Coefficient in front of the derivative of the other patch. E.g. Patch1: coeff\_deriv\_other\_patch = a \| Patch2: coeff\_deriv\_other\_patch = b 
* `weights_patch1` Weights in front of the function values on patch 1. 
* `weights_patch2` Weights in front of the function values on patch 2. 




        

<hr>



### function set\_coefficients\_non\_uniform\_case 

_Compute the coefficients a, b and the weights omega applying the recursive formula._ 
```C++
inline void DerivReconstructionCoeffs::set_coefficients_non_uniform_case () 
```




<hr>



### function set\_coefficients\_uniform\_per\_patch\_case 

_Compute the coefficients a, b and the weights omega applying the explicit formula for the uniform case._ 
```C++
inline void DerivReconstructionCoeffs::set_coefficients_uniform_per_patch_case () 
```




<hr>
## Protected Static Functions Documentation




### function build\_array 

_Build an array with the persistent\_memory at the requested input\_position._ 
```C++
template<class IdxRangeType>
static inline std::array< DField< IdxRangeType, Kokkos::HostSpace >, 3 > DerivReconstructionCoeffs::build_array (
    int input_position,
    DField< IdxRangeType, Kokkos::HostSpace > persistent_memory,
    DField< IdxRangeType, Kokkos::HostSpace > buffer_1,
    DField< IdxRangeType, Kokkos::HostSpace > buffer_2
) 
```





**Parameters:**


* `input_position` The position of the persistent\_memory in the output array. 
* `persistent_memory` The memory that is saved in the class. 
* `buffer_1` A temporary empty buffer to be saved in the array. 
* `buffer_2` A temporary empty buffer to be saved in the array. 




        

<hr>



### function cycle\_buffers [1/2]

_Cycle the roles of the doubles in the array of values._ 
```C++
static inline void DerivReconstructionCoeffs::cycle_buffers (
    std::array< double, 3 > & buffer
) 
```



Cycle the doubles in the array to change the indices ready for the next iteration.




**Parameters:**


* `buffer` The array of values for the indices n-2, n-1 and n. 




        

<hr>



### function cycle\_buffers [2/2]

_Cycle the memory pointers in the array._ 
```C++
template<class IdxRangeType, class MemorySpace>
static inline void DerivReconstructionCoeffs::cycle_buffers (
    std::array< DField< IdxRangeType, MemorySpace >, 3 > & buffer
) 
```



Cycle the memory pointers in the array to change the indices ready for the next iteration.




**Parameters:**


* `buffer` The array of fields for the indices n-2, n-1 and n. 




        

<hr>



### function get\_alpha 

_Compute the coefficient alpha\_i._ 
```C++
static inline double DerivReconstructionCoeffs::get_alpha (
    double const cell_length_left,
    double const cell_length_right
) 
```




<hr>



### function get\_beta 

_Compute the coefficient beta\_i._ 
```C++
static inline double DerivReconstructionCoeffs::get_beta (
    double const cell_length_left,
    double const cell_length_right
) 
```




<hr>



### function set\_local\_uniform\_weights 

_Apply the explicit formula for the uniform case to weights on one side of the interface._ 
```C++
template<class EdgePerpGrid>
static inline void DerivReconstructionCoeffs::set_local_uniform_weights (
    host_t< DField< IdxRange< EdgePerpGrid > > > weights,
    double factor,
    double un_other_patch,
    Idx< EdgePerpGrid > idx_interface
) 
```





**Parameters:**


* `weights` The weights on the patch which should be filled (except at the interface). 
* `factor` The multiplicative factor associated with this patch. 
* `un_other_patch` The value \(u_N\) calculated for the other patch. 
* `idx_interface` The index of the interface point on this patch. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/interface_derivatives/deriv_reconstruction_coeffs.hpp`

