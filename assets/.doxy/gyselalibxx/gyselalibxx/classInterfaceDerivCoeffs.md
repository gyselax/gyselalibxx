

# Class InterfaceDerivCoeffs

**template &lt;class InterfaceType&gt;**



[**ClassList**](annotated.md) **>** [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md)



_Compute the coefficients a, b and c of the interface derivative reconstruction method for a given_ [_**Interface**_](structInterface.md) _between two patches._[More...](#detailed-description)

* `#include <interface_derivative_coefficients.hpp>`



Inherits the following classes: [DerivReconstructionCoeffs](classDerivReconstructionCoeffs.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef InterfaceType | [**associated\_interface**](#typedef-associated_interface)  <br>[_**Interface**_](structInterface.md) _between the two involved patches._ |








































## Public Functions

| Type | Name |
| ---: | :--- |
|  double | [**get\_approx\_cross\_deriv**](#function-get_approx_cross_deriv-12) (DConstField&lt; IdxRange1DPerp\_1, Kokkos::HostSpace, Layout1 &gt; const & derivs\_1, DConstField&lt; IdxRange1DPerp\_2, Kokkos::HostSpace, Layout2 &gt; const & derivs\_2) const<br>_Get the linear combination of the derivatives values (c)._  |
|  double | [**get\_approx\_cross\_deriv**](#function-get_approx_cross_deriv-22) (DConstField&lt; IdxRange1DPerp\_2, Kokkos::HostSpace, Layout2 &gt; const & derivs\_2, DConstField&lt; IdxRange1DPerp\_1, Kokkos::HostSpace, Layout1 &gt; const & derivs\_1) const<br>_Get the linear combination of the derivatives values (c). See_ get\_approx\_cross\_deriv _._ |
|  double | [**get\_coeff\_deriv\_on\_patch**](#function-get_coeff_deriv_on_patch) () const<br>_Get the coefficient (a) or (b) in front of the derivative on the given patch of the given_ [_**Interface**_](structInterface.md) _._ |


## Public Functions inherited from DerivReconstructionCoeffs

See [DerivReconstructionCoeffs](classDerivReconstructionCoeffs.md)

| Type | Name |
| ---: | :--- |
|   | [**DerivReconstructionCoeffs**](classDerivReconstructionCoeffs.md#function-derivreconstructioncoeffs-14) ([**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) const & idx\_range\_1d\_1, [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) const & idx\_range\_1d\_2, ddc::SplineBuilderClosure const & Closure1=ddc::SplineBuilderClosure::HERMITE, ddc::SplineBuilderClosure const & Closure2=ddc::SplineBuilderClosure::HERMITE) <br>_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _._ |
|   | [**DerivReconstructionCoeffs**](classDerivReconstructionCoeffs.md#function-derivreconstructioncoeffs-24) (IdxRangeA const & idx\_range\_a, IdxRangeB const & idx\_range\_b, ddc::SplineBuilderClosure const & Closure1=ddc::SplineBuilderClosure::HERMITE, ddc::SplineBuilderClosure const & Closure2=ddc::SplineBuilderClosure::HERMITE) <br>_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _._ |
|   | [**DerivReconstructionCoeffs**](classDerivReconstructionCoeffs.md#function-derivreconstructioncoeffs-34) ([**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) const & idx\_range\_1d\_1, [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) const & idx\_range\_1d\_2, std::size\_t const number\_chosen\_cells) <br>_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _. This constructor calculates an approximation of the formula._ |
|   | [**DerivReconstructionCoeffs**](classDerivReconstructionCoeffs.md#function-derivreconstructioncoeffs-44) (IdxRangeA const & idx\_range\_a, IdxRangeB const & idx\_range\_b, std::size\_t const number\_chosen\_cells) <br>_Instantiate_ [_**DerivReconstructionCoeffs**_](classDerivReconstructionCoeffs.md) _. See_InterfaceDerivCoeffsInstantiator _. This constructor calculates an approximation of the formula._ |
|  double | [**get\_approx\_cross\_deriv**](classDerivReconstructionCoeffs.md#function-get_approx_cross_deriv) (DConstField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1), Kokkos::HostSpace, Layout1 &gt; const & derivs\_1, DConstField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2), Kokkos::HostSpace, Layout2 &gt; const & derivs\_2, int const sign\_deriv\_2) const<br>_Get the linear combination of the derivatives values (c)._  |
|  double | [**get\_approx\_deriv**](classDerivReconstructionCoeffs.md#function-get_approx_deriv-12) (DConstField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1), Kokkos::HostSpace, Layout1 &gt; const & function\_1, DConstField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2), Kokkos::HostSpace, Layout2 &gt; const & function\_2) const<br>_Get the linear combination of the function values (c)._  |
|  double | [**get\_approx\_deriv**](classDerivReconstructionCoeffs.md#function-get_approx_deriv-22) (DConstField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2), Kokkos::HostSpace, Layout2 &gt; const & function\_2, DConstField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1), Kokkos::HostSpace, Layout1 &gt; const & function\_1) const<br>_Get the linear combination of the function values (c). See_ get\_approx\_deriv _._ |
|  double | [**get\_coeff\_deriv\_patch\_1**](classDerivReconstructionCoeffs.md#function-get_coeff_deriv_patch_1) () const<br>_Get the coefficient (a) in front of the derivative on the patch 1 of the given_ [_**Interface**_](structInterface.md) _._ |
|  double | [**get\_coeff\_deriv\_patch\_2**](classDerivReconstructionCoeffs.md#function-get_coeff_deriv_patch_2) () const<br>_Get the coefficient (b) in front of the derivative on the patch 2 of the given_ [_**Interface**_](structInterface.md) _._ |








## Protected Types inherited from DerivReconstructionCoeffs

See [DerivReconstructionCoeffs](classDerivReconstructionCoeffs.md)

| Type | Name |
| ---: | :--- |
| typedef BSplinesPerp1Type | [**BSplinesPerp1**](classDerivReconstructionCoeffs.md#typedef-bsplinesperp1)  <br>_The B-splines in the direction perpendicular to the interface on patch 1._  |
| typedef BSplinesPerp2Type | [**BSplinesPerp2**](classDerivReconstructionCoeffs.md#typedef-bsplinesperp2)  <br>_The B-splines in the direction perpendicular to the interface on patch 2._  |
| typedef ddc::knot\_discrete\_dimension\_t&lt; [**BSplinesPerp1**](classDerivReconstructionCoeffs.md#typedef-bsplinesperp1) &gt; | [**GridBreakPt1**](classDerivReconstructionCoeffs.md#typedef-gridbreakpt1)  <br>_The grid of break points of the B-splines on patch 1._  |
| typedef ddc::knot\_discrete\_dimension\_t&lt; [**BSplinesPerp2**](classDerivReconstructionCoeffs.md#typedef-bsplinesperp2) &gt; | [**GridBreakPt2**](classDerivReconstructionCoeffs.md#typedef-gridbreakpt2)  <br>_The grid of break points of the B-splines on patch 2._  |
| typedef EdgePerpGrid1Type | [**GridEdgePerp1**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp1)  <br>_The grid of interpolation points perpendicular to the interface on patch 1._  |
| typedef EdgePerpGrid2Type | [**GridEdgePerp2**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp2)  <br>_The grid of interpolation points perpendicular to the interface on patch 2._  |
| typedef Idx&lt; [**GridEdgePerp1**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp1) &gt; | [**IdxEdgePerp1**](classDerivReconstructionCoeffs.md#typedef-idxedgeperp1)  <br>_An index of_ [_**GridEdgePerp1**_](classDerivReconstructionCoeffs.md#typedef-gridedgeperp1) _._ |
| typedef Idx&lt; [**GridEdgePerp2**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp2) &gt; | [**IdxEdgePerp2**](classDerivReconstructionCoeffs.md#typedef-idxedgeperp2)  <br>_An index of_ [_**GridEdgePerp2**_](classDerivReconstructionCoeffs.md#typedef-gridedgeperp2) _._ |
| typedef IdxRange&lt; [**GridEdgePerp1**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp1) &gt; | [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1)  <br>_The index range over the interpolation points perpendicular to the interface on patch 1._  |
| typedef IdxRange&lt; [**GridEdgePerp2**](classDerivReconstructionCoeffs.md#typedef-gridedgeperp2) &gt; | [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2)  <br>_The index range over the interpolation points perpendicular to the interface on patch 2._  |








## Protected Attributes inherited from DerivReconstructionCoeffs

See [DerivReconstructionCoeffs](classDerivReconstructionCoeffs.md)

| Type | Name |
| ---: | :--- |
|  double | [**m\_coeff\_deriv\_patch\_1**](classDerivReconstructionCoeffs.md#variable-m_coeff_deriv_patch_1)  <br>_The coefficient (a) in front of the derivative on patch 1._  |
|  double | [**m\_coeff\_deriv\_patch\_2**](classDerivReconstructionCoeffs.md#variable-m_coeff_deriv_patch_2)  <br>_The coefficient (b) in front of the derivative on patch 2._  |
|  [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) const | [**m\_idx\_range\_perp\_1**](classDerivReconstructionCoeffs.md#variable-m_idx_range_perp_1)  <br>_The index range of the interpolation points used on patch 1._  |
|  [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) const | [**m\_idx\_range\_perp\_2**](classDerivReconstructionCoeffs.md#variable-m_idx_range_perp_2)  <br>_The index range of the interpolation points used on patch 2._  |
|  const bool | [**m\_is\_cell\_bound\_1\_with\_extra\_interpol\_pt**](classDerivReconstructionCoeffs.md#variable-m_is_cell_bound_1_with_extra_interpol_pt)  <br>_True if patch 1 has an additional interpolation point in its boundary cell (Greville closure)._  |
|  const bool | [**m\_is\_cell\_bound\_2\_with\_extra\_interpol\_pt**](classDerivReconstructionCoeffs.md#variable-m_is_cell_bound_2_with_extra_interpol_pt)  <br>_True if patch 2 has an additional interpolation point in its boundary cell (Greville closure)._  |
|  host\_t&lt; DField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) &gt; &gt; | [**m\_weights\_patch\_1**](classDerivReconstructionCoeffs.md#variable-m_weights_patch_1)  <br>_The weights_ \(\omega_k\) _of the function values on patch 1 in the coefficient (c)._ |
|  host\_t&lt; DFieldMem&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) &gt; &gt; | [**m\_weights\_patch\_1\_alloc**](classDerivReconstructionCoeffs.md#variable-m_weights_patch_1_alloc)  <br>_The memory allocated for the weights of the function values on patch 1._  |
|  host\_t&lt; DField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) &gt; &gt; | [**m\_weights\_patch\_2**](classDerivReconstructionCoeffs.md#variable-m_weights_patch_2)  <br>_The weights_ \(\omega_k\) _of the function values on patch 2 in the coefficient (c)._ |
|  host\_t&lt; DFieldMem&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) &gt; &gt; | [**m\_weights\_patch\_2\_alloc**](classDerivReconstructionCoeffs.md#variable-m_weights_patch_2_alloc)  <br>_The memory allocated for the weights of the function values on patch 2._  |




## Protected Static Attributes inherited from DerivReconstructionCoeffs

See [DerivReconstructionCoeffs](classDerivReconstructionCoeffs.md)

| Type | Name |
| ---: | :--- |
|  Extremity constexpr | [**s\_extremity\_1**](classDerivReconstructionCoeffs.md#variable-s_extremity_1)   = `extremity\_1`<br>_The extremity of patch 1 where the interface is located._  |
|  Extremity constexpr | [**s\_extremity\_2**](classDerivReconstructionCoeffs.md#variable-s_extremity_2)   = `extremity\_2`<br>_The extremity of patch 2 where the interface is located._  |




























## Protected Functions inherited from DerivReconstructionCoeffs

See [DerivReconstructionCoeffs](classDerivReconstructionCoeffs.md)

| Type | Name |
| ---: | :--- |
|  void | [**check\_additional\_interpolation\_points\_location**](classDerivReconstructionCoeffs.md#function-check_additional_interpolation_points_location) (IdxRange&lt; EdgePerpGrid &gt; const & idx\_range\_perp, Extremity const extremity) <br>_Check that the additional interpolation point is located in the boundary cell of the patch opposite the interface. An exception is raised if this is not the case._  |
|  void | [**check\_break\_points\_are\_interpolation\_points**](classDerivReconstructionCoeffs.md#function-check_break_points_are_interpolation_points) (bool const is\_cell\_bound\_with\_extra\_interpol\_pt, IdxRange&lt; EdgePerpGrid &gt; const & idx\_range\_perp, Extremity const extremity) const<br>_Check that the break points of the B-splines are interpolation points. An exception is raised if this is not the case._  |
|  void | [**correction\_boundary**](classDerivReconstructionCoeffs.md#function-correction_boundary) (int const n\_points, Idx1D const interface\_idx, std::array&lt; double, 3 &gt; & coeff\_deriv\_patch, std::array&lt; double, 3 &gt; & coeff\_deriv\_other\_patch, std::array&lt; host\_t&lt; DField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) &gt; &gt;, 3 &gt; & weights\_patch1, std::array&lt; host\_t&lt; DField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) &gt; &gt;, 3 &gt; & weights\_patch2) <br>_Correction of the local coefficients on the boundary for additional interpolation points in the middle of a cell._  |
|  double | [**get\_approx\_with\_change\_of\_sign**](classDerivReconstructionCoeffs.md#function-get_approx_with_change_of_sign) (DConstField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1), Kokkos::HostSpace, Layout1 &gt; const & function\_1, DConstField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2), Kokkos::HostSpace, Layout2 &gt; const & function\_2, int const sign\_function\_2) const<br>_Get the linear combination of the function values (c), changing the sign of the values on patch 2 if required._  |
|  std::tuple&lt; double, double &gt; | [**get\_cell\_lengths**](classDerivReconstructionCoeffs.md#function-get_cell_lengths) (Idx1D const & idx) const<br>_Get the lengths of the cell on the right side and the left side of the given index._  |
|  KOKKOS\_INLINE\_FUNCTION Idx&lt; Grid1D &gt; | [**get\_extremity\_idx**](classDerivReconstructionCoeffs.md#function-get_extremity_idx) (Extremity const extremity, IdxRange&lt; Grid1D &gt; const & idx\_range) const<br>_Get the index at the given extremity of an index range._  |
|  std::array&lt; double, 3 &gt; | [**get\_gammas**](classDerivReconstructionCoeffs.md#function-get_gammas) (double const cell\_length\_left, double const cell\_length\_right) const<br>_Compute the coefficients gammas\_i._  |
|  void | [**recursion**](classDerivReconstructionCoeffs.md#function-recursion) (int const n\_cells, Idx1D const interface\_idx, std::array&lt; double, 3 &gt; & coeff\_deriv\_patch, std::array&lt; double, 3 &gt; & coeff\_deriv\_other\_patch, std::array&lt; host\_t&lt; DField&lt; [**IdxRange1DPerp\_1**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_1) &gt; &gt;, 3 &gt; & weights\_patch1, std::array&lt; host\_t&lt; DField&lt; [**IdxRange1DPerp\_2**](classDerivReconstructionCoeffs.md#typedef-idxrange1dperp_2) &gt; &gt;, 3 &gt; & weights\_patch2) <br>_Recursion to compute the coefficients a, b and the weights omega._  |
|  void | [**set\_coefficients\_non\_uniform\_case**](classDerivReconstructionCoeffs.md#function-set_coefficients_non_uniform_case) () <br>_Compute the coefficients a, b and the weights omega applying the recursive formula._  |
|  void | [**set\_coefficients\_uniform\_per\_patch\_case**](classDerivReconstructionCoeffs.md#function-set_coefficients_uniform_per_patch_case) () <br>_Compute the coefficients a, b and the weights omega applying the explicit formula for the uniform case._  |




## Protected Static Functions inherited from DerivReconstructionCoeffs

See [DerivReconstructionCoeffs](classDerivReconstructionCoeffs.md)

| Type | Name |
| ---: | :--- |
|  double | [**get\_alpha**](classDerivReconstructionCoeffs.md#function-get_alpha) (double const cell\_length\_left, double const cell\_length\_right) <br>_Compute the coefficient alpha\_i._  |
|  double | [**get\_beta**](classDerivReconstructionCoeffs.md#function-get_beta) (double const cell\_length\_left, double const cell\_length\_right) <br>_Compute the coefficient beta\_i._  |
|  void | [**set\_local\_uniform\_weights**](classDerivReconstructionCoeffs.md#function-set_local_uniform_weights) (host\_t&lt; DField&lt; IdxRange&lt; EdgePerpGrid &gt; &gt; &gt; weights, double factor, double un\_other\_patch, Idx&lt; EdgePerpGrid &gt; idx\_interface) <br>_Apply the explicit formula for the uniform case to weights on one side of the interface._  |


## Detailed Description


This class specialises [**DerivReconstructionCoeffs**](classDerivReconstructionCoeffs.md) using the information stored in the [**Interface**](structInterface.md) (perpendicular grids, B-splines and extremities of the two Edges). It also provides patch-based accessors and handles the change of sign of the cross-derivatives when the orientations of the two patches do not agree.




**Template parameters:**


* `InterfaceType` The interface between two patches where we want to compute the derivatives. 




    
## Public Types Documentation




### typedef associated\_interface 

[_**Interface**_](structInterface.md) _between the two involved patches._
```C++
using InterfaceDerivCoeffs< InterfaceType >::associated_interface =  InterfaceType;
```




<hr>
## Public Functions Documentation




### function get\_approx\_cross\_deriv [1/2]

_Get the linear combination of the derivatives values (c)._ 
```C++
template<class Layout1, class Layout2>
inline double InterfaceDerivCoeffs::get_approx_cross_deriv (
    DConstField< IdxRange1DPerp_1, Kokkos::HostSpace, Layout1 > const & derivs_1,
    DConstField< IdxRange1DPerp_2, Kokkos::HostSpace, Layout2 > const & derivs_2
) const
```



This linear combination corresponds to an approximation of the interface cross-derivative.


It is more appropriate than get\_approx\_deriv to compute the cross-derivatives. The computation is the same. However, in addition, the sign of the derivatives on the patch 2 is changed in the case where the orientations of the two patches do not match. E.g. dx(f(-x)) = - dxf(-x)




**Parameters:**


* `derivs_1` Derivatives at the interpolation points on patch 1. 
* `derivs_2` Derivatives at the interpolation points on patch 2. 



**Returns:**

the linear combination of the function values (c). 





        

<hr>



### function get\_approx\_cross\_deriv [2/2]

_Get the linear combination of the derivatives values (c). See_ get\_approx\_cross\_deriv _._
```C++
template<class Layout1, class Layout2>
inline double InterfaceDerivCoeffs::get_approx_cross_deriv (
    DConstField< IdxRange1DPerp_2, Kokkos::HostSpace, Layout2 > const & derivs_2,
    DConstField< IdxRange1DPerp_1, Kokkos::HostSpace, Layout1 > const & derivs_1
) const
```





**Parameters:**


* `derivs_1` Derivatives at the interpolation points on patch 1. 
* `derivs_2` Derivatives at the interpolation points on patch 2. 



**Returns:**

the linear combination of the function values (c). 





        

<hr>



### function get\_coeff\_deriv\_on\_patch 

_Get the coefficient (a) or (b) in front of the derivative on the given patch of the given_ [_**Interface**_](structInterface.md) _._
```C++
template<class Patch>
inline double InterfaceDerivCoeffs::get_coeff_deriv_on_patch () const
```





**Returns:**

The value of the coefficient (a) or (b). 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/interface_derivative_coefficients.hpp`

