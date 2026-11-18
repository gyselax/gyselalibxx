

# Class InterfacesDerivativeCalculator

**template &lt;class Connectivity, class Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**InterfacesDerivativeCalculator**](classInterfacesDerivativeCalculator.md)



_Class to compute the interface derivatives along a given direction on all the interfaces with an approximation formula._ [More...](#detailed-description)

* `#include <interfaces_derivative_calculator.hpp>`



Inherits the following classes: [DerivativeCalculatorCommon](classDerivativeCalculatorCommon.md)














## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::Deriv&lt; std::conditional\_t&lt; std::is\_same\_v&lt; typename patch\_direction::Grid1, [**grid\_direction**](classInterfacesDerivativeCalculator.md#typedef-grid_direction) &gt;, typename patch\_direction::Dim2, typename patch\_direction::Dim1 &gt; &gt; | [**deriv\_par\_dim**](#typedef-deriv_par_dim)  <br>_Derivative dimension of the direction patch parallel to the mesh line._  |
| typedef ddc::Deriv&lt; typename Grid1D::continuous\_dimension\_type &gt; | [**deriv\_perp\_dim**](#typedef-deriv_perp_dim)  <br>_Derivative dimension of the direction patch perpendicular to the mesh line._  |
| typedef Grid1D | [**grid\_direction**](#typedef-grid_direction)  <br>_Grid 1D given to the class indicating the direction of the mesh line along which the derivatives are computed._  |
| typedef Idx&lt; [**deriv\_perp\_dim**](classInterfacesDerivativeCalculator.md#typedef-deriv_perp_dim), [**deriv\_par\_dim**](classInterfacesDerivativeCalculator.md#typedef-deriv_par_dim) &gt; | [**idx\_cross\_deriv**](#typedef-idx_cross_deriv)  <br>_Index on the two derivative dimensions of_ [_**patch\_direction**_](classInterfacesDerivativeCalculator.md#typedef-patch_direction) _._ |
| typedef Idx&lt; [**deriv\_perp\_dim**](classInterfacesDerivativeCalculator.md#typedef-deriv_perp_dim) &gt; | [**idx\_deriv\_perp**](#typedef-idx_deriv_perp)  <br>_Index on the derivative dimension perpendicular to the mesh line._  |
| typedef find\_patch\_t&lt; Grid1D, all\_patches &gt; | [**patch\_direction**](#typedef-patch_direction)  <br>[_**Patch**_](structPatch.md) _of the given direction grid._ |
| typedef PatchFirst | [**patch\_first**](#typedef-patch_first)  <br>_First patch of the collection of patches in the direction of the given Grid1D._  |


## Public Types inherited from DerivativeCalculatorCommon

See [DerivativeCalculatorCommon](classDerivativeCalculatorCommon.md)

| Type | Name |
| ---: | :--- |
| typedef MultipatchFieldFromSeq&lt; DerivConstFieldOnPatch\_host, relevant\_patches &gt; | [**MultipatchDerivConstField\_host**](classDerivativeCalculatorCommon.md#typedef-multipatchderivconstfield_host)  <br>_Type containing a DerivConstField per patch defined on CPU._  |
| typedef MultipatchFieldFromSeq&lt; DerivFieldOnPatch\_host, relevant\_patches &gt; | [**MultipatchDerivField\_host**](classDerivativeCalculatorCommon.md#typedef-multipatchderivfield_host)  <br>_Type containing a_ [_**DerivField**_](classDerivField.md) _per patch defined on CPU._ |
| typedef MultipatchTypeFromSeq&lt; IdxRangeOnPatch, relevant\_patches &gt; | [**MultipatchIdxRange**](classDerivativeCalculatorCommon.md#typedef-multipatchidxrange)  <br>_Type of the set of index ranges used in the operator._  |
| typedef ddc::Deriv&lt; std::conditional\_t&lt; std::is\_same\_v&lt; typename patch\_direction::Grid1, [**grid\_direction**](classDerivativeCalculatorCommon.md#typedef-grid_direction) &gt;, typename patch\_direction::Dim2, typename patch\_direction::Dim1 &gt; &gt; | [**deriv\_par\_dim**](classDerivativeCalculatorCommon.md#typedef-deriv_par_dim)  <br>_Derivative dimension of the direction patch parallel to the mesh line._  |
| typedef ddc::Deriv&lt; typename Grid1D::continuous\_dimension\_type &gt; | [**deriv\_perp\_dim**](classDerivativeCalculatorCommon.md#typedef-deriv_perp_dim)  <br>_Derivative dimension of the direction patch perpendicular to the mesh line._  |
| typedef Grid1D | [**grid\_direction**](classDerivativeCalculatorCommon.md#typedef-grid_direction)  <br>_Grid 1D given to the class indicating the direction of the mesh line along which the derivatives are computed._  |
| typedef Idx&lt; [**deriv\_perp\_dim**](classDerivativeCalculatorCommon.md#typedef-deriv_perp_dim), [**deriv\_par\_dim**](classDerivativeCalculatorCommon.md#typedef-deriv_par_dim) &gt; | [**idx\_cross\_deriv**](classDerivativeCalculatorCommon.md#typedef-idx_cross_deriv)  <br>_Index on the two derivative dimensions of_ [_**patch\_direction**_](classDerivativeCalculatorCommon.md#typedef-patch_direction) _._ |
| typedef Idx&lt; [**deriv\_perp\_dim**](classDerivativeCalculatorCommon.md#typedef-deriv_perp_dim) &gt; | [**idx\_deriv\_perp**](classDerivativeCalculatorCommon.md#typedef-idx_deriv_perp)  <br>_Index on the derivative dimension perpendicular to the mesh line._  |
| typedef find\_patch\_t&lt; Grid1D, all\_patches &gt; | [**patch\_direction**](classDerivativeCalculatorCommon.md#typedef-patch_direction)  <br>[_**Patch**_](structPatch.md) _of the given direction grid._ |
| typedef PatchFirst | [**patch\_first**](classDerivativeCalculatorCommon.md#typedef-patch_first)  <br>_First patch of the collection of patches in the direction of the given Grid1D._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**InterfacesDerivativeCalculator**](#function-interfacesderivativecalculator) (MultipatchIdxRange const & idx\_ranges, DerivCoeffs const &... deriv\_coeffs) <br>_Instantiate_ [_**InterfacesDerivativeCalculator**_](classInterfacesDerivativeCalculator.md) _._ |
|  void | [**operator()**](#function-operator) (IdxDeriv const & deriv\_order, MultipatchDerivField\_host functions\_and\_derivs, IdxPar const & idx\_par) <br>_Compute the vector C\_trunc to determine all the interface derivatives or cross-derivatives in the Grid1D direction at a given index._  |
|  void | [**operator()**](#function-operator_1) (IdxDeriv const & deriv\_order, MultipatchDerivField\_host functions\_and\_derivs) <br>_Compute the vector C\_trunc to determine all the interface derivatives in the Grid1D direction at all the indices of the first patch, or all the corner cross-derivatives in the Grid1D direction at all the last and the first indices of the first patch._  |


## Public Functions inherited from DerivativeCalculatorCommon

See [DerivativeCalculatorCommon](classDerivativeCalculatorCommon.md)

| Type | Name |
| ---: | :--- |
|   | [**DerivativeCalculatorCommon**](classDerivativeCalculatorCommon.md#function-derivativecalculatorcommon) ([**MultipatchIdxRange**](classDerivativeCalculatorCommon.md#typedef-multipatchidxrange) const & idx\_ranges, DerivCoeffs const &... deriv\_coeffs) <br>_Instantiate_ [_**DerivativeCalculatorCommon**_](classDerivativeCalculatorCommon.md) _._ |








## Protected Types inherited from DerivativeCalculatorCommon

See [DerivativeCalculatorCommon](classDerivativeCalculatorCommon.md)

| Type | Name |
| ---: | :--- |
| typedef collect\_grids\_on\_dim\_t&lt; find\_patch\_t&lt; Grid1D, all\_patches &gt;, Grid1D, [**interface\_sorted\_collection**](classDerivativeCalculatorCommon.md#typedef-interface_sorted_collection) &gt; | [**Grid1DSeq**](classDerivativeCalculatorCommon.md#typedef-grid1dseq)  <br>_Sequence of 1D grids in the direction of the given Grid1D over the patches._  |
| typedef typename ddc::Deriv&lt; typename [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)&lt; I &gt;::continuous\_dimension\_type &gt; | [**deriv\_par1**](classDerivativeCalculatorCommon.md#typedef-deriv_par1)  <br>_Tag the derivative dimension along parallel grid of the first patch._  |
| typedef typename ddc::Deriv&lt; typename [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)&lt; I &gt;::continuous\_dimension\_type &gt; | [**deriv\_par2**](classDerivativeCalculatorCommon.md#typedef-deriv_par2)  <br>_Tag the derivative dimension along parallel grid of the second patch._  |
| typedef typename ddc::Deriv&lt; typename [**grid\_perp1**](classDerivativeCalculatorCommon.md#typedef-grid_perp1)&lt; I &gt;::continuous\_dimension\_type &gt; | [**deriv\_perp1**](classDerivativeCalculatorCommon.md#typedef-deriv_perp1)  <br>_Tag the derivative dimension along perpendicular grid of the first patch._  |
| typedef typename ddc::Deriv&lt; typename [**grid\_perp2**](classDerivativeCalculatorCommon.md#typedef-grid_perp2)&lt; I &gt;::continuous\_dimension\_type &gt; | [**deriv\_perp2**](classDerivativeCalculatorCommon.md#typedef-deriv_perp2)  <br>_Tag the derivative dimension along perpendicular grid of the second patch._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge1::parallel\_grid | [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)  <br>_Tag the parallel grid of the first patch to the interface._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge2::parallel\_grid | [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)  <br>_Tag the parallel grid of the second patch to the interface._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge1::perpendicular\_grid | [**grid\_perp1**](classDerivativeCalculatorCommon.md#typedef-grid_perp1)  <br>_Tag the perpendicular grid of the first patch to the interface._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge2::perpendicular\_grid | [**grid\_perp2**](classDerivativeCalculatorCommon.md#typedef-grid_perp2)  <br>_Tag the perpendicular grid of the second patch to the interface._  |
| typedef Idx&lt; ddc::Deriv&lt; typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Dim1 &gt;, typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Grid1, ddc::Deriv&lt; typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Dim2 &gt;, typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Grid2 &gt; | [**idx\_cross\_deriv1**](classDerivativeCalculatorCommon.md#typedef-idx_cross_deriv1)  <br>_Tag the index type for cross-derivative on the first patch._  |
| typedef Idx&lt; ddc::Deriv&lt; typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Dim1 &gt;, typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Grid1, ddc::Deriv&lt; typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Dim2 &gt;, typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Grid2 &gt; | [**idx\_cross\_deriv2**](classDerivativeCalculatorCommon.md#typedef-idx_cross_deriv2)  <br>_Tag the index type for cross-derivatives on the second patch._  |
| typedef Idx&lt; [**deriv\_perp1**](classDerivativeCalculatorCommon.md#typedef-deriv_perp1)&lt; I &gt;, typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Grid1, typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Grid2 &gt; | [**idx\_deriv\_perp1**](classDerivativeCalculatorCommon.md#typedef-idx_deriv_perp1)  <br>_Tag the index type for derivatives along the perpendicular grid on the first patch._  |
| typedef Idx&lt; [**deriv\_perp2**](classDerivativeCalculatorCommon.md#typedef-deriv_perp2)&lt; I &gt;, typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Grid1, typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Grid2 &gt; | [**idx\_deriv\_perp2**](classDerivativeCalculatorCommon.md#typedef-idx_deriv_perp2)  <br>_Tag the index type for derivatives along the perpendicular grid on the second patch._  |
| typedef ddc::type\_seq\_remove\_t&lt; [**interface\_sorted\_collection**](classDerivativeCalculatorCommon.md#typedef-interface_sorted_collection), outer\_interface\_collection &gt; | [**inner\_interface\_collection**](classDerivativeCalculatorCommon.md#typedef-inner_interface_collection)  <br>_Tag for all the interfaces in the Grid1D direction without the interfaces with_ [_**OutsideEdge**_](structOutsideEdge.md) _._ |
| typedef ddc::type\_seq\_element\_t&lt; I, [**inner\_interface\_collection**](classDerivativeCalculatorCommon.md#typedef-inner_interface_collection) &gt; | [**interface**](classDerivativeCalculatorCommon.md#typedef-interface)  <br>_Tag the interface with the correct patch order._  |
| typedef typename Connectivity::template get\_all\_interfaces\_along\_direction\_t&lt; Grid1D &gt; | [**interface\_sorted\_collection**](classDerivativeCalculatorCommon.md#typedef-interface_sorted_collection)  <br>_Tag for all the sorted interfaces with the correct patch order._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge1::associated\_patch | [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)  <br>_Tag the first patch of the stored interface._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge2::associated\_patch | [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)  <br>_Tag the second patch of the stored interface._  |
| typedef find\_associated\_interface\_t&lt; typename [**interface**](classDerivativeCalculatorCommon.md#typedef-interface)&lt; I &gt;::Edge1, all\_interface\_collection &gt; | [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)  <br>_Tag the equivalent interface defined in the main file, i.e. given to_ [_**MultipatchConnectivity**_](classMultipatchConnectivity.md) _._ |








## Protected Attributes inherited from DerivativeCalculatorCommon

See [DerivativeCalculatorCommon](classDerivativeCalculatorCommon.md)

| Type | Name |
| ---: | :--- |
|  [**MultipatchIdxRange**](classDerivativeCalculatorCommon.md#typedef-multipatchidxrange) | [**m\_idx\_ranges**](classDerivativeCalculatorCommon.md#variable-m_idx_ranges)  <br>_MultipatchType on IdxRangeOnPatch to collect the index ranges on the patches._  |




## Protected Static Attributes inherited from DerivativeCalculatorCommon

See [DerivativeCalculatorCommon](classDerivativeCalculatorCommon.md)

| Type | Name |
| ---: | :--- |
|  Extremity | [**extremity\_1**](classDerivativeCalculatorCommon.md#variable-extremity_1)   = `[**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt;I&gt;::Edge1::extremity`<br>_Tag the connection extremity of the first patch of the interface._  |
|  Extremity | [**extremity\_2**](classDerivativeCalculatorCommon.md#variable-extremity_2)   = `[**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt;I&gt;::Edge2::extremity`<br>_Tag the connection extremity of the second patch of the interface._  |
|  bool | [**is\_periodic**](classDerivativeCalculatorCommon.md#variable-is_periodic)   = `/* multi line expression */`<br>_Boolean to indicate the periodicty of the global domain._  |
|  std::size\_t | [**n\_inner\_interfaces**](classDerivativeCalculatorCommon.md#variable-n_inner_interfaces)   = `ddc::type\_seq\_size\_v&lt;[**inner\_interface\_collection**](classDerivativeCalculatorCommon.md#typedef-inner_interface_collection)&gt;`<br>_Number of inner interfaces (fixing the number of elements in the matrix and vectors)._  |
|  bool | [**patch\_order\_matches\_interfaces**](classDerivativeCalculatorCommon.md#variable-patch_order_matches_interfaces)   = `/* multi line expression */`<br>_Boolean if the order of Patch1\|Patch2 of the interface matches with the sorted 1D grid sequence._  |




























## Protected Functions inherited from DerivativeCalculatorCommon

See [DerivativeCalculatorCommon](classDerivativeCalculatorCommon.md)

| Type | Name |
| ---: | :--- |
|  double | [**get\_approximated\_deriv**](classDerivativeCalculatorCommon.md#function-get_approximated_deriv-12) (Idx&lt; [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)&lt; I &gt; &gt; const & idx\_to\_slice\_1, Idx&lt; [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)&lt; I &gt; &gt; const & idx\_to\_slice\_2, [**MultipatchDerivConstField\_host**](classDerivativeCalculatorCommon.md#typedef-multipatchderivconstfield_host) functions\_and\_derivs) const<br>_Compute the approximated derivative (or coefficient "c") at the interface along the mesh line indexed at the given indices._  |
|  double | [**get\_approximated\_deriv**](classDerivativeCalculatorCommon.md#function-get_approximated_deriv-12) (Idx&lt; [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)&lt; I &gt; &gt; const & idx\_to\_slice\_1, Idx&lt; [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)&lt; I &gt; &gt; const & idx\_to\_slice\_2, [**MultipatchDerivConstField\_host**](classDerivativeCalculatorCommon.md#typedef-multipatchderivconstfield_host) functions\_and\_derivs) const<br>_Compute the approximated cross-derivative (or coefficient "c") at the interface along the mesh line indexed at the given indices. Work only for indices in the corners of the patches (where the cross-derivatives are supposed to be computed)._  |
|  [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md)&lt; [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt; &gt; const & | [**get\_deriv\_coeffs**](classDerivativeCalculatorCommon.md#function-get_deriv_coeffs) () const<br>_Get the_ [_**InterfaceDerivCoeffs**_](classInterfaceDerivCoeffs.md) _defined on the I-th inner interface._ |
|  std::conditional\_t&lt; is\_patch1, Idx&lt; [**grid\_perp1**](classDerivativeCalculatorCommon.md#typedef-grid_perp1)&lt; I &gt; &gt;, Idx&lt; [**grid\_perp2**](classDerivativeCalculatorCommon.md#typedef-grid_perp2)&lt; I &gt; &gt; &gt; | [**get\_idx\_interface**](classDerivativeCalculatorCommon.md#function-get_idx_interface) () const<br>_Get the index of the given index range slice corresponding to one of the extremities of the 1D grid. If is\_on\_interface is true, the extremity corresponds to the interface. If is\_on\_interface is false, the extremity corresponds to the other bound of the interface (next or previous)._  |
|  void | [**set\_slice\_indexes**](classDerivativeCalculatorCommon.md#function-set_slice_indexes) (int & value\_idx\_to\_slice, Idx&lt; [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)&lt; I &gt; &gt; & idx\_to\_slice\_1, Idx&lt; [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)&lt; I &gt; &gt; & idx\_to\_slice\_2) const<br>_Get the indices on Patch1 and Patch2 of the given interface. The index on Patch1 corresponds to the value given in input (if the interface is well oriented, otherwise it is Patch2). The index on Patch2 is the equivalent index of the one on Patch1. It also updates the value given in input with the value of the index on Patch2 for the next interface (if the next interface is well oriented, otherwise it is Patch1)._  |




## Protected Static Functions inherited from DerivativeCalculatorCommon

See [DerivativeCalculatorCommon](classDerivativeCalculatorCommon.md)

| Type | Name |
| ---: | :--- |
|  int | [**to\_sign**](classDerivativeCalculatorCommon.md#function-to_sign) (bool agree) <br>_Transform a boolean into a integer. If true, the sign is positive. If false, the sign is negative._  |


## Detailed Description


This operator is implemented for conforming equivalent global meshes.


When we call operator(), the operator loops over the interface to compute the vector C\_trunc to determine all the interface derivatives and update the [**DerivField**](classDerivField.md) given in input.


The vector C\_trunc is given by the coefficients (see README), 
\[\{\sum_{k = - N_{reduc}}^{N_{reduc}} \omega_{k, N_{reduc}, N_{reduc}}^{i_I} f_{i_I+k} \}_I\]



with
* \(N_{reduc}\) the number of chosen cells for the approximation,
* \(i_I\) the index of the grid at the Ith interface,
* \(f_{i_I+k}\) the function values if we compute the first derivatives, or the first derivatives if we compute the cross-derivatives.






**Template parameters:**


* `Connectivity` A [**MultipatchConnectivity**](classMultipatchConnectivity.md) class describing all the patch connections. 
* `Grid1D` A given direction. 




    
## Public Types Documentation




### typedef deriv\_par\_dim 

_Derivative dimension of the direction patch parallel to the mesh line._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::deriv_par_dim =  ddc::Deriv<std::conditional_t<
            std::is_same_v<typename patch_direction::Grid1, grid_direction>,
            typename patch_direction::Dim2,
            typename patch_direction::Dim1>>;
```




<hr>



### typedef deriv\_perp\_dim 

_Derivative dimension of the direction patch perpendicular to the mesh line._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::deriv_perp_dim =  ddc::Deriv<typename Grid1D::continuous_dimension_type>;
```




<hr>



### typedef grid\_direction 

_Grid 1D given to the class indicating the direction of the mesh line along which the derivatives are computed._ 
```C++
using InterfacesDerivativeCalculator< Connectivity, Grid1D >::grid_direction =  Grid1D;
```




<hr>



### typedef idx\_cross\_deriv 

_Index on the two derivative dimensions of_ [_**patch\_direction**_](classInterfacesDerivativeCalculator.md#typedef-patch_direction) _._
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::idx_cross_deriv =  Idx<deriv_perp_dim, deriv_par_dim>;
```




<hr>



### typedef idx\_deriv\_perp 

_Index on the derivative dimension perpendicular to the mesh line._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::idx_deriv_perp =  Idx<deriv_perp_dim>;
```




<hr>



### typedef patch\_direction 

[_**Patch**_](structPatch.md) _of the given direction grid._
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::patch_direction =  find_patch_t<Grid1D, all_patches>;
```




<hr>



### typedef patch\_first 

_First patch of the collection of patches in the direction of the given Grid1D._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::patch_first =  PatchFirst;
```




<hr>
## Public Functions Documentation




### function InterfacesDerivativeCalculator 

_Instantiate_ [_**InterfacesDerivativeCalculator**_](classInterfacesDerivativeCalculator.md) _._
```C++
template<class... DerivCoeffs>
inline explicit InterfacesDerivativeCalculator::InterfacesDerivativeCalculator (
    MultipatchIdxRange const & idx_ranges,
    DerivCoeffs const &... deriv_coeffs
) 
```



The [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md) are stored by reference. They must outlive this object.




**Parameters:**


* `idx_ranges` MultipatchType collection of index ranges defined on the given list of patches. 
* `deriv_coeffs` The [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md) for each inner interface in the given Grid1D direction (i.e. excluding interfaces with an [**OutsideEdge**](structOutsideEdge.md)). They can be given in any order. 




        

<hr>



### function operator() 

_Compute the vector C\_trunc to determine all the interface derivatives or cross-derivatives in the Grid1D direction at a given index._ 
```C++
template<class IdxDeriv, class IdxPar>
inline void InterfacesDerivativeCalculator::operator() (
    IdxDeriv const & deriv_order,
    MultipatchDerivField_host functions_and_derivs,
    IdxPar const & idx_par
) 
```



For the derivatives, it uses the function values to compute the first-derivatives perpendicular to the interfaces. For the cross-derivatives, it uses the first-derivatives to compute the cross-derivatives in a perpendicular direction to the interfaces.




**Template parameters:**


* `IdxDeriv` Type for the given derivative order, either on the perpendicular direction for the first derivatives or both for the cross-derivatives. 
* `IdxPar` Type for the given index. It can be on the first or the second grid of the first patch. 



**Parameters:**


* `deriv_order` Indicate if we compute the first derivatives or the cross-derivatives. Its value is not taken into account. 
* `functions_and_derivs` Collection fields with its first derivatives and cross-derivatives. 
* `idx_par` Index of the line in the Grid1D direction where we compute all the interface derivatives. 




        

<hr>



### function operator() 

_Compute the vector C\_trunc to determine all the interface derivatives in the Grid1D direction at all the indices of the first patch, or all the corner cross-derivatives in the Grid1D direction at all the last and the first indices of the first patch._ 
```C++
template<class IdxDeriv>
inline void InterfacesDerivativeCalculator::operator() (
    IdxDeriv const & deriv_order,
    MultipatchDerivField_host functions_and_derivs
) 
```



For the derivatives, it uses the function values to compute the first-derivatives perpendicular to the interfaces. For the cross-derivatives, it uses the first-derivatives to compute the cross-derivatives in a perpendicular direction to the interfaces.




**Template parameters:**


* `IdxDeriv` Type for the given derivative order, either on the perpendicular direction for the first derivatives or both for the cross-derivatives. 



**Parameters:**


* `deriv_order` Indicate if we compute the first derivatives or the cross-derivatives. Its value is not taken into account. 
* `functions_and_derivs` Collection fields with its first derivatives and cross-derivatives we want to compute. 



**Warning:**

The global mesh is supposed to be conforming. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/interfaces_derivative_calculator.hpp`

