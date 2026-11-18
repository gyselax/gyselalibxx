

# Class DerivativeCalculatorCommon

**template &lt;class Connectivity, class Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**DerivativeCalculatorCommon**](classDerivativeCalculatorCommon.md)



_Abstract class to hold the interface derivative coefficients and the index ranges. This is the superclass for interface derivative calculators defined in this folder._ [More...](#detailed-description)

* `#include <derivative_calculator_common.hpp>`





Inherited by the following classes: [InterfacesDerivativeCalculator](classInterfacesDerivativeCalculator.md),  [InterfacesExactDerivativeCalculator](classInterfacesExactDerivativeCalculator.md)












## Public Types

| Type | Name |
| ---: | :--- |
| typedef MultipatchFieldFromSeq&lt; DerivConstFieldOnPatch\_host, relevant\_patches &gt; | [**MultipatchDerivConstField\_host**](#typedef-multipatchderivconstfield_host)  <br>_Type containing a DerivConstField per patch defined on CPU._  |
| typedef MultipatchFieldFromSeq&lt; DerivFieldOnPatch\_host, relevant\_patches &gt; | [**MultipatchDerivField\_host**](#typedef-multipatchderivfield_host)  <br>_Type containing a_ [_**DerivField**_](classDerivField.md) _per patch defined on CPU._ |
| typedef MultipatchTypeFromSeq&lt; IdxRangeOnPatch, relevant\_patches &gt; | [**MultipatchIdxRange**](#typedef-multipatchidxrange)  <br>_Type of the set of index ranges used in the operator._  |
| typedef ddc::Deriv&lt; std::conditional\_t&lt; std::is\_same\_v&lt; typename patch\_direction::Grid1, [**grid\_direction**](classDerivativeCalculatorCommon.md#typedef-grid_direction) &gt;, typename patch\_direction::Dim2, typename patch\_direction::Dim1 &gt; &gt; | [**deriv\_par\_dim**](#typedef-deriv_par_dim)  <br>_Derivative dimension of the direction patch parallel to the mesh line._  |
| typedef ddc::Deriv&lt; typename Grid1D::continuous\_dimension\_type &gt; | [**deriv\_perp\_dim**](#typedef-deriv_perp_dim)  <br>_Derivative dimension of the direction patch perpendicular to the mesh line._  |
| typedef Grid1D | [**grid\_direction**](#typedef-grid_direction)  <br>_Grid 1D given to the class indicating the direction of the mesh line along which the derivatives are computed._  |
| typedef Idx&lt; [**deriv\_perp\_dim**](classDerivativeCalculatorCommon.md#typedef-deriv_perp_dim), [**deriv\_par\_dim**](classDerivativeCalculatorCommon.md#typedef-deriv_par_dim) &gt; | [**idx\_cross\_deriv**](#typedef-idx_cross_deriv)  <br>_Index on the two derivative dimensions of_ [_**patch\_direction**_](classDerivativeCalculatorCommon.md#typedef-patch_direction) _._ |
| typedef Idx&lt; [**deriv\_perp\_dim**](classDerivativeCalculatorCommon.md#typedef-deriv_perp_dim) &gt; | [**idx\_deriv\_perp**](#typedef-idx_deriv_perp)  <br>_Index on the derivative dimension perpendicular to the mesh line._  |
| typedef find\_patch\_t&lt; Grid1D, all\_patches &gt; | [**patch\_direction**](#typedef-patch_direction)  <br>[_**Patch**_](structPatch.md) _of the given direction grid._ |
| typedef PatchFirst | [**patch\_first**](#typedef-patch_first)  <br>_First patch of the collection of patches in the direction of the given Grid1D._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**DerivativeCalculatorCommon**](#function-derivativecalculatorcommon) ([**MultipatchIdxRange**](classDerivativeCalculatorCommon.md#typedef-multipatchidxrange) const & idx\_ranges, DerivCoeffs const &... deriv\_coeffs) <br>_Instantiate_ [_**DerivativeCalculatorCommon**_](classDerivativeCalculatorCommon.md) _._ |




## Protected Types

| Type | Name |
| ---: | :--- |
| typedef collect\_grids\_on\_dim\_t&lt; find\_patch\_t&lt; Grid1D, all\_patches &gt;, Grid1D, [**interface\_sorted\_collection**](classDerivativeCalculatorCommon.md#typedef-interface_sorted_collection) &gt; | [**Grid1DSeq**](#typedef-grid1dseq)  <br>_Sequence of 1D grids in the direction of the given Grid1D over the patches._  |
| typedef typename ddc::Deriv&lt; typename [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)&lt; I &gt;::continuous\_dimension\_type &gt; | [**deriv\_par1**](#typedef-deriv_par1)  <br>_Tag the derivative dimension along parallel grid of the first patch._  |
| typedef typename ddc::Deriv&lt; typename [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)&lt; I &gt;::continuous\_dimension\_type &gt; | [**deriv\_par2**](#typedef-deriv_par2)  <br>_Tag the derivative dimension along parallel grid of the second patch._  |
| typedef typename ddc::Deriv&lt; typename [**grid\_perp1**](classDerivativeCalculatorCommon.md#typedef-grid_perp1)&lt; I &gt;::continuous\_dimension\_type &gt; | [**deriv\_perp1**](#typedef-deriv_perp1)  <br>_Tag the derivative dimension along perpendicular grid of the first patch._  |
| typedef typename ddc::Deriv&lt; typename [**grid\_perp2**](classDerivativeCalculatorCommon.md#typedef-grid_perp2)&lt; I &gt;::continuous\_dimension\_type &gt; | [**deriv\_perp2**](#typedef-deriv_perp2)  <br>_Tag the derivative dimension along perpendicular grid of the second patch._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge1::parallel\_grid | [**grid\_par1**](#typedef-grid_par1)  <br>_Tag the parallel grid of the first patch to the interface._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge2::parallel\_grid | [**grid\_par2**](#typedef-grid_par2)  <br>_Tag the parallel grid of the second patch to the interface._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge1::perpendicular\_grid | [**grid\_perp1**](#typedef-grid_perp1)  <br>_Tag the perpendicular grid of the first patch to the interface._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge2::perpendicular\_grid | [**grid\_perp2**](#typedef-grid_perp2)  <br>_Tag the perpendicular grid of the second patch to the interface._  |
| typedef Idx&lt; ddc::Deriv&lt; typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Dim1 &gt;, typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Grid1, ddc::Deriv&lt; typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Dim2 &gt;, typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Grid2 &gt; | [**idx\_cross\_deriv1**](#typedef-idx_cross_deriv1)  <br>_Tag the index type for cross-derivative on the first patch._  |
| typedef Idx&lt; ddc::Deriv&lt; typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Dim1 &gt;, typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Grid1, ddc::Deriv&lt; typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Dim2 &gt;, typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Grid2 &gt; | [**idx\_cross\_deriv2**](#typedef-idx_cross_deriv2)  <br>_Tag the index type for cross-derivatives on the second patch._  |
| typedef Idx&lt; [**deriv\_perp1**](classDerivativeCalculatorCommon.md#typedef-deriv_perp1)&lt; I &gt;, typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Grid1, typename [**patch1**](classDerivativeCalculatorCommon.md#typedef-patch1)&lt; I &gt;::Grid2 &gt; | [**idx\_deriv\_perp1**](#typedef-idx_deriv_perp1)  <br>_Tag the index type for derivatives along the perpendicular grid on the first patch._  |
| typedef Idx&lt; [**deriv\_perp2**](classDerivativeCalculatorCommon.md#typedef-deriv_perp2)&lt; I &gt;, typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Grid1, typename [**patch2**](classDerivativeCalculatorCommon.md#typedef-patch2)&lt; I &gt;::Grid2 &gt; | [**idx\_deriv\_perp2**](#typedef-idx_deriv_perp2)  <br>_Tag the index type for derivatives along the perpendicular grid on the second patch._  |
| typedef ddc::type\_seq\_remove\_t&lt; [**interface\_sorted\_collection**](classDerivativeCalculatorCommon.md#typedef-interface_sorted_collection), outer\_interface\_collection &gt; | [**inner\_interface\_collection**](#typedef-inner_interface_collection)  <br>_Tag for all the interfaces in the Grid1D direction without the interfaces with_ [_**OutsideEdge**_](structOutsideEdge.md) _._ |
| typedef ddc::type\_seq\_element\_t&lt; I, [**inner\_interface\_collection**](classDerivativeCalculatorCommon.md#typedef-inner_interface_collection) &gt; | [**interface**](#typedef-interface)  <br>_Tag the interface with the correct patch order._  |
| typedef typename Connectivity::template get\_all\_interfaces\_along\_direction\_t&lt; Grid1D &gt; | [**interface\_sorted\_collection**](#typedef-interface_sorted_collection)  <br>_Tag for all the sorted interfaces with the correct patch order._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge1::associated\_patch | [**patch1**](#typedef-patch1)  <br>_Tag the first patch of the stored interface._  |
| typedef typename [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt;::Edge2::associated\_patch | [**patch2**](#typedef-patch2)  <br>_Tag the second patch of the stored interface._  |
| typedef find\_associated\_interface\_t&lt; typename [**interface**](classDerivativeCalculatorCommon.md#typedef-interface)&lt; I &gt;::Edge1, all\_interface\_collection &gt; | [**stored\_interface**](#typedef-stored_interface)  <br>_Tag the equivalent interface defined in the main file, i.e. given to_ [_**MultipatchConnectivity**_](classMultipatchConnectivity.md) _._ |




## Protected Attributes

| Type | Name |
| ---: | :--- |
|  [**MultipatchIdxRange**](classDerivativeCalculatorCommon.md#typedef-multipatchidxrange) | [**m\_idx\_ranges**](#variable-m_idx_ranges)  <br>_MultipatchType on IdxRangeOnPatch to collect the index ranges on the patches._  |


## Protected Static Attributes

| Type | Name |
| ---: | :--- |
|  Extremity | [**extremity\_1**](#variable-extremity_1)   = `[**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt;I&gt;::Edge1::extremity`<br>_Tag the connection extremity of the first patch of the interface._  |
|  Extremity | [**extremity\_2**](#variable-extremity_2)   = `[**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt;I&gt;::Edge2::extremity`<br>_Tag the connection extremity of the second patch of the interface._  |
|  bool | [**is\_periodic**](#variable-is_periodic)   = `/* multi line expression */`<br>_Boolean to indicate the periodicty of the global domain._  |
|  std::size\_t | [**n\_inner\_interfaces**](#variable-n_inner_interfaces)   = `ddc::type\_seq\_size\_v&lt;[**inner\_interface\_collection**](classDerivativeCalculatorCommon.md#typedef-inner_interface_collection)&gt;`<br>_Number of inner interfaces (fixing the number of elements in the matrix and vectors)._  |
|  bool | [**patch\_order\_matches\_interfaces**](#variable-patch_order_matches_interfaces)   = `/* multi line expression */`<br>_Boolean if the order of Patch1\|Patch2 of the interface matches with the sorted 1D grid sequence._  |














## Protected Functions

| Type | Name |
| ---: | :--- |
|  double | [**get\_approximated\_deriv**](#function-get_approximated_deriv-12) (Idx&lt; [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)&lt; I &gt; &gt; const & idx\_to\_slice\_1, Idx&lt; [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)&lt; I &gt; &gt; const & idx\_to\_slice\_2, [**MultipatchDerivConstField\_host**](classDerivativeCalculatorCommon.md#typedef-multipatchderivconstfield_host) functions\_and\_derivs) const<br>_Compute the approximated derivative (or coefficient "c") at the interface along the mesh line indexed at the given indices._  |
|  double | [**get\_approximated\_deriv**](#function-get_approximated_deriv-12) (Idx&lt; [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)&lt; I &gt; &gt; const & idx\_to\_slice\_1, Idx&lt; [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)&lt; I &gt; &gt; const & idx\_to\_slice\_2, [**MultipatchDerivConstField\_host**](classDerivativeCalculatorCommon.md#typedef-multipatchderivconstfield_host) functions\_and\_derivs) const<br>_Compute the approximated cross-derivative (or coefficient "c") at the interface along the mesh line indexed at the given indices. Work only for indices in the corners of the patches (where the cross-derivatives are supposed to be computed)._  |
|  [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md)&lt; [**stored\_interface**](classDerivativeCalculatorCommon.md#typedef-stored_interface)&lt; I &gt; &gt; const & | [**get\_deriv\_coeffs**](#function-get_deriv_coeffs) () const<br>_Get the_ [_**InterfaceDerivCoeffs**_](classInterfaceDerivCoeffs.md) _defined on the I-th inner interface._ |
|  std::conditional\_t&lt; is\_patch1, Idx&lt; [**grid\_perp1**](classDerivativeCalculatorCommon.md#typedef-grid_perp1)&lt; I &gt; &gt;, Idx&lt; [**grid\_perp2**](classDerivativeCalculatorCommon.md#typedef-grid_perp2)&lt; I &gt; &gt; &gt; | [**get\_idx\_interface**](#function-get_idx_interface) () const<br>_Get the index of the given index range slice corresponding to one of the extremities of the 1D grid. If is\_on\_interface is true, the extremity corresponds to the interface. If is\_on\_interface is false, the extremity corresponds to the other bound of the interface (next or previous)._  |
|  void | [**set\_slice\_indexes**](#function-set_slice_indexes) (int & value\_idx\_to\_slice, Idx&lt; [**grid\_par1**](classDerivativeCalculatorCommon.md#typedef-grid_par1)&lt; I &gt; &gt; & idx\_to\_slice\_1, Idx&lt; [**grid\_par2**](classDerivativeCalculatorCommon.md#typedef-grid_par2)&lt; I &gt; &gt; & idx\_to\_slice\_2) const<br>_Get the indices on Patch1 and Patch2 of the given interface. The index on Patch1 corresponds to the value given in input (if the interface is well oriented, otherwise it is Patch2). The index on Patch2 is the equivalent index of the one on Patch1. It also updates the value given in input with the value of the index on Patch2 for the next interface (if the next interface is well oriented, otherwise it is Patch1)._  |


## Protected Static Functions

| Type | Name |
| ---: | :--- |
|  int | [**to\_sign**](#function-to_sign) (bool agree) <br>_Transform a boolean into a integer. If true, the sign is positive. If false, the sign is negative._  |


## Detailed Description




**Template parameters:**


* `Connectivity` A [**MultipatchConnectivity**](classMultipatchConnectivity.md) class describing all the patch connections. 
* `Grid1D` A given direction. 




    
## Public Types Documentation




### typedef MultipatchDerivConstField\_host 

_Type containing a DerivConstField per patch defined on CPU._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::MultipatchDerivConstField_host =  MultipatchFieldFromSeq<DerivConstFieldOnPatch_host, relevant_patches>;
```




<hr>



### typedef MultipatchDerivField\_host 

_Type containing a_ [_**DerivField**_](classDerivField.md) _per patch defined on CPU._
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::MultipatchDerivField_host =  MultipatchFieldFromSeq<DerivFieldOnPatch_host, relevant_patches>;
```




<hr>



### typedef MultipatchIdxRange 

_Type of the set of index ranges used in the operator._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::MultipatchIdxRange =  MultipatchTypeFromSeq<IdxRangeOnPatch, relevant_patches>;
```




<hr>



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
using DerivativeCalculatorCommon< Connectivity, Grid1D >::grid_direction =  Grid1D;
```




<hr>



### typedef idx\_cross\_deriv 

_Index on the two derivative dimensions of_ [_**patch\_direction**_](classDerivativeCalculatorCommon.md#typedef-patch_direction) _._
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




### function DerivativeCalculatorCommon 

_Instantiate_ [_**DerivativeCalculatorCommon**_](classDerivativeCalculatorCommon.md) _._
```C++
template<class... DerivCoeffs>
inline explicit DerivativeCalculatorCommon::DerivativeCalculatorCommon (
    MultipatchIdxRange const & idx_ranges,
    DerivCoeffs const &... deriv_coeffs
) 
```



The [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md) are stored by reference. They must outlive this object.




**Parameters:**


* `idx_ranges` MultipatchType collection of index ranges defined on the given list of patches. 
* `deriv_coeffs` The [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md) for each inner interface in the given Grid1D direction (i.e. excluding interfaces with an [**OutsideEdge**](structOutsideEdge.md)). They can be given in any order. 




        

<hr>
## Protected Types Documentation




### typedef Grid1DSeq 

_Sequence of 1D grids in the direction of the given Grid1D over the patches._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::Grid1DSeq =  collect_grids_on_dim_t<
            find_patch_t<Grid1D, all_patches>,
            Grid1D,
            interface_sorted_collection>;
```




<hr>



### typedef deriv\_par1 

_Tag the derivative dimension along parallel grid of the first patch._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::deriv_par1 =  typename ddc::Deriv<typename grid_par1<I>::continuous_dimension_type>;
```




<hr>



### typedef deriv\_par2 

_Tag the derivative dimension along parallel grid of the second patch._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::deriv_par2 =  typename ddc::Deriv<typename grid_par2<I>::continuous_dimension_type>;
```




<hr>



### typedef deriv\_perp1 

_Tag the derivative dimension along perpendicular grid of the first patch._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::deriv_perp1 =  typename ddc::Deriv<typename grid_perp1<I>::continuous_dimension_type>;
```




<hr>



### typedef deriv\_perp2 

_Tag the derivative dimension along perpendicular grid of the second patch._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::deriv_perp2 =  typename ddc::Deriv<typename grid_perp2<I>::continuous_dimension_type>;
```




<hr>



### typedef grid\_par1 

_Tag the parallel grid of the first patch to the interface._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::grid_par1 =  typename stored_interface<I>::Edge1::parallel_grid;
```




<hr>



### typedef grid\_par2 

_Tag the parallel grid of the second patch to the interface._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::grid_par2 =  typename stored_interface<I>::Edge2::parallel_grid;
```




<hr>



### typedef grid\_perp1 

_Tag the perpendicular grid of the first patch to the interface._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::grid_perp1 =  typename stored_interface<I>::Edge1::perpendicular_grid;
```




<hr>



### typedef grid\_perp2 

_Tag the perpendicular grid of the second patch to the interface._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::grid_perp2 =  typename stored_interface<I>::Edge2::perpendicular_grid;
```




<hr>



### typedef idx\_cross\_deriv1 

_Tag the index type for cross-derivative on the first patch._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::idx_cross_deriv1 =  Idx<ddc::Deriv<typename patch1<I>::Dim1>,
                  typename patch1<I>::Grid1,
                  ddc::Deriv<typename patch1<I>::Dim2>,
                  typename patch1<I>::Grid2>;
```




<hr>



### typedef idx\_cross\_deriv2 

_Tag the index type for cross-derivatives on the second patch._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::idx_cross_deriv2 =  Idx<ddc::Deriv<typename patch2<I>::Dim1>,
                  typename patch2<I>::Grid1,
                  ddc::Deriv<typename patch2<I>::Dim2>,
                  typename patch2<I>::Grid2>;
```




<hr>



### typedef idx\_deriv\_perp1 

_Tag the index type for derivatives along the perpendicular grid on the first patch._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::idx_deriv_perp1 =  Idx<deriv_perp1<I>, typename patch1<I>::Grid1, typename patch1<I>::Grid2>;
```




<hr>



### typedef idx\_deriv\_perp2 

_Tag the index type for derivatives along the perpendicular grid on the second patch._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::idx_deriv_perp2 =  Idx<deriv_perp2<I>, typename patch2<I>::Grid1, typename patch2<I>::Grid2>;
```




<hr>



### typedef inner\_interface\_collection 

_Tag for all the interfaces in the Grid1D direction without the interfaces with_ [_**OutsideEdge**_](structOutsideEdge.md) _._
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::inner_interface_collection =  ddc::type_seq_remove_t<interface_sorted_collection, outer_interface_collection>;
```




<hr>



### typedef interface 

_Tag the interface with the correct patch order._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::interface =  ddc::type_seq_element_t<I, inner_interface_collection>;
```




<hr>



### typedef interface\_sorted\_collection 

_Tag for all the sorted interfaces with the correct patch order._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::interface_sorted_collection = 
            typename Connectivity::template get_all_interfaces_along_direction_t<Grid1D>;
```




<hr>



### typedef patch1 

_Tag the first patch of the stored interface._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::patch1 =  typename stored_interface<I>::Edge1::associated_patch;
```




<hr>



### typedef patch2 

_Tag the second patch of the stored interface._ 
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::patch2 =  typename stored_interface<I>::Edge2::associated_patch;
```




<hr>



### typedef stored\_interface 

_Tag the equivalent interface defined in the main file, i.e. given to_ [_**MultipatchConnectivity**_](classMultipatchConnectivity.md) _._
```C++
using DerivativeCalculatorCommon< Connectivity, Grid1D >::stored_interface =  find_associated_interface_t<typename interface<I>::Edge1, all_interface_collection>;
```




<hr>
## Protected Attributes Documentation




### variable m\_idx\_ranges 

_MultipatchType on IdxRangeOnPatch to collect the index ranges on the patches._ 
```C++
MultipatchIdxRange DerivativeCalculatorCommon< Connectivity, Grid1D >::m_idx_ranges;
```




<hr>
## Protected Static Attributes Documentation




### variable extremity\_1 

_Tag the connection extremity of the first patch of the interface._ 
```C++
Extremity DerivativeCalculatorCommon< Connectivity, Grid1D >::extremity_1;
```




<hr>



### variable extremity\_2 

_Tag the connection extremity of the second patch of the interface._ 
```C++
Extremity DerivativeCalculatorCommon< Connectivity, Grid1D >::extremity_2;
```




<hr>



### variable is\_periodic 

_Boolean to indicate the periodicty of the global domain._ 
```C++
bool DerivativeCalculatorCommon< Connectivity, Grid1D >::is_periodic;
```




<hr>



### variable n\_inner\_interfaces 

_Number of inner interfaces (fixing the number of elements in the matrix and vectors)._ 
```C++
std::size_t DerivativeCalculatorCommon< Connectivity, Grid1D >::n_inner_interfaces;
```




<hr>



### variable patch\_order\_matches\_interfaces 

_Boolean if the order of Patch1\|Patch2 of the interface matches with the sorted 1D grid sequence._ 
```C++
bool DerivativeCalculatorCommon< Connectivity, Grid1D >::patch_order_matches_interfaces;
```




<hr>
## Protected Functions Documentation




### function get\_approximated\_deriv [1/2]

_Compute the approximated derivative (or coefficient "c") at the interface along the mesh line indexed at the given indices._ 
```C++
template<std::size_t I, bool is_first_deriv, std::enable_if_t< is_first_deriv, bool >>
inline double DerivativeCalculatorCommon::get_approximated_deriv (
    Idx< grid_par1 < I > > const & idx_to_slice_1,
    Idx< grid_par2 < I > > const & idx_to_slice_2,
    MultipatchDerivConstField_host functions_and_derivs
) const
```





**Template parameters:**


* `I` Index of the interface. 



**Parameters:**


* `idx_to_slice_1` Index on the parallel grid of the first patch to select the mesh line along which the derivative is computed. 
* `idx_to_slice_2` Index on the parallel grid of the second patch to select the mesh line along which the derivative is computed. 
* `functions_and_derivs` Collection fields with its first derivatives and cross-derivatives. 



**Returns:**

The approximated derivative at the interface. 





        

<hr>



### function get\_approximated\_deriv [1/2]

_Compute the approximated cross-derivative (or coefficient "c") at the interface along the mesh line indexed at the given indices. Work only for indices in the corners of the patches (where the cross-derivatives are supposed to be computed)._ 
```C++
template<std::size_t I, bool is_first_deriv, std::enable_if_t<!is_first_deriv, bool >>
inline double DerivativeCalculatorCommon::get_approximated_deriv (
    Idx< grid_par1 < I > > const & idx_to_slice_1,
    Idx< grid_par2 < I > > const & idx_to_slice_2,
    MultipatchDerivConstField_host functions_and_derivs
) const
```





**Template parameters:**


* `I` Index of the interface. 



**Parameters:**


* `idx_to_slice_1` Index on the parallel grid of the first patch to select the mesh line along which the derivative is computed. 
* `idx_to_slice_2` Index on the parallel grid of the second patch to select the mesh line along which the derivative is computed. 
* `functions_and_derivs` Collection fields with its first derivatives and cross-derivatives. 



**Returns:**

The approximated cross-derivative at the interface. 





        

<hr>



### function get\_deriv\_coeffs 

_Get the_ [_**InterfaceDerivCoeffs**_](classInterfaceDerivCoeffs.md) _defined on the I-th inner interface._
```C++
template<std::size_t I>
inline InterfaceDerivCoeffs < stored_interface < I > > const & DerivativeCalculatorCommon::get_deriv_coeffs () const
```





**Template parameters:**


* `I` Index of the interface in [**inner\_interface\_collection**](classDerivativeCalculatorCommon.md#typedef-inner_interface_collection). 



**Returns:**

A constant reference to the [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md). 





        

<hr>



### function get\_idx\_interface 

_Get the index of the given index range slice corresponding to one of the extremities of the 1D grid. If is\_on\_interface is true, the extremity corresponds to the interface. If is\_on\_interface is false, the extremity corresponds to the other bound of the interface (next or previous)._ 
```C++
template<bool is_on_interface, std::size_t I, class Patch, bool is_patch1>
inline std::conditional_t< is_patch1, Idx< grid_perp1 < I > >, Idx< grid_perp2 < I > > > DerivativeCalculatorCommon::get_idx_interface () const
```




<hr>



### function set\_slice\_indexes 

_Get the indices on Patch1 and Patch2 of the given interface. The index on Patch1 corresponds to the value given in input (if the interface is well oriented, otherwise it is Patch2). The index on Patch2 is the equivalent index of the one on Patch1. It also updates the value given in input with the value of the index on Patch2 for the next interface (if the next interface is well oriented, otherwise it is Patch1)._ 
```C++
template<std::size_t I>
inline void DerivativeCalculatorCommon::set_slice_indexes (
    int & value_idx_to_slice,
    Idx< grid_par1 < I > > & idx_to_slice_1,
    Idx< grid_par2 < I > > & idx_to_slice_2
) const
```



Example for well oriented interfaces:


\| Patch1 \| Patch2 Patch1 \| Patch2 \| 
 \| idx\_1 &gt; idx\_2 idx\_1 &gt; idx\_2 \| \| ^ \| v ^ \| v \| \|idx\_val idx\_val \| idx\_val idx\_val \| idx\_val \|




**Parameters:**


* `value_idx_to_slice` Integer containing the value of the index on the parallel grid of the last interface. 
* `idx_to_slice_1` Index on the parallel grid of the first patch to select the mesh line along which the derivative will be computed. 
* `idx_to_slice_2` Index on the parallel grid of the second patch to select the mesh line along which the derivative will be computed. 




        

<hr>
## Protected Static Functions Documentation




### function to\_sign 

_Transform a boolean into a integer. If true, the sign is positive. If false, the sign is negative._ 
```C++
static inline int DerivativeCalculatorCommon::to_sign (
    bool agree
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/derivative_calculator_common.hpp`

