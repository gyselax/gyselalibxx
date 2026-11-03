

# Class InterfacesDerivativeCalculator&lt; Connectivity, Grid1D, ddc::TypeSeq&lt; Patches... &gt;, DerivCoeffCollection &gt;

**template &lt;class Connectivity, class Grid1D, class... Patches, class DerivCoeffCollection&gt;**



[**ClassList**](annotated.md) **>** [**InterfacesDerivativeCalculator&lt; Connectivity, Grid1D, ddc::TypeSeq&lt; Patches... &gt;, DerivCoeffCollection &gt;**](classInterfacesDerivativeCalculator_3_01Connectivity_00_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patches992a465ddaf32ad2e1b39ede53a467ad.md)



_Class to compute the interface derivatives along a given direction on all the interfaces with an approximation formula._ [More...](#detailed-description)

* `#include <interfaces_derivative_calculator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef PatchFirst | [**first\_patch**](#typedef-first_patch)  <br>_First patch of the collection of patches in the direction of the given Grid1D._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**InterfacesDerivativeCalculator**](#function-interfacesderivativecalculator) (MultipatchType&lt; IdxRangeOnPatch, Patches... &gt; const & idx\_ranges, DerivCoeffCollection const & deriv\_coeffs\_collect) <br>_Instantiate_ [_**InterfacesDerivativeCalculator**_](classInterfacesDerivativeCalculator.md) _._ |
|  void | [**solve\_cross\_deriv**](#function-solve_cross_deriv-12) (MultipatchField&lt; DerivFieldOnPatch\_host, Patches... &gt; functions\_and\_derivs, IdxPar const & idx\_par) <br>_Compute the vector C\_trunc to determine all the interface cross-derivatives in the Grid1D direction at a given index._  |
|  void | [**solve\_cross\_deriv**](#function-solve_cross_deriv-22) (MultipatchField&lt; DerivFieldOnPatch\_host, Patches... &gt; functions\_and\_derivs) <br>_Compute the vector C\_trunc to determine all the corner cross-derivatives in the Grid1D direction at all the last and the first indices of the first patch._  |
|  void | [**solve\_deriv**](#function-solve_deriv-12) (MultipatchField&lt; DerivFieldOnPatch\_host, Patches... &gt; functions\_and\_derivs, IdxPar const & idx\_par) <br>_Compute the vector C\_trunc to determine all the interface derivatives in the Grid1D direction at a given index._  |
|  void | [**solve\_deriv**](#function-solve_deriv-22) (MultipatchField&lt; DerivFieldOnPatch\_host, Patches... &gt; functions\_and\_derivs) <br>_Compute the vector C\_trunc to determine all the interface derivatives in the Grid1D direction at all the indices of the first patch._  |




























## Detailed Description


This operator is implemented for conforming equivalent global meshes.


When we call the operator .[**solve\_deriv()**](classInterfacesDerivativeCalculator_3_01Connectivity_00_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patches992a465ddaf32ad2e1b39ede53a467ad.md#function-solve_deriv-12) or .[**solve\_cross\_deriv()**](classInterfacesDerivativeCalculator_3_01Connectivity_00_01Grid1D_00_01ddc_1_1TypeSeq_3_01Patches992a465ddaf32ad2e1b39ede53a467ad.md#function-solve_cross_deriv-12), the operator loops over the interface to compute the vector C\_trunc to determine all the interface derivatives and update the [**DerivField**](classDerivField.md) given in input.


The vector C\_trunc is given by the coefficients (see README), 
\[\{\sum_{k = - N_{reduc}}^{N_{reduc}} \omega_{k, N_{reduc}, N_{reduc}}^{i_I} f_{i_I+k} \}_I\]



with
* \(N_{reduc}\) the number of chosen cells for the approximation,
* \(i_I\) the index of the grid at the Ith interface,
* \(f_{i_I+k}\) the function values if we compute the first derivatives, or the first derivatives if we compute the cross-derivatives.






**Template parameters:**


* `Connectivity` A [**MultipatchConnectivity**](classMultipatchConnectivity.md) class describing all the patch connections. 
* `Grid1D` A given direction. 
* `Patches` List of patches containing all the involved patches in the given direction. 
* `DerivCoeffCollection` A [**InterfaceDerivCoeffsCollection**](classInterfaceDerivCoeffsCollection.md) that stores the [**InterfaceDerivCoeffs**](classInterfaceDerivCoeffs.md) needed to compute the interface derivatives. 




    
## Public Types Documentation




### typedef first\_patch 

_First patch of the collection of patches in the direction of the given Grid1D._ 
```C++
using InterfacesDerivativeCalculator< Connectivity, Grid1D, ddc::TypeSeq< Patches... >, DerivCoeffCollection >::first_patch =  PatchFirst;
```




<hr>
## Public Functions Documentation




### function InterfacesDerivativeCalculator 

_Instantiate_ [_**InterfacesDerivativeCalculator**_](classInterfacesDerivativeCalculator.md) _._
```C++
inline InterfacesDerivativeCalculator< Connectivity, Grid1D, ddc::TypeSeq< Patches... >, DerivCoeffCollection >::InterfacesDerivativeCalculator (
    MultipatchType< IdxRangeOnPatch, Patches... > const & idx_ranges,
    DerivCoeffCollection const & deriv_coeffs_collect
) 
```





**Parameters:**


* `idx_ranges` MultipatchType collection of index ranges defined on the given list of patches. 
* `deriv_coeffs_collect` [**InterfaceDerivCoeffsCollection**](classInterfaceDerivCoeffsCollection.md) containing all the interface derivative calculator for each interface in the given Grid1D direction. 




        

<hr>



### function solve\_cross\_deriv [1/2]

_Compute the vector C\_trunc to determine all the interface cross-derivatives in the Grid1D direction at a given index._ 
```C++
template<class IdxPar>
inline void InterfacesDerivativeCalculator< Connectivity, Grid1D, ddc::TypeSeq< Patches... >, DerivCoeffCollection >::solve_cross_deriv (
    MultipatchField< DerivFieldOnPatch_host, Patches... > functions_and_derivs,
    IdxPar const & idx_par
) 
```



It uses the first-derivatives s to compute the cross-derivatives in a perpendicular direction to the interfaces.




**Template parameters:**


* `IdxPar` Type for the given index. It can be on the first or the second grid of the first patch. 



**Parameters:**


* `functions_and_derivs` Collection fields with its first derivatives and cross-derivatives. 
* `idx_par` Index of the line in the Grid1D direction where we compute all the interface cross-derivatives. 




        

<hr>



### function solve\_cross\_deriv [2/2]

_Compute the vector C\_trunc to determine all the corner cross-derivatives in the Grid1D direction at all the last and the first indices of the first patch._ 
```C++
inline void InterfacesDerivativeCalculator< Connectivity, Grid1D, ddc::TypeSeq< Patches... >, DerivCoeffCollection >::solve_cross_deriv (
    MultipatchField< DerivFieldOnPatch_host, Patches... > functions_and_derivs
) 
```



It uses the first-derivatives s to compute the cross-derivatives in a perpendicular direction to the interfaces.




**Parameters:**


* `functions_and_derivs` Collection fields with its first derivatives and cross-derivatives. 



**Warning:**

The global mesh is supposed to be conforming. 





        

<hr>



### function solve\_deriv [1/2]

_Compute the vector C\_trunc to determine all the interface derivatives in the Grid1D direction at a given index._ 
```C++
template<class IdxPar>
inline void InterfacesDerivativeCalculator< Connectivity, Grid1D, ddc::TypeSeq< Patches... >, DerivCoeffCollection >::solve_deriv (
    MultipatchField< DerivFieldOnPatch_host, Patches... > functions_and_derivs,
    IdxPar const & idx_par
) 
```



It uses the function values to compute the first-derivatives perpendicular to the interfaces.




**Template parameters:**


* `IdxPar` Type for the given index. It can be on the first or the second grid of the first patch. 



**Parameters:**


* `functions_and_derivs` Collection fields with its first derivatives and cross-derivatives. 
* `idx_par` Index of the line in the Grid1D direction where we compute all the interface derivatives. 




        

<hr>



### function solve\_deriv [2/2]

_Compute the vector C\_trunc to determine all the interface derivatives in the Grid1D direction at all the indices of the first patch._ 
```C++
inline void InterfacesDerivativeCalculator< Connectivity, Grid1D, ddc::TypeSeq< Patches... >, DerivCoeffCollection >::solve_deriv (
    MultipatchField< DerivFieldOnPatch_host, Patches... > functions_and_derivs
) 
```



It uses the function values to compute the first-derivatives perpendicular to the interfaces.




**Parameters:**


* `functions_and_derivs` Collection fields with its first derivatives and cross-derivatives we want to compute. 



**Warning:**

The global mesh is supposed to be conforming. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/interfaces_derivative_calculator.hpp`

