

# Class SingleInterfaceDerivativesCalculator

**template &lt;class InterfaceType&gt;**



[**ClassList**](annotated.md) **>** [**SingleInterfaceDerivativesCalculator**](classSingleInterfaceDerivativesCalculator.md)



_Compute the derivative of an equivalent global spline at the interface between two patches._ [More...](#detailed-description)

* `#include <single_interface_derivatives_calculator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef InterfaceType | [**associated\_interface**](#typedef-associated_interface)  <br>[_**Interface**_](structInterface.md) _between the two involved patches._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SingleInterfaceDerivativesCalculator**](#function-singleinterfacederivativescalculator-14) (IdxRange1DPerp\_1 const & idx\_range\_1d\_1, IdxRange1DPerp\_2 const & idx\_range\_1d\_2, ddc::BoundCond const & Bound1=ddc::BoundCond::HERMITE, ddc::BoundCond const & Bound2=ddc::BoundCond::HERMITE) <br>_Instantiate_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _._ |
|   | [**SingleInterfaceDerivativesCalculator**](#function-singleinterfacederivativescalculator-24) (IdxRangeA const & idx\_range\_a, IdxRangeB const & idx\_range\_b, ddc::BoundCond const & Bound1=ddc::BoundCond::HERMITE, ddc::BoundCond const & Bound2=ddc::BoundCond::HERMITE) <br>_Instantiate_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _. See_SingleInterfaceDerivativesCalculatorInstantiator _._ |
|   | [**SingleInterfaceDerivativesCalculator**](#function-singleinterfacederivativescalculator-34) (IdxRange1DPerp\_1 const & idx\_range\_1d\_1, IdxRange1DPerp\_2 const & idx\_range\_1d\_2, std::size\_t const number\_chosen\_cells) <br>_Instantiate_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _. See_SingleInterfaceDerivativesCalculatorInstantiator _. This constructor calculates an approximation of the formula._ |
|   | [**SingleInterfaceDerivativesCalculator**](#function-singleinterfacederivativescalculator-44) (IdxRangeA const & idx\_range\_a, IdxRangeB const & idx\_range\_b, std::size\_t const number\_chosen\_cells) <br>_Instantiate_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _. See_SingleInterfaceDerivativesCalculatorInstantiator _. This constructor calculates an approximation of the formula._ |
|  double | [**get\_coeff\_deriv\_on\_patch**](#function-get_coeff_deriv_on_patch) () const<br>_Get the coefficient (a) or (b) in front of the derivative on the given patch of the given_ [_**Interface**_](structInterface.md) _._ |
|  double | [**get\_coeff\_deriv\_patch\_1**](#function-get_coeff_deriv_patch_1) () const<br>_Get the coefficient (a) in front of the derivative on the patch 1 of the given_ [_**Interface**_](structInterface.md) _._ |
|  double | [**get\_coeff\_deriv\_patch\_2**](#function-get_coeff_deriv_patch_2) () const<br>_Get the coefficient (b) in front of the derivative on the patch 2 of the given_ [_**Interface**_](structInterface.md) _._ |
|  double | [**get\_function\_coefficients**](#function-get_function_coefficients-12) (DConstField&lt; IdxRange1DPerp\_1, Kokkos::HostSpace, Layout1 &gt; const & function\_1, DConstField&lt; IdxRange1DPerp\_2, Kokkos::HostSpace, Layout2 &gt; const & function\_2) const<br>_Get the linear combination of the function values (c)._  |
|  double | [**get\_function\_coefficients**](#function-get_function_coefficients-22) (DConstField&lt; IdxRange1DPerp\_2, Kokkos::HostSpace, Layout2 &gt; const & function\_2, DConstField&lt; IdxRange1DPerp\_1, Kokkos::HostSpace, Layout1 &gt; const & function\_1) const<br>_Get the linear combination of the function values (c). See_ get\_function\_coefficients _._ |




























## Detailed Description


For a given [**Interface**](structInterface.md), this operator computes the coefficients a, b and c of the following relation: \(s'(X_I) = c + a s'(X_{I+1}) + b s'(X_{I+1})\),


with
* s'(X\_I) the derivative at the [**Interface**](structInterface.md) of the two local splines on the patches.
* s'(X\_{I+1}) the derivative at the next (or right) [**Interface**](structInterface.md) of the local splines on patch 1.
* s'(X\_{I-1}) the derivative at the previous (or left) [**Interface**](structInterface.md) of the local splines on patch 2.
* a,b scalars depending on the meshes of the two patches.
* c a linear combination of the function values on the two patches. The weights in the linear combination depend on the meshes of the two patches. It can be written as a sum: \(c = \sum_{k = -N2}^{N1} \omega_k f_k\), with \(\omega_k\) the weights, \(f_k\) the function values, \(N1\) the number of cells in patch 1, and \(N2\) the number of cells in patch 2.




Scheme of the two patches: X\_{I-1} X\_I X\_{I+1} 
 \| \| \| \| 1 \| 2 \| \|\_\_\_\_\_\_\_\|\_\_\_\_\_\_\_\|
* ← → +
* → ← - Orientation of an equivalent global mesh : →




All the formulae and more details are given in the README.md.




**Template parameters:**


* [**Interface**](structInterface.md) The interface between two patches where we want to compute the derivatives. 




**Warning:**

The applied method only works for interpolation points located on the break points. In the case where "ddc::BoundCond::GREVILLE" is specified, additional interpolation points are also placed in the first or last cell of the patch. Please be sure to initialise the discrete space of your Grid on the break points (especially in the non-uniform case). 





    
## Public Types Documentation




### typedef associated\_interface 

[_**Interface**_](structInterface.md) _between the two involved patches._
```C++
using SingleInterfaceDerivativesCalculator< InterfaceType >::associated_interface =  InterfaceType;
```




<hr>
## Public Functions Documentation




### function SingleInterfaceDerivativesCalculator [1/4]

_Instantiate_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _._
```C++
inline SingleInterfaceDerivativesCalculator::SingleInterfaceDerivativesCalculator (
    IdxRange1DPerp_1 const & idx_range_1d_1,
    IdxRange1DPerp_2 const & idx_range_1d_2,
    ddc::BoundCond const & Bound1=ddc::BoundCond::HERMITE,
    ddc::BoundCond const & Bound2=ddc::BoundCond::HERMITE
) 
```



 It computes the coefficients a and b, and the weights \(\omega_k\) and stores the values in the class. If you want to compute the derivatives with an exact formula, please provide the index ranges on the full domain of the patches. If you want to use an approximation (e.g. use only 5 cells), please provide the index ranges on a reduced number of cells (e.g. 5 cells).


If the interpolation points are uniform, it computes the coefficients with an explicit formula. Otherwise, it uses the recursive formula (see README.md).




**Parameters:**


* `idx_range_1d_1` 1D index range perpendicular to the [**Interface**](structInterface.md), on the patch 1. 
* `idx_range_1d_2` 1D index range perpendicular to the [**Interface**](structInterface.md), on the patch 2. 
* `Bound1` The boundary condition type on the opposite edge of the interface on the patch 1. By default, the value is set to ddc::BoundCond::HERMITE. If ddc::BoundCond::GREVILLE is given, a treatment will be applied to consider the additional interpolation point. Giving ddc::BoundCond::PERIODIC does not make sense. 
* `Bound2` The boundary condition type on the opposite edge of the interface on the patch 2. By default, the value is set to ddc::BoundCond::HERMITE. If ddc::BoundCond::GREVILLE is given, a treatment will be applied to consider the additional interpolation point. Giving ddc::BoundCond::PERIODIC does not make sense. 
 




        

<hr>



### function SingleInterfaceDerivativesCalculator [2/4]

_Instantiate_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _. See_SingleInterfaceDerivativesCalculatorInstantiator _._
```C++
template<class IdxRangeA, class IdxRangeB>
inline SingleInterfaceDerivativesCalculator::SingleInterfaceDerivativesCalculator (
    IdxRangeA const & idx_range_a,
    IdxRangeB const & idx_range_b,
    ddc::BoundCond const & Bound1=ddc::BoundCond::HERMITE,
    ddc::BoundCond const & Bound2=ddc::BoundCond::HERMITE
) 
```





**Parameters:**


* `idx_range_a` Index range on one patch. 
* `idx_range_b` Index range on the other patch. 
* `Bound1` The boundary condition type on the opposite edge of the interface on the patch 1. By default, the value is set to ddc::BoundCond::HERMITE. 
* `Bound2` The boundary condition type on the opposite edge of the interface on the patch 2. By default, the value is set to ddc::BoundCond::HERMITE. 




        

<hr>



### function SingleInterfaceDerivativesCalculator [3/4]

_Instantiate_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _. See_SingleInterfaceDerivativesCalculatorInstantiator _. This constructor calculates an approximation of the formula._
```C++
inline SingleInterfaceDerivativesCalculator::SingleInterfaceDerivativesCalculator (
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



### function SingleInterfaceDerivativesCalculator [4/4]

_Instantiate_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _. See_SingleInterfaceDerivativesCalculatorInstantiator _. This constructor calculates an approximation of the formula._
```C++
template<class IdxRangeA, class IdxRangeB>
inline SingleInterfaceDerivativesCalculator::SingleInterfaceDerivativesCalculator (
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



### function get\_coeff\_deriv\_on\_patch 

_Get the coefficient (a) or (b) in front of the derivative on the given patch of the given_ [_**Interface**_](structInterface.md) _._
```C++
template<class Patch>
inline double SingleInterfaceDerivativesCalculator::get_coeff_deriv_on_patch () const
```





**Returns:**

The value of the coefficient (a) or (b). 





        

<hr>



### function get\_coeff\_deriv\_patch\_1 

_Get the coefficient (a) in front of the derivative on the patch 1 of the given_ [_**Interface**_](structInterface.md) _._
```C++
inline double SingleInterfaceDerivativesCalculator::get_coeff_deriv_patch_1 () const
```





**Returns:**

The value of the coefficient a. 





        

<hr>



### function get\_coeff\_deriv\_patch\_2 

_Get the coefficient (b) in front of the derivative on the patch 2 of the given_ [_**Interface**_](structInterface.md) _._
```C++
inline double SingleInterfaceDerivativesCalculator::get_coeff_deriv_patch_2 () const
```





**Returns:**

The value of the coefficient b. 





        

<hr>



### function get\_function\_coefficients [1/2]

_Get the linear combination of the function values (c)._ 
```C++
template<class Layout1, class Layout2>
inline double SingleInterfaceDerivativesCalculator::get_function_coefficients (
    DConstField< IdxRange1DPerp_1, Kokkos::HostSpace, Layout1 > const & function_1,
    DConstField< IdxRange1DPerp_2, Kokkos::HostSpace, Layout2 > const & function_2
) const
```



 

**Parameters:**


* `function_1` Function values at the interpolation points on patch 1. 
* `function_2` Function values at the interpolation points on patch 2. 



**Returns:**

the linear combination of the function values (c).


Remark: for the approximation case, the cells were selected by the given index ranges in the instantiation. Here, function\_1 and function\_2 can be given on the whole domain and the operator will select the correct values. But, it could be more optimised to directly give a slice on the selected cells. E.g. get\_function\_coefficients(function\_1[idx\_interface\_1][selected\_cells\_idx\_range\_1], 
                          function\_2[idx\_interface\_2][selected\_cells\_idx\_range\_2]); and get\_function\_coefficients(function\_1[idx\_interface\_1], 
                          function\_2[idx\_interface\_2]); will return the same value. 


        

<hr>



### function get\_function\_coefficients [2/2]

_Get the linear combination of the function values (c). See_ get\_function\_coefficients _._
```C++
template<class Layout1, class Layout2>
inline double SingleInterfaceDerivativesCalculator::get_function_coefficients (
    DConstField< IdxRange1DPerp_2, Kokkos::HostSpace, Layout2 > const & function_2,
    DConstField< IdxRange1DPerp_1, Kokkos::HostSpace, Layout1 > const & function_1
) const
```





**Parameters:**


* `function_1` Function values at the interpolation points on patch 1. 
* `function_2` Function values at the interpolation points on patch 2. 



**Returns:**

the linear combination of the function values (c). 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/single_interface_derivatives_calculator.hpp`

