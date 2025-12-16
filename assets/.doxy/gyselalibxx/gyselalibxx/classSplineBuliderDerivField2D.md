

# Class SplineBuliderDerivField2D

**template &lt;class ExecSpace, class BSplines1, class BSplines2, class Grid1, class Grid2, ddc::BoundCond BoundCond1min, ddc::BoundCond BoundCond1max, ddc::BoundCond BoundCond2min, ddc::BoundCond BoundCond2max&gt;**



[**ClassList**](annotated.md) **>** [**SplineBuliderDerivField2D**](classSplineBuliderDerivField2D.md)



_[Temporary] Apply a SplineBuilder2D to a_ [_**DerivField**_](classDerivField.md) _._[More...](#detailed-description)

* `#include <spline_builder_deriv_field_2d.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplineBuliderDerivField2D**](#function-splinebuliderderivfield2d) (Builder2D const & builder) <br>_Instantiate the class by storing a reference to the spline builder we can to use._  |
|  void | [**fill\_in\_cross\_deriv**](#function-fill_in_cross_deriv) (CrossDerivField cross\_deriv, [**DerivFieldType**](classDerivField.md) function\_and\_derivs, Idx&lt; Grid1 &gt; idx\_slice\_1, Idx&lt; Grid2 &gt; idx\_slice\_2) const<br>_Fill in the cross\_deriv field with the cross derivatives stored in the function\_and\_derivs._  |
|  void | [**fill\_in\_deriv1**](#function-fill_in_deriv1) (Deriv1Field deriv1, [**DerivFieldType**](classDerivField.md) function\_and\_derivs, Idx&lt; Grid1 &gt; idx\_slice) const<br>_Fill in the deriv1 field with the derivatives along the first dimension stored in the function\_and\_derivs._  |
|  void | [**fill\_in\_deriv2**](#function-fill_in_deriv2) (Deriv2Field deriv2, [**DerivFieldType**](classDerivField.md) function\_and\_derivs, Idx&lt; Grid2 &gt; idx\_slice) const<br>_Fill in the deriv2 field with the derivatives along the second dimension stored in the function\_and\_derivs._  |
|  void | [**fill\_in\_function**](#function-fill_in_function) (FunctField function, [**DerivFieldType**](classDerivField.md) function\_and\_derivs) const<br>_Fill in the function field with the values stored in the function\_and\_derivs._  |
|  void | [**operator()**](#function-operator) (SplineType spline, [**DerivFieldType**](classDerivField.md) function\_and\_derivs) const<br>_Build the spline representation of the function given in a_ [_**DerivField**_](classDerivField.md) _applying the referenced spline builder stored in the class._ |




























## Detailed Description


[**DerivField**](classDerivField.md) stores the values and the derivatives of a function. The inputs of the SplineBuilder2D are on different layouts than the fields we can get from the [**DerivField**](classDerivField.md). [**SplineBuliderDerivField2D**](classSplineBuliderDerivField2D.md) allows to directly apply a stored SplineBuilder2D to a [**DerivField**](classDerivField.md) by copying data in fields with correct layout.


Implemented only for 2D case. 


    
## Public Functions Documentation




### function SplineBuliderDerivField2D 

_Instantiate the class by storing a reference to the spline builder we can to use._ 
```C++
inline explicit SplineBuliderDerivField2D::SplineBuliderDerivField2D (
    Builder2D const & builder
) 
```





**Parameters:**


* `builder` A reference to SplineBuilder2D from DDC that we store in the class. 




        

<hr>



### function fill\_in\_cross\_deriv 

_Fill in the cross\_deriv field with the cross derivatives stored in the function\_and\_derivs._ 
```C++
inline void SplineBuliderDerivField2D::fill_in_cross_deriv (
    CrossDerivField cross_deriv,
    DerivFieldType function_and_derivs,
    Idx< Grid1 > idx_slice_1,
    Idx< Grid2 > idx_slice_2
) const
```





**Parameters:**


* `cross_deriv` Field with layout\_right where we copy the cross-derivatives. 
* `function_and_derivs` [**DerivField**](classDerivField.md) from where the cross-derivatives are copied. 
* `idx_slice_1` Index to determine which bound (mon/max) we select for the derivative field on the first dimension. 
* `idx_slice_2` Index to determine which bound (mon/max) we select for the derivative field on the second dimension. 




        

<hr>



### function fill\_in\_deriv1 

_Fill in the deriv1 field with the derivatives along the first dimension stored in the function\_and\_derivs._ 
```C++
inline void SplineBuliderDerivField2D::fill_in_deriv1 (
    Deriv1Field deriv1,
    DerivFieldType function_and_derivs,
    Idx< Grid1 > idx_slice
) const
```





**Parameters:**


* `deriv1` Field with layout\_right where we copy the derivatives. 
* `function_and_derivs` [**DerivField**](classDerivField.md) from where the derivatives are copied. 
* `idx_slice` Index to determine which bound (mon/max) we select for the derivative field. 




        

<hr>



### function fill\_in\_deriv2 

_Fill in the deriv2 field with the derivatives along the second dimension stored in the function\_and\_derivs._ 
```C++
inline void SplineBuliderDerivField2D::fill_in_deriv2 (
    Deriv2Field deriv2,
    DerivFieldType function_and_derivs,
    Idx< Grid2 > idx_slice
) const
```





**Parameters:**


* `deriv2` Field with layout\_right where we copy the derivatives. 
* `function_and_derivs` [**DerivField**](classDerivField.md) from where the derivatives are copied. 
* `idx_slice` Index to determine which bound (mon/max) we select for the derivative field. 




        

<hr>



### function fill\_in\_function 

_Fill in the function field with the values stored in the function\_and\_derivs._ 
```C++
inline void SplineBuliderDerivField2D::fill_in_function (
    FunctField function,
    DerivFieldType function_and_derivs
) const
```





**Parameters:**


* `function` Field with layout\_right where we copy the function values. 
* `function_and_derivs` [**DerivField**](classDerivField.md) from where the function values are copied. 




        

<hr>



### function operator() 

_Build the spline representation of the function given in a_ [_**DerivField**_](classDerivField.md) _applying the referenced spline builder stored in the class._
```C++
inline void SplineBuliderDerivField2D::operator() (
    SplineType spline,
    DerivFieldType function_and_derivs
) const
```





**Parameters:**


* `spline` Spline coefficients on a 2D grid. 
* `function_and_derivs` Data defining the function on a 2D grid. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/spline_builder_deriv_field_2d.hpp`

