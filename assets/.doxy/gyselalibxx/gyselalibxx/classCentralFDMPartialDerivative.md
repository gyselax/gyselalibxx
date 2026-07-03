

# Class CentralFDMPartialDerivative

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**CentralFDMPartialDerivative**](classCentralFDMPartialDerivative.md)



_A class which implements a partial derivative operator using a finite differences calculation of order two. A decentered scheme is used at the boundary, whereas centred finite difference are used inside the domain._ [More...](#detailed-description)

* `#include <central_fdm_partial_derivatives.hpp>`



Inherits the following classes: [IPartialDerivative](classIPartialDerivative.md)
















## Public Types inherited from IPartialDerivative

See [IPartialDerivative](classIPartialDerivative.md)

| Type | Name |
| ---: | :--- |
| typedef DConstField&lt; IdxRangeFull &gt; | [**DConstFieldType**](classIPartialDerivative.md#typedef-dconstfieldtype)  <br>_The type of a constant reference to the field to be differentiated._  |
| typedef DField&lt; IdxRangeFull &gt; | [**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype)  <br>_The type of a reference to the field to be differentiated._  |
| typedef find\_grid\_t&lt; DerivativeDimension, ddc::to\_type\_seq\_t&lt; IdxRangeFull &gt; &gt; | [**GridDerivativeDimension**](classIPartialDerivative.md#typedef-gridderivativedimension)  <br>_The type of the grid on the dimension on which the partial derivative is calculated._  |
| typedef ddc::remove\_dims\_of\_t&lt; IdxRangeFull, [**GridDerivativeDimension**](classIPartialDerivative.md#typedef-gridderivativedimension) &gt; | [**IdxRangeBatch**](classIPartialDerivative.md#typedef-idxrangebatch)  <br>_The index range of all dimensions except DerivativeDimension._  |
| typedef IdxRange&lt; [**GridDerivativeDimension**](classIPartialDerivative.md#typedef-gridderivativedimension) &gt; | [**IdxRangeDeriv**](classIPartialDerivative.md#typedef-idxrangederiv)  <br>_The index range of the dimension on which the partial derivative is calculated._  |






































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CentralFDMPartialDerivative**](#function-centralfdmpartialderivative) (DConstFieldType const field\_ref) <br>_Construct an instance of the class_ [_**CentralFDMPartialDerivative**_](classCentralFDMPartialDerivative.md) _._ |
| virtual void | [**operator()**](#function-operator) (DFieldType differentiated\_field) const<br>_Compute the partial derivative of a field in a given direction using a finite difference scheme. For more information about the coefficients, see_ `./README.md` __ |


## Public Functions inherited from IPartialDerivative

See [IPartialDerivative](classIPartialDerivative.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIPartialDerivative.md#function-operator) ([**DFieldType**](classIPartialDerivative.md#typedef-dfieldtype) differentiated\_field) const = 0<br>_Compute the partial derivative of a field in a given direction._  |
| virtual  | [**~IPartialDerivative**](classIPartialDerivative.md#function-ipartialderivative) () = default<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION double | [**get\_derivative**](#function-get_derivative) (std::array&lt; Coord&lt; DerivativeDimension &gt;, 3 &gt; const & field\_positions, std::array&lt; double, 3 &gt; const & field\_values) <br>_Compute the partial derivative locally from the values at 3 known positions._  |
|  constexpr int | [**n\_local\_points**](#function-n_local_points) () <br>_The number of points required to calculate the derivative locally._  |




















































## Detailed Description




**Template parameters:**


* `IdxRangeFull` The index range of the field on which the operator acts (with all dimensions, batched and dimension of interest). 
* `DerivativeDimension` The dimension on which the partial derivative is calculated. 




    
## Public Functions Documentation




### function CentralFDMPartialDerivative 

_Construct an instance of the class_ [_**CentralFDMPartialDerivative**_](classCentralFDMPartialDerivative.md) _._
```C++
inline explicit CentralFDMPartialDerivative::CentralFDMPartialDerivative (
    DConstFieldType const field_ref
) 
```





**Parameters:**


* `field_ref` The field to be differentiated. 




        

<hr>



### function operator() 

_Compute the partial derivative of a field in a given direction using a finite difference scheme. For more information about the coefficients, see_ `./README.md` __
```C++
inline virtual void CentralFDMPartialDerivative::operator() (
    DFieldType differentiated_field
) const
```





**Parameters:**


* `differentiated_field` On output, contains values of the differentiated field. 




        
Implements [*IPartialDerivative::operator()*](classIPartialDerivative.md#function-operator)


<hr>
## Public Static Functions Documentation




### function get\_derivative 

_Compute the partial derivative locally from the values at 3 known positions._ 
```C++
static inline KOKKOS_INLINE_FUNCTION double CentralFDMPartialDerivative::get_derivative (
    std::array< Coord< DerivativeDimension >, 3 > const & field_positions,
    std::array< double, 3 > const & field_values
) 
```





**Parameters:**


* `field_positions` The positions where the field is provided, in ascending order. 
* `field_values` The values of the field at the provided positions. 



**Returns:**

The derivative at the central point. 





        

<hr>



### function n\_local\_points 

_The number of points required to calculate the derivative locally._ 
```C++
static inline constexpr int CentralFDMPartialDerivative::n_local_points () 
```





**Returns:**

The number of points. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/central_fdm_partial_derivatives.hpp`

