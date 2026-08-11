

# Class DiscreteMapping

**template &lt;class StartCoord, class EndCoord, concepts::InterpolationEvaluator NDEvaluator&gt;**



[**ClassList**](annotated.md) **>** [**DiscreteMapping**](classDiscreteMapping.md)



_A mapping whose values are only known at the mesh points of a grid, evaluated elsewhere using an interpolating function._ [More...](#detailed-description)

* `#include <discrete_mapping.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::to\_type\_seq\_t&lt; StartCoord &gt; | [**ArgBasis**](#typedef-argbasis)  <br>_The type sequence of dimensions describing the argument of the function described by this mapping._  |
| typedef StartCoord | [**CoordArg**](#typedef-coordarg)  <br>_The type of the argument of the function described by this mapping._  |
| typedef StartCoord | [**CoordJacobian**](#typedef-coordjacobian)  <br>_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._  |
| typedef EndCoord | [**CoordResult**](#typedef-coordresult)  <br>_The type of the result of the function described by this mapping._  |
| typedef typename [**InterpolationEvaluatorTraits**](structInterpolationEvaluatorTraits.md)&lt; NDEvaluator &gt;::data\_type | [**DataType**](#typedef-datatype)  <br>_The data type of the coefficients and the values calculated by this mapping._  |
| typedef ddc::to\_type\_seq\_t&lt; EndCoord &gt; | [**ResultBasis**](#typedef-resultbasis)  <br>_The type sequence of dimensions describing the result of the function described by this mapping._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**DiscreteMapping**](#function-discretemapping-12) ([**CoeffField**](classVectorField.md) coeff\_representation, NDEvaluator const & evaluator) <br>_Instantiate a_ [_**DiscreteMapping**_](classDiscreteMapping.md) _from the coefficients of an interpolating function which approximates the mapping._ |
|  KOKKOS\_DEFAULTED\_FUNCTION | [**DiscreteMapping**](#function-discretemapping-22) ([**DiscreteMapping**](classDiscreteMapping.md) const & other) = default<br>_Copy-construct a_ [_**DiscreteMapping**_](classDiscreteMapping.md) _._ |
|  KOKKOS\_FUNCTION double | [**jacobian**](#function-jacobian) ([**CoordJacobian**](classDiscreteMapping.md#typedef-coordjacobian) const & coord) const<br>_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**jacobian\_component**](#function-jacobian_component) ([**CoordJacobian**](classDiscreteMapping.md#typedef-coordjacobian) coord) const<br>_Compute the (i,j) coefficient of the Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**Tensor**](classTensor.md)&lt; [**DataType**](classDiscreteMapping.md#typedef-datatype), [**ResultBasis**](classDiscreteMapping.md#typedef-resultbasis), get\_covariant\_dims\_t&lt; [**ArgBasis**](classDiscreteMapping.md#typedef-argbasis) &gt; &gt; | [**jacobian\_matrix**](#function-jacobian_matrix) ([**CoordJacobian**](classDiscreteMapping.md#typedef-coordjacobian) const & coord) const<br>_Compute full Jacobian matrix._  |
|  KOKKOS\_FUNCTION [**CoordResult**](classDiscreteMapping.md#typedef-coordresult) | [**operator()**](#function-operator) ([**CoordArg**](classDiscreteMapping.md#typedef-coordarg) const & coord) const<br>_Compute the target coordinates from the start coordinates._  |
|  void | [**operator()**](#function-operator_1) (ExecSpace exec\_space, Field&lt; [**CoordResult**](classDiscreteMapping.md#typedef-coordresult), IdxRange&lt; GridType... &gt; &gt; coords) <br>_Compute the physical coordinates at the points on the logical grid._  |




























## Detailed Description


A discrete mapping is a mapping whose values are known only at the mesh points of the grid. To interpolate the mapping, we use an interpolation on a basis. The [**DiscreteMapping**](classDiscreteMapping.md) is initialised from the coefficients in front of the basis functions which arise when we approximate the functions \(E(S)\) (with \(E\) the coordinates in the End vector space and S the coordinates in the Start vector space). Then to interpolate the mapping, we evaluate the decomposed functions using the chosen interpolating function (see DiscreteMapping::operator()).




**Template parameters:**


* `StartCoord` The type of the coordinates in the domain of definition of the mapping. 
* `EndCoord` The type of the coordinates in the image of the mapping. 
* `NDEvaluator` The type of the evaluator used to evaluate the interpolating function. 




    
## Public Types Documentation




### typedef ArgBasis 

_The type sequence of dimensions describing the argument of the function described by this mapping._ 
```C++
using DiscreteMapping< StartCoord, EndCoord, NDEvaluator >::ArgBasis =  ddc::to_type_seq_t<StartCoord>;
```




<hr>



### typedef CoordArg 

_The type of the argument of the function described by this mapping._ 
```C++
using DiscreteMapping< StartCoord, EndCoord, NDEvaluator >::CoordArg =  StartCoord;
```




<hr>



### typedef CoordJacobian 

_The type of the coordinate that can be used to evaluate the Jacobian of this mapping._ 
```C++
using DiscreteMapping< StartCoord, EndCoord, NDEvaluator >::CoordJacobian =  StartCoord;
```




<hr>



### typedef CoordResult 

_The type of the result of the function described by this mapping._ 
```C++
using DiscreteMapping< StartCoord, EndCoord, NDEvaluator >::CoordResult =  EndCoord;
```




<hr>



### typedef DataType 

_The data type of the coefficients and the values calculated by this mapping._ 
```C++
using DiscreteMapping< StartCoord, EndCoord, NDEvaluator >::DataType =  typename InterpolationEvaluatorTraits<NDEvaluator>::data_type;
```




<hr>



### typedef ResultBasis 

_The type sequence of dimensions describing the result of the function described by this mapping._ 
```C++
using DiscreteMapping< StartCoord, EndCoord, NDEvaluator >::ResultBasis =  ddc::to_type_seq_t<EndCoord>;
```




<hr>
## Public Functions Documentation




### function DiscreteMapping [1/2]

_Instantiate a_ [_**DiscreteMapping**_](classDiscreteMapping.md) _from the coefficients of an interpolating function which approximates the mapping._
```C++
inline DiscreteMapping::DiscreteMapping (
    CoeffField coeff_representation,
    NDEvaluator const & evaluator
) 
```



A discrete mapping is a mapping whose values are known only at the mesh points of the grid. To interpolate the mapping, we use an interpolation on a basis. The [**DiscreteMapping**](classDiscreteMapping.md) is initialised from the coefficients in front of the basis functions which arise when we approximate the functions \(E(S)\) (with \(E\) the coordinates in the End vector space and S the coordinates in the Start vector space). Then to interpolate the mapping, we will evaluate the decomposed functions on the chosen interpolating function (see DiscreteMapping::operator()).


Here, the evaluator is given as input.




**Parameters:**


* `coeff_representation` The coefficients of the interpolating function representing the coords in the target domain.
* `evaluator` The evaluator used to evaluate the mapping. 




        

<hr>



### function DiscreteMapping [2/2]

_Copy-construct a_ [_**DiscreteMapping**_](classDiscreteMapping.md) _._
```C++
KOKKOS_DEFAULTED_FUNCTION DiscreteMapping::DiscreteMapping (
    DiscreteMapping const & other
) = default
```





**Parameters:**


* `other` The [**DiscreteMapping**](classDiscreteMapping.md) being copied. 




        

<hr>



### function jacobian 

_Compute the Jacobian, the determinant of the Jacobian matrix of the mapping._ 
```C++
inline KOKKOS_FUNCTION double DiscreteMapping::jacobian (
    CoordJacobian const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian.



**Returns:**

A double with the value of the determinant of the Jacobian matrix. 





        

<hr>



### function jacobian\_component 

_Compute the (i,j) coefficient of the Jacobian matrix._ 
```C++
template<class IndexTag1, class IndexTag2>
inline KOKKOS_INLINE_FUNCTION double DiscreteMapping::jacobian_component (
    CoordJacobian coord
) const
```



For a mapping given by \(\mathcal{F} : {q_i}\mapsto {x_i}\), with \({q_i}\) the curvilinear coordinates and \({x_i}\) the Cartesian coordinates, the (i,j) coefficient of the Jacobian matrix is given by \(J^i_j\frac{\partial x_i}{\partial q_j}\).




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix.



**Returns:**

A double with the value of the (i,j) coefficient of the Jacobian matrix.




**See also:** SplineEvaluator2D 



        

<hr>



### function jacobian\_matrix 

_Compute full Jacobian matrix._ 
```C++
inline KOKKOS_FUNCTION Tensor < DataType , ResultBasis , get_covariant_dims_t< ArgBasis > > DiscreteMapping::jacobian_matrix (
    CoordJacobian const & coord
) const
```



For some computations, we need the complete Jacobian matrix or just the coefficients. The coefficients can be given independently with the function jacobian\_component.




**Parameters:**


* `coord` The coordinate where we evaluate the Jacobian matrix. 



**Returns:**

The Jacobian matrix. 





        

<hr>



### function operator() 

_Compute the target coordinates from the start coordinates._ 
```C++
inline KOKKOS_FUNCTION CoordResult DiscreteMapping::operator() (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinates in the start domain.



**Returns:**

The coordinates of the mapping in the end domain. 





        

<hr>



### function operator() 

_Compute the physical coordinates at the points on the logical grid._ 
```C++
template<class ExecSpace, class... GridType>
inline void DiscreteMapping::operator() (
    ExecSpace exec_space,
    Field< CoordResult , IdxRange< GridType... > > coords
) 
```



It evaluates the decomposed mapping on B-splines at the coordinate point with a SplineEvaluator2D.




**Parameters:**


* `exec_space` The execution space where the calculation should be carried out. 
* `coords` The coordinates of the mapping in the physical domain at the grid points.



**See also:** SplineEvaluator2D 



        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/discrete_mapping.hpp`

