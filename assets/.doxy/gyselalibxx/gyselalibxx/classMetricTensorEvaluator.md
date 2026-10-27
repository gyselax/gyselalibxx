

# Class MetricTensorEvaluator

**template &lt;concepts::MappingWithJacobian Mapping, class PositionCoordinate&gt;**



[**ClassList**](annotated.md) **>** [**MetricTensorEvaluator**](classMetricTensorEvaluator.md)



_An operator for calculating the metric tensor._ [More...](#detailed-description)

* `#include <metric_tensor_evaluator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef DTensor&lt; Dims &gt; | [**ContravariantVectorType**](#typedef-contravariantvectortype)  <br>_The type of a contravariant vector associated with this mapping._  |
| typedef PositionCoordinate | [**CoordArg**](#typedef-coordarg)  <br>_The type of a coordinate associated with this mapping._  |
| typedef DTensor&lt; vector\_index\_set\_dual\_t&lt; Dims &gt; &gt; | [**CovariantVectorType**](#typedef-covariantvectortype)  <br>_The type of a covariant vector associated with this mapping._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_FUNCTION | [**MetricTensorEvaluator**](#function-metrictensorevaluator) (Mapping mapping) <br>_A constructor for the metric tensor operator._  |
|  KOKKOS\_FUNCTION DTensor&lt; Dims, Dims &gt; | [**inverse**](#function-inverse-12) ([**CoordArg**](classMetricTensorEvaluator.md#typedef-coordarg) const & coord) const<br>_Compute the inverse metric tensor associated with the mapping at a given position in space._  |
|  KOKKOS\_FUNCTION DTensor&lt; Dims, Dims &gt; | [**inverse**](#function-inverse-22) (IdxType const & idx) const<br>_Compute the inverse metric tensor associated with the mapping at the coordinate of an index._  |
|  KOKKOS\_FUNCTION DTensor&lt; Dims\_cov, Dims\_cov &gt; | [**operator()**](#function-operator) ([**CoordArg**](classMetricTensorEvaluator.md#typedef-coordarg) const & coord) const<br>_Compute the metric tensor associated with the mapping at a given position in space._  |
|  KOKKOS\_FUNCTION DTensor&lt; Dims\_cov, Dims\_cov &gt; | [**operator()**](#function-operator_1) (IdxType const & idx) const<br>_Compute the metric tensor associated with the mapping at the coordinate of an index._  |




























## Detailed Description




**Template parameters:**


* `Mapping` The mapping providing the Jacobian operator. 
* `PositionCoordinate` The coordinate type where the metric tensor can be evaluated. 




    
## Public Types Documentation




### typedef ContravariantVectorType 

_The type of a contravariant vector associated with this mapping._ 
```C++
using MetricTensorEvaluator< Mapping, PositionCoordinate >::ContravariantVectorType =  DTensor<Dims>;
```




<hr>



### typedef CoordArg 

_The type of a coordinate associated with this mapping._ 
```C++
using MetricTensorEvaluator< Mapping, PositionCoordinate >::CoordArg =  PositionCoordinate;
```




<hr>



### typedef CovariantVectorType 

_The type of a covariant vector associated with this mapping._ 
```C++
using MetricTensorEvaluator< Mapping, PositionCoordinate >::CovariantVectorType =  DTensor<vector_index_set_dual_t<Dims>>;
```




<hr>
## Public Functions Documentation




### function MetricTensorEvaluator 

_A constructor for the metric tensor operator._ 
```C++
inline explicit KOKKOS_FUNCTION MetricTensorEvaluator::MetricTensorEvaluator (
    Mapping mapping
) 
```





**Parameters:**


* `mapping` The mapping which can be used to calculate the Jacobian. 




        

<hr>



### function inverse [1/2]

_Compute the inverse metric tensor associated with the mapping at a given position in space._ 
```C++
inline KOKKOS_FUNCTION DTensor< Dims, Dims > MetricTensorEvaluator::inverse (
    CoordArg const & coord
) const
```





**Parameters:**


* `coord` The coordinate where we evaluate the metric tensor. 



**Returns:**

inverse\_metric\_tensor A DTensor object containing the value of the inverse of the metric tensor. 





        

<hr>



### function inverse [2/2]

_Compute the inverse metric tensor associated with the mapping at the coordinate of an index._ 
```C++
template<ddc::concepts::discrete_element IdxType>
inline KOKKOS_FUNCTION DTensor< Dims, Dims > MetricTensorEvaluator::inverse (
    IdxType const & idx
) const
```





**Parameters:**


* `idx` The index where we evaluate the metric tensor. 



**Returns:**

inverse\_metric\_tensor A DTensor object containing the value of the inverse of the metric tensor. 





        

<hr>



### function operator() 

_Compute the metric tensor associated with the mapping at a given position in space._ 
```C++
inline KOKKOS_FUNCTION DTensor< Dims_cov, Dims_cov > MetricTensorEvaluator::operator() (
    CoordArg const & coord
) const
```



The metric tensor matrix is defined as: \(G = (J_{\mathcal{F}})^T J_{\mathcal{F}}\). with \(J_{\mathcal{F}}\) the Jacobian matrix.




**Parameters:**


* `coord` The coordinate where we evaluate the metric tensor. 



**Returns:**

metric\_tensor A DTensor object containing the value of the metric tensor. 





        

<hr>



### function operator() 

_Compute the metric tensor associated with the mapping at the coordinate of an index._ 
```C++
template<ddc::concepts::discrete_element IdxType>
inline KOKKOS_FUNCTION DTensor< Dims_cov, Dims_cov > MetricTensorEvaluator::operator() (
    IdxType const & idx
) const
```





**Parameters:**


* `idx` The index where we evaluate the metric tensor. 



**Returns:**

metric\_tensor A DTensor object containing the value of the metric tensor. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/coord_transformations/metric_tensor_evaluator.hpp`

