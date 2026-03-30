

# Class NDLagrangeEvaluator

**template &lt;class HeadEvaluator, class... Evaluators1D&gt;**



[**ClassList**](annotated.md) **>** [**NDLagrangeEvaluator**](classNDLagrangeEvaluator.md)



_Evaluates an ND Lagrange polynomial via a tensor product of 1D evaluations._ [More...](#detailed-description)

* `#include <nd_lagrange_evaluator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::detail::convert\_type\_seq\_to\_discrete\_domain\_t&lt; ddc::type\_seq\_remove\_t&lt; ddc::to\_type\_seq\_t&lt; BatchedInterpolationIdxRange &gt;, ddc::to\_type\_seq\_t&lt; [**evaluation\_idx\_range\_type**](classNDLagrangeEvaluator.md#typedef-evaluation_idx_range_type) &gt; &gt; &gt; | [**batch\_idx\_range\_type**](#typedef-batch_idx_range_type)  <br>_The type of the batch index range (obtained by removing the dimensions of interest from the whole domain)._  |
| typedef ddc::detail::convert\_type\_seq\_to\_discrete\_domain\_t&lt; ddc::type\_seq\_replace\_t&lt; ddc::to\_type\_seq\_t&lt; BatchedInterpolationIdxRange &gt;, ddc::to\_type\_seq\_t&lt; [**evaluation\_idx\_range\_type**](classNDLagrangeEvaluator.md#typedef-evaluation_idx_range_type) &gt;, ddc::to\_type\_seq\_t&lt; [**coeff\_idx\_range\_type**](classNDLagrangeEvaluator.md#typedef-coeff_idx_range_type) &gt; &gt; &gt; | [**batched\_coeff\_idx\_range\_type**](#typedef-batched_coeff_idx_range_type)  <br>_The type of the ND index range on which the Lagrange coefficients are defined plus any batch dimensions._  |
| typedef BatchedInterpolationIdxRange | [**batched\_evaluation\_idx\_range\_type**](#typedef-batched_evaluation_idx_range_type)  <br>_The type of the whole index range representing evaluation points including any batch dimensions._  |
| typedef IdxRange&lt; HeadCoeffGrid, typename Evaluators1D::coeff\_grid\_type... &gt; | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The type of the ND index range on which the Lagrange coefficients are defined._  |
| typedef typename HeadEvaluator::data\_type | [**data\_type**](#typedef-data_type)  <br>_The data type._  |
| typedef ddc::detail::convert\_type\_seq\_to\_discrete\_domain\_t&lt; type\_seq\_cat\_t&lt; ddc::to\_type\_seq\_t&lt; IdxRangeHeadEval &gt;, ddc::to\_type\_seq\_t&lt; typename Evaluators1D::evaluation\_idx\_range\_type &gt;... &gt; &gt; | [**evaluation\_idx\_range\_type**](#typedef-evaluation_idx_range_type)  <br>_The type of the index range for the ND evaluation mesh used by this class._  |
| typedef typename HeadEvaluator::exec\_space | [**exec\_space**](#typedef-exec_space)  <br>_The type of the Kokkos execution space._  |
| typedef typename HeadEvaluator::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The type of the Kokkos memory space._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**NDLagrangeEvaluator**](#function-ndlagrangeevaluator-13) (HeadEvaluator const & head\_eval, Evaluators1D const &... tail\_evals) <br>_Construct from a sequence of 1D LagrangeEvaluators, one per dimension._  |
|   | [**NDLagrangeEvaluator**](#function-ndlagrangeevaluator-23) ([**NDLagrangeEvaluator**](classNDLagrangeEvaluator.md) const & x) = default<br>_Copy-construct._  |
|   | [**NDLagrangeEvaluator**](#function-ndlagrangeevaluator-33) ([**NDLagrangeEvaluator**](classNDLagrangeEvaluator.md) && x) = default<br>_Move-construct._  |
|  KOKKOS\_FUNCTION double | [**deriv**](#function-deriv-13) (IdxDerivDims const & deriv\_order, Coord&lt; CoordsDims... &gt; const & coord\_eval, ConstField&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), BatchedLagrangeIdxRange, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout &gt; const lagrange\_coef) const<br>_Differentiate 1D Lagrange function at a given coordinate._  |
|  void | [**deriv**](#function-deriv-23) (IdxDerivDims const & deriv\_order, Field&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), IdxRangeBatchedInterpolation, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; const lagrange\_eval, ConstField&lt; Coord&lt; CoordsDims... &gt;, IdxRangeBatchedInterpolation, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; const coords\_eval, ConstField&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), [**batched\_coeff\_idx\_range\_type**](classNDLagrangeEvaluator.md#typedef-batched_coeff_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt;, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout3 &gt; const lagrange\_coef) const<br>_Differentiate 1D spline function (described by its spline coefficients) on a mesh._  |
|  void | [**deriv**](#function-deriv-33) (IdxDerivDims const & deriv\_order, Field&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), IdxRangeBatchedInterpolation, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; const lagrange\_eval, ConstField&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), [**batched\_coeff\_idx\_range\_type**](classNDLagrangeEvaluator.md#typedef-batched_coeff_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt;, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; const lagrange\_coef) const<br>_Differentiate 1D Lagrange function on a mesh._  |
|  KOKKOS\_FUNCTION [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type) | [**operator()**](#function-operator) (Coord&lt; CoordsDims... &gt; const & coord, ConstField&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), NdCoeffIdxRange, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout &gt; const & coeff) const<br>_Evaluate the ND Lagrange polynomial at a single coordinate._  |
|  void | [**operator()**](#function-operator_1) (Field&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), IdxRangeBatched, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; lagrange\_eval, ConstField&lt; Coord&lt; CoordsDims... &gt;, IdxRangeBatched, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; coords\_eval, ConstField&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), [**batched\_coeff\_idx\_range\_type**](classNDLagrangeEvaluator.md#typedef-batched_coeff_idx_range_type)&lt; IdxRangeBatched &gt;, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout3 &gt; lagrange\_coef) const<br>_Evaluate the ND Lagrange polynomial on a mesh (with explicit coordinates)._  |
|  void | [**operator()**](#function-operator_2) (Field&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), IdxRangeBatched, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; lagrange\_eval, ConstField&lt; [**data\_type**](classNDLagrangeEvaluator.md#typedef-data_type), [**batched\_coeff\_idx\_range\_type**](classNDLagrangeEvaluator.md#typedef-batched_coeff_idx_range_type)&lt; IdxRangeBatched &gt;, [**memory\_space**](classNDLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; lagrange\_coef) const<br>_Evaluate the ND Lagrange polynomial on a mesh (coordinates from the grid)._  |
|  [**NDLagrangeEvaluator**](classNDLagrangeEvaluator.md) & | [**operator=**](#function-operator_3) ([**NDLagrangeEvaluator**](classNDLagrangeEvaluator.md) const & x) = default<br>_Copy-assign._  |
|  [**NDLagrangeEvaluator**](classNDLagrangeEvaluator.md) & | [**operator=**](#function-operator_4) ([**NDLagrangeEvaluator**](classNDLagrangeEvaluator.md) && x) = default<br>_Move-assign._  |
|   | [**~NDLagrangeEvaluator**](#function-ndlagrangeevaluator) () = default<br>_Destruct._  |




























## Detailed Description


The ND Lagrange polynomial is defined as:  \(\[
  L(x_1, \ldots, x_N)
  = \sum_{i_1=0}^{d_1}\dots\sum_{i_N=0}^{d_N} \prod_{j=1}^N f(x_1, \ldots, x_N) l_{j,i_j}(x_j)
\]\) The evaluation is recursive:  \(\[
  L(x_1, \ldots, x_N)
  = \sum_{i_1=0}^{d_1} l_{1,i_1}(x_1) L(x_2, \dots, x_N)
\]\)


The recursion terminates when only one evaluator remains: in that case the tail type is the [**LagrangeEvaluator**](classLagrangeEvaluator.md) itself and the 1D evaluation (including its configured extrapolation rules) is invoked directly.




**Template parameters:**


* `HeadEvaluator` The 1D [**LagrangeEvaluator**](classLagrangeEvaluator.md) for the first dimension. 
* `Evaluators1D` The 1D LagrangeEvaluators for the remaining dimensions. 




    
## Public Types Documentation




### typedef batch\_idx\_range\_type 

_The type of the batch index range (obtained by removing the dimensions of interest from the whole domain)._ 
```C++
using NDLagrangeEvaluator< HeadEvaluator, Evaluators1D >::batch_idx_range_type =  ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_remove_t< ddc::to_type_seq_t<BatchedInterpolationIdxRange>, ddc::to_type_seq_t<evaluation_idx_range_type> >>;
```





**Template parameters:**


* `The` batched index range on which the interpolation points are defined. 




        

<hr>



### typedef batched\_coeff\_idx\_range\_type 

_The type of the ND index range on which the Lagrange coefficients are defined plus any batch dimensions._ 
```C++
using NDLagrangeEvaluator< HeadEvaluator, Evaluators1D >::batched_coeff_idx_range_type =  ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t< ddc::to_type_seq_t<BatchedInterpolationIdxRange>, ddc::to_type_seq_t<evaluation_idx_range_type>, ddc::to_type_seq_t<coeff_idx_range_type> >>;
```





**Template parameters:**


* `The` batched index range on which the interpolation points are defined. 




        

<hr>



### typedef batched\_evaluation\_idx\_range\_type 

_The type of the whole index range representing evaluation points including any batch dimensions._ 
```C++
using NDLagrangeEvaluator< HeadEvaluator, Evaluators1D >::batched_evaluation_idx_range_type =  BatchedInterpolationIdxRange;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef coeff\_idx\_range\_type 

_The type of the ND index range on which the Lagrange coefficients are defined._ 
```C++
using NDLagrangeEvaluator< HeadEvaluator, Evaluators1D >::coeff_idx_range_type =  IdxRange<HeadCoeffGrid, typename Evaluators1D::coeff_grid_type...>;
```




<hr>



### typedef data\_type 

_The data type._ 
```C++
using NDLagrangeEvaluator< HeadEvaluator, Evaluators1D >::data_type =  typename HeadEvaluator::data_type;
```




<hr>



### typedef evaluation\_idx\_range\_type 

_The type of the index range for the ND evaluation mesh used by this class._ 
```C++
using NDLagrangeEvaluator< HeadEvaluator, Evaluators1D >::evaluation_idx_range_type =  ddc::detail::convert_type_seq_to_discrete_domain_t<type_seq_cat_t< ddc::to_type_seq_t<IdxRangeHeadEval>, ddc::to_type_seq_t<typename Evaluators1D::evaluation_idx_range_type>...> >;
```




<hr>



### typedef exec\_space 

_The type of the Kokkos execution space._ 
```C++
using NDLagrangeEvaluator< HeadEvaluator, Evaluators1D >::exec_space =  typename HeadEvaluator::exec_space;
```




<hr>



### typedef memory\_space 

_The type of the Kokkos memory space._ 
```C++
using NDLagrangeEvaluator< HeadEvaluator, Evaluators1D >::memory_space =  typename HeadEvaluator::memory_space;
```




<hr>
## Public Functions Documentation




### function NDLagrangeEvaluator [1/3]

_Construct from a sequence of 1D LagrangeEvaluators, one per dimension._ 
```C++
inline explicit NDLagrangeEvaluator::NDLagrangeEvaluator (
    HeadEvaluator const & head_eval,
    Evaluators1D const &... tail_evals
) 
```





**Parameters:**


* `head_eval` The evaluator for the first dimension. 
* `tail_evals` The evaluators for the remaining dimensions. 




        

<hr>



### function NDLagrangeEvaluator [2/3]

_Copy-construct._ 
```C++
NDLagrangeEvaluator::NDLagrangeEvaluator (
    NDLagrangeEvaluator const & x
) = default
```




<hr>



### function NDLagrangeEvaluator [3/3]

_Move-construct._ 
```C++
NDLagrangeEvaluator::NDLagrangeEvaluator (
    NDLagrangeEvaluator && x
) = default
```




<hr>



### function deriv [1/3]

_Differentiate 1D Lagrange function at a given coordinate._ 
```C++
template<class IdxDerivDims, class Layout, class BatchedLagrangeIdxRange, class... CoordsDims>
inline KOKKOS_FUNCTION double NDLagrangeEvaluator::deriv (
    IdxDerivDims const & deriv_order,
    Coord< CoordsDims... > const & coord_eval,
    ConstField< data_type , BatchedLagrangeIdxRange, memory_space , Layout > const lagrange_coef
) const
```





**Parameters:**


* `deriv_order` An Idx containing the order of derivation for the dimension of interest. If the dimension is not present, the order of derivation is considered to be 0. 
* `coord_eval` The coordinate where the polynomial is differentiated. Note that only the component along the dimension of interest is used. 
* `lagrange_coef` A Field storing the 1D Lagrange coefficients.



**Returns:**

The derivative of the spline function at the desired coordinate. 





        

<hr>



### function deriv [2/3]

_Differentiate 1D spline function (described by its spline coefficients) on a mesh._ 
```C++
template<class IdxDerivDims, class Layout1, class Layout2, class Layout3, class IdxRangeBatchedInterpolation, class... CoordsDims>
inline void NDLagrangeEvaluator::deriv (
    IdxDerivDims const & deriv_order,
    Field< data_type , IdxRangeBatchedInterpolation, memory_space , Layout1 > const lagrange_eval,
    ConstField< Coord< CoordsDims... >, IdxRangeBatchedInterpolation, memory_space , Layout2 > const coords_eval,
    ConstField< data_type , batched_coeff_idx_range_type < IdxRangeBatchedInterpolation >, memory_space , Layout3 > const lagrange_coef
) const
```



The spline coefficients represent a spline function defined on a cartesian product of batch\_domain and B-splines (basis splines). They can be obtained via various methods, such as using a SplineBuilder.


The derivation is not performed in a multidimensional way (in any sense). This is a batched 1D derivation. This means that for each slice of coordinates identified by a batch\_domain\_type::discrete\_element\_type, the derivation is performed with the 1D set of spline coefficients identified by the same batch\_domain\_type::discrete\_element\_type.




**Parameters:**


* `deriv_order` An Idx containing the order of derivation for the dimension of interest. If the dimension is not present, the order of derivation is considered to be 0. 
* `lagrange_eval` The values of the derivatives of the Lagrange polynomials at the desired coordinates. 
* `coords_eval` The coordinates where the Lagrange polynomials are evaluated. Those are stored in a Field defined on a BatchedInterpolationIdxRange. Note that the coordinates of the points represented by this index range are unused and irrelevant. 
* `lagrange_coef` A Field storing the 1D Lagrange coefficients. 




        

<hr>



### function deriv [3/3]

_Differentiate 1D Lagrange function on a mesh._ 
```C++
template<class IdxDerivDims, class Layout1, class Layout2, class IdxRangeBatchedInterpolation>
inline void NDLagrangeEvaluator::deriv (
    IdxDerivDims const & deriv_order,
    Field< data_type , IdxRangeBatchedInterpolation, memory_space , Layout1 > const lagrange_eval,
    ConstField< data_type , batched_coeff_idx_range_type < IdxRangeBatchedInterpolation >, memory_space , Layout2 > const lagrange_coef
) const
```



The differentiation is not performed in a multidimensional way (in any sense). This is a batched 1D derivation. This means that for each slice of lagrange\_eval the evaluation is performed with the relevant 1D set of Lagrange coefficients.




**Parameters:**


* `deriv_order` An Idx containing the order of derivation for the dimension of interest. If the dimension is not present, the order of derivation is considered to be 0. 
* `lagrange_eval` The derivatives of the spline function at the coordinates. 
* `lagrange_coef` A ChunkSpan storing the spline coefficients. 




        

<hr>



### function operator() 

_Evaluate the ND Lagrange polynomial at a single coordinate._ 
```C++
template<class Layout, class... CoordsDims, class NdCoeffIdxRange>
inline KOKKOS_FUNCTION data_type NDLagrangeEvaluator::operator() (
    Coord< CoordsDims... > const & coord,
    ConstField< data_type , NdCoeffIdxRange, memory_space , Layout > const & coeff
) const
```





**Parameters:**


* `coord` The ND evaluation coordinate. 
* `coeff` The ND Lagrange coefficient field (no batch dimensions). 



**Returns:**

The interpolated value. 





        

<hr>



### function operator() 

_Evaluate the ND Lagrange polynomial on a mesh (with explicit coordinates)._ 
```C++
template<class Layout1, class Layout2, class Layout3, class IdxRangeBatched, class... CoordsDims>
inline void NDLagrangeEvaluator::operator() (
    Field< data_type , IdxRangeBatched, memory_space , Layout1 > lagrange_eval,
    ConstField< Coord< CoordsDims... >, IdxRangeBatched, memory_space , Layout2 > coords_eval,
    ConstField< data_type , batched_coeff_idx_range_type < IdxRangeBatched >, memory_space , Layout3 > lagrange_coef
) const
```



The evaluation is parallelised over the full (batch + ND evaluation) domain. For each point the single-point operator() is called with the corresponding coordinate and batch-sliced coefficient field.




**Parameters:**


* `lagrange_eval` The interpolated values on the full batched domain. 
* `coords_eval` The evaluation coordinates on the full batched domain. 
* `lagrange_coef` The Lagrange coefficients on the batched coefficient domain. 




        

<hr>



### function operator() 

_Evaluate the ND Lagrange polynomial on a mesh (coordinates from the grid)._ 
```C++
template<class Layout1, class Layout2, class IdxRangeBatched>
inline void NDLagrangeEvaluator::operator() (
    Field< data_type , IdxRangeBatched, memory_space , Layout1 > lagrange_eval,
    ConstField< data_type , batched_coeff_idx_range_type < IdxRangeBatched >, memory_space , Layout2 > lagrange_coef
) const
```



Coordinates are derived from the evaluation grid indices via ddc::coordinate. This variant requires all evaluation grids to be genuine discrete approximations of continuous dimensions.




**Parameters:**


* `lagrange_eval` The interpolated values on the full batched domain. 
* `lagrange_coef` The Lagrange coefficients on the batched coefficient domain. 




        

<hr>



### function operator= 

_Copy-assign._ 
```C++
NDLagrangeEvaluator & NDLagrangeEvaluator::operator= (
    NDLagrangeEvaluator const & x
) = default
```




<hr>



### function operator= 

_Move-assign._ 
```C++
NDLagrangeEvaluator & NDLagrangeEvaluator::operator= (
    NDLagrangeEvaluator && x
) = default
```




<hr>



### function ~NDLagrangeEvaluator 

_Destruct._ 
```C++
NDLagrangeEvaluator::~NDLagrangeEvaluator () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/nd_lagrange_evaluator.hpp`

