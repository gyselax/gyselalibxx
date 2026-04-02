

# Class LagrangeEvaluator

**template &lt;class ExecSpace, class MemorySpace, class DataType, class LagrangeBasis, class InterpolationGrid, class LowerExtrapolationRule, class UpperExtrapolationRule&gt;**



[**ClassList**](annotated.md) **>** [**LagrangeEvaluator**](classLagrangeEvaluator.md)



_A class to evaluate, differentiate or integrate a Lagrange function._ [More...](#detailed-description)

* `#include <lagrange_evaluator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::remove\_dims\_of\_t&lt; BatchedInterpolationIdxRange, InterpolationGrid &gt; | [**batch\_idx\_range\_type**](#typedef-batch_idx_range_type)  <br>_The type of the batch domain (obtained by removing the dimension of interest from the whole domain)._  |
| typedef ddc::replace\_dim\_of\_t&lt; BatchedInterpolationIdxRange, InterpolationGrid, [**coeff\_grid\_type**](classLagrangeEvaluator.md#typedef-coeff_grid_type) &gt; | [**batched\_coeff\_idx\_range\_type**](#typedef-batched_coeff_idx_range_type)  <br>_The type of the whole Lagrange domain (cartesian product of 1D Lagrange domain and batch domain) preserving the order of dimensions._  |
| typedef BatchedInterpolationIdxRange | [**batched\_evaluation\_idx\_range\_type**](#typedef-batched_evaluation_idx_range_type)  <br>_The type of the whole domain representing evaluation points._  |
| typedef typename LagrangeBasis::template Impl&lt; LagrangeBasis, MemorySpace &gt;::knot\_grid | [**coeff\_grid\_type**](#typedef-coeff_grid_type)  <br>_The grid on which the interpolation coefficients should be provided._  |
| typedef IdxRange&lt; [**coeff\_grid\_type**](classLagrangeEvaluator.md#typedef-coeff_grid_type) &gt; | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The type of the 1D Lagrange domain corresponding to the dimension of interest._  |
| typedef typename LagrangeBasis::continuous\_dimension\_type | [**continuous\_dimension\_type**](#typedef-continuous_dimension_type)  <br>_The type of the evaluation continuous dimension (continuous dimension of interest) used by this class._  |
| typedef DataType | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef IdxRange&lt; InterpolationGrid &gt; | [**evaluation\_idx\_range\_type**](#typedef-evaluation_idx_range_type)  <br>_The type of the domain for the 1D evaluation mesh used by this class._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The type of the Kokkos execution space used by this class._  |
| typedef LagrangeBasis | [**lagrange\_basis\_type**](#typedef-lagrange_basis_type)  <br>_The discrete dimension representing the Lagrange basis._  |
| typedef LowerExtrapolationRule | [**lower\_extrapolation\_rule\_type**](#typedef-lower_extrapolation_rule_type)  <br>_The type of the extrapolation rule at the lower boundary._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The type of the Kokkos memory space used by this class._  |
| typedef UpperExtrapolationRule | [**upper\_extrapolation\_rule\_type**](#typedef-upper_extrapolation_rule_type)  <br>_The type of the extrapolation rule at the upper boundary._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**LagrangeEvaluator**](#function-lagrangeevaluator-13) (LowerExtrapolationRule const & lower\_extrap\_rule, UpperExtrapolationRule const & upper\_extrap\_rule) <br>_Build a_ [_**LagrangeEvaluator**_](classLagrangeEvaluator.md) _acting on batched\_lagrange\_domain._ |
|   | [**LagrangeEvaluator**](#function-lagrangeevaluator-23) ([**LagrangeEvaluator**](classLagrangeEvaluator.md) const & x) = default<br>_Copy-constructs._  |
|   | [**LagrangeEvaluator**](#function-lagrangeevaluator-33) ([**LagrangeEvaluator**](classLagrangeEvaluator.md) && x) = default<br>_Move-constructs._  |
|  KOKKOS\_FUNCTION DataType | [**deriv**](#function-deriv-13) (IdxDerivDims const & deriv\_order, Coord&lt; CoordsDims... &gt; const & coord\_eval, ConstField&lt; DataType, BatchedLagrangeIdxRange, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout &gt; const lagrange\_coef) const<br>_Differentiate 1D Lagrange function at a given coordinate._  |
|  void | [**deriv**](#function-deriv-23) (IdxDerivDims const & deriv\_order, Field&lt; DataType, IdxRangeBatchedInterpolation, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; const lagrange\_eval, ConstField&lt; Coord&lt; CoordsDims... &gt;, IdxRangeBatchedInterpolation, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; const coords\_eval, ConstField&lt; DataType, [**batched\_coeff\_idx\_range\_type**](classLagrangeEvaluator.md#typedef-batched_coeff_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt;, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout3 &gt; const lagrange\_coef) const<br>_Differentiate 1D spline function (described by its spline coefficients) on a mesh._  |
|  void | [**deriv**](#function-deriv-33) (IdxDerivDims const & deriv\_order, Field&lt; DataType, IdxRangeBatchedInterpolation, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; const lagrange\_eval, ConstField&lt; DataType, [**batched\_coeff\_idx\_range\_type**](classLagrangeEvaluator.md#typedef-batched_coeff_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt;, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; const lagrange\_coef) const<br>_Differentiate 1D Lagrange function on a mesh._  |
|  [**lower\_extrapolation\_rule\_type**](classLagrangeEvaluator.md#typedef-lower_extrapolation_rule_type) | [**lower\_extrapolation\_rule**](#function-lower_extrapolation_rule) () const<br>_Get the lower extrapolation rule._  |
|  KOKKOS\_FUNCTION DataType | [**operator()**](#function-operator) (Coord&lt; CoordsDims... &gt; const & coord\_eval, ConstField&lt; DataType, BatchedLagrangeIdxRange, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout &gt; const lagrange\_coef) const<br>_Evaluate 1D Lagrange polynomial (described by its Lagrange coefficients) at a given coordinate._  |
|  void | [**operator()**](#function-operator_1) (Field&lt; DataType, IdxRangeBatchedInterpolation, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; const lagrange\_eval, ConstField&lt; Coord&lt; CoordsDims... &gt;, IdxRangeBatchedInterpolation, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; const coords\_eval, ConstField&lt; DataType, [**batched\_coeff\_idx\_range\_type**](classLagrangeEvaluator.md#typedef-batched_coeff_idx_range_type)&lt; IdxRangeBatchedInterpolation &gt;, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout3 &gt; const lagrange\_coef) const<br>_Evaluate Lagrange polynomials (described by their Lagrange coefficients) on a mesh._  |
|  void | [**operator()**](#function-operator_2) (Field&lt; DataType, BatchedInterpolationIdxRange, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; const lagrange\_eval, ConstField&lt; DataType, [**batched\_coeff\_idx\_range\_type**](classLagrangeEvaluator.md#typedef-batched_coeff_idx_range_type)&lt; BatchedInterpolationIdxRange &gt;, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; const lagrange\_coef) const<br>_Evaluate Lagrange polynomials (described by their Lagrange coefficients) on a mesh._  |
|  [**LagrangeEvaluator**](classLagrangeEvaluator.md) & | [**operator=**](#function-operator_3) ([**LagrangeEvaluator**](classLagrangeEvaluator.md) const & x) = default<br>_Copy-assigns._  |
|  [**LagrangeEvaluator**](classLagrangeEvaluator.md) & | [**operator=**](#function-operator_4) ([**LagrangeEvaluator**](classLagrangeEvaluator.md) && x) = default<br>_Move-assigns._  |
|  [**upper\_extrapolation\_rule\_type**](classLagrangeEvaluator.md#typedef-upper_extrapolation_rule_type) | [**upper\_extrapolation\_rule**](#function-upper_extrapolation_rule) () const<br>_Get the upper extrapolation rule._  |
|   | [**~LagrangeEvaluator**](#function-lagrangeevaluator) () = default<br>_Destructs._  |




























## Detailed Description


A class which contains an operator () which can be used to evaluate, or differentiate a Lagrange polynomial.




**Template parameters:**


* `ExecSpace` The Kokkos execution space on which the Lagrange evaluation is performed. 
* `MemorySpace` The Kokkos memory space on which the data (Lagrange coefficients and evaluation) is stored. 
* `DataType` The data type on which calculations are made. 
* `LagrangeBasis` The discrete dimension representing the Lagrange basis. 
* `InterpolationGrid` The discrete dimension on which evaluation points are defined. 
* `LowerExtrapolationRule` The lower extrapolation rule type. 
* `UpperExtrapolationRule` The upper extrapolation rule type. 




    
## Public Types Documentation




### typedef batch\_idx\_range\_type 

_The type of the batch domain (obtained by removing the dimension of interest from the whole domain)._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::batch_idx_range_type =  ddc::remove_dims_of_t<BatchedInterpolationIdxRange, InterpolationGrid>;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef batched\_coeff\_idx\_range\_type 

_The type of the whole Lagrange domain (cartesian product of 1D Lagrange domain and batch domain) preserving the order of dimensions._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::batched_coeff_idx_range_type =  ddc:: replace_dim_of_t<BatchedInterpolationIdxRange, InterpolationGrid, coeff_grid_type>;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef batched\_evaluation\_idx\_range\_type 

_The type of the whole domain representing evaluation points._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::batched_evaluation_idx_range_type =  BatchedInterpolationIdxRange;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef coeff\_grid\_type 

_The grid on which the interpolation coefficients should be provided._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::coeff_grid_type =  typename LagrangeBasis::template Impl<LagrangeBasis, MemorySpace>::knot_grid;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The type of the 1D Lagrange domain corresponding to the dimension of interest._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::coeff_idx_range_type =  IdxRange<coeff_grid_type>;
```




<hr>



### typedef continuous\_dimension\_type 

_The type of the evaluation continuous dimension (continuous dimension of interest) used by this class._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::continuous_dimension_type =  typename LagrangeBasis::continuous_dimension_type;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::data_type =  DataType;
```




<hr>



### typedef evaluation\_idx\_range\_type 

_The type of the domain for the 1D evaluation mesh used by this class._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::evaluation_idx_range_type =  IdxRange<InterpolationGrid>;
```




<hr>



### typedef exec\_space 

_The type of the Kokkos execution space used by this class._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::exec_space =  ExecSpace;
```




<hr>



### typedef lagrange\_basis\_type 

_The discrete dimension representing the Lagrange basis._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::lagrange_basis_type =  LagrangeBasis;
```




<hr>



### typedef lower\_extrapolation\_rule\_type 

_The type of the extrapolation rule at the lower boundary._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::lower_extrapolation_rule_type =  LowerExtrapolationRule;
```




<hr>



### typedef memory\_space 

_The type of the Kokkos memory space used by this class._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::memory_space =  MemorySpace;
```




<hr>



### typedef upper\_extrapolation\_rule\_type 

_The type of the extrapolation rule at the upper boundary._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::upper_extrapolation_rule_type =  UpperExtrapolationRule;
```




<hr>
## Public Functions Documentation




### function LagrangeEvaluator [1/3]

_Build a_ [_**LagrangeEvaluator**_](classLagrangeEvaluator.md) _acting on batched\_lagrange\_domain._
```C++
inline explicit LagrangeEvaluator::LagrangeEvaluator (
    LowerExtrapolationRule const & lower_extrap_rule,
    UpperExtrapolationRule const & upper_extrap_rule
) 
```





**Parameters:**


* `lower_extrap_rule` The extrapolation rule at the lower boundary. 
* `upper_extrap_rule` The extrapolation rule at the upper boundary.



**See also:** [**NullExtrapolationRule**](structNullExtrapolationRule.md) ConstantExtrapolationRule PeriodicExtrapolationRule 



        

<hr>



### function LagrangeEvaluator [2/3]

_Copy-constructs._ 
```C++
LagrangeEvaluator::LagrangeEvaluator (
    LagrangeEvaluator const & x
) = default
```





**Parameters:**


* `x` A reference to another [**LagrangeEvaluator**](classLagrangeEvaluator.md). 




        

<hr>



### function LagrangeEvaluator [3/3]

_Move-constructs._ 
```C++
LagrangeEvaluator::LagrangeEvaluator (
    LagrangeEvaluator && x
) = default
```





**Parameters:**


* `x` An rvalue to another [**LagrangeEvaluator**](classLagrangeEvaluator.md). 




        

<hr>



### function deriv [1/3]

_Differentiate 1D Lagrange function at a given coordinate._ 
```C++
template<class IdxDerivDims, class Layout, class BatchedLagrangeIdxRange, class... CoordsDims>
inline KOKKOS_FUNCTION DataType LagrangeEvaluator::deriv (
    IdxDerivDims const & deriv_order,
    Coord< CoordsDims... > const & coord_eval,
    ConstField< DataType, BatchedLagrangeIdxRange, memory_space , Layout > const lagrange_coef
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
inline void LagrangeEvaluator::deriv (
    IdxDerivDims const & deriv_order,
    Field< DataType, IdxRangeBatchedInterpolation, memory_space , Layout1 > const lagrange_eval,
    ConstField< Coord< CoordsDims... >, IdxRangeBatchedInterpolation, memory_space , Layout2 > const coords_eval,
    ConstField< DataType, batched_coeff_idx_range_type < IdxRangeBatchedInterpolation >, memory_space , Layout3 > const lagrange_coef
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
inline void LagrangeEvaluator::deriv (
    IdxDerivDims const & deriv_order,
    Field< DataType, IdxRangeBatchedInterpolation, memory_space , Layout1 > const lagrange_eval,
    ConstField< DataType, batched_coeff_idx_range_type < IdxRangeBatchedInterpolation >, memory_space , Layout2 > const lagrange_coef
) const
```



The differentiation is not performed in a multidimensional way (in any sense). This is a batched 1D derivation. This means that for each slice of lagrange\_eval the evaluation is performed with the relevant 1D set of Lagrange coefficients.




**Parameters:**


* `deriv_order` An Idx containing the order of derivation for the dimension of interest. If the dimension is not present, the order of derivation is considered to be 0. 
* `lagrange_eval` The derivatives of the spline function at the coordinates. 
* `lagrange_coef` A ChunkSpan storing the spline coefficients. 




        

<hr>



### function lower\_extrapolation\_rule 

_Get the lower extrapolation rule._ 
```C++
inline lower_extrapolation_rule_type LagrangeEvaluator::lower_extrapolation_rule () const
```



Extrapolation rules are functors used to define the behaviour of the [**LagrangeEvaluator**](classLagrangeEvaluator.md) outside the domain where the break points of the Lagrange basis are defined.




**Returns:**

The lower extrapolation rule.




**See also:** [**NullExtrapolationRule**](structNullExtrapolationRule.md) ConstantExtrapolationRule PeriodicExtrapolationRule 



        

<hr>



### function operator() 

_Evaluate 1D Lagrange polynomial (described by its Lagrange coefficients) at a given coordinate._ 
```C++
template<class Layout, class... CoordsDims, class BatchedLagrangeIdxRange>
inline KOKKOS_FUNCTION DataType LagrangeEvaluator::operator() (
    Coord< CoordsDims... > const & coord_eval,
    ConstField< DataType, BatchedLagrangeIdxRange, memory_space , Layout > const lagrange_coef
) const
```



The Lagrange coefficients describe 1D Lagrange polynomials defined on a Lagrange basis. They can be obtained using the [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md).




**Parameters:**


* `coord_eval` The coordinate where the Lagrange polynomial is evaluated. Note that only the component along the dimension of interest is used. 
* `lagrange_coef` A Field storing the 1D Lagrange coefficients.



**Returns:**

The value of the Lagrange polynomial at the desired coordinate. 





        

<hr>



### function operator() 

_Evaluate Lagrange polynomials (described by their Lagrange coefficients) on a mesh._ 
```C++
template<class Layout1, class Layout2, class Layout3, class IdxRangeBatchedInterpolation, class... CoordsDims>
inline void LagrangeEvaluator::operator() (
    Field< DataType, IdxRangeBatchedInterpolation, memory_space , Layout1 > const lagrange_eval,
    ConstField< Coord< CoordsDims... >, IdxRangeBatchedInterpolation, memory_space , Layout2 > const coords_eval,
    ConstField< DataType, batched_coeff_idx_range_type < IdxRangeBatchedInterpolation >, memory_space , Layout3 > const lagrange_coef
) const
```



The Lagrange coefficients describe Lagrange polynomials defined on a cartesian product of batch\_idx\_range and Lagrange basis. They can be obtained using the [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md).


This is not a multidimensional evaluation. This is a batched 1D evaluation. This means that for each slice of coordinates identified by a batch\_idx\_range\_type::discrete\_element\_type, the evaluation is performed with the 1D set of Lagrange coefficients identified by the same batch\_idx\_range\_type::discrete\_element\_type.




**Parameters:**


* `lagrange_eval` The values of the Lagrange polynomials at the desired coordinates. 
* `coords_eval` The coordinates where the Lagrange polynomials are evaluated. Those are stored in a Field defined on a BatchedInterpolationIdxRange. Note that the coordinates of the points represented by this index range are unused and irrelevant. 
* `lagrange_coef` A Field storing the Lagrange coefficients. 




        

<hr>



### function operator() 

_Evaluate Lagrange polynomials (described by their Lagrange coefficients) on a mesh._ 
```C++
template<class Layout1, class Layout2, class BatchedInterpolationIdxRange>
inline void LagrangeEvaluator::operator() (
    Field< DataType, BatchedInterpolationIdxRange, memory_space , Layout1 > const lagrange_eval,
    ConstField< DataType, batched_coeff_idx_range_type < BatchedInterpolationIdxRange >, memory_space , Layout2 > const lagrange_coef
) const
```



The Lagrange coefficients describe Lagrange polynomials defined on a cartesian product of batch\_idx\_range and Lagrange basis. They can be obtained using the [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md).


This is not a multidimensional evaluation. This is a batched 1D evaluation. This means that for each slice of coordinates identified by a batch\_idx\_range\_type::discrete\_element\_type, the evaluation is performed with the 1D set of Lagrange coefficients identified by the same batch\_idx\_range\_type::discrete\_element\_type.




**Parameters:**


* `lagrange_eval` The values of the Lagrange polynomials at the mesh points. 
* `lagrange_coef` A Field storing the Lagrange coefficients. 




        

<hr>



### function operator= 

_Copy-assigns._ 
```C++
LagrangeEvaluator & LagrangeEvaluator::operator= (
    LagrangeEvaluator const & x
) = default
```





**Parameters:**


* `x` A reference to another [**LagrangeEvaluator**](classLagrangeEvaluator.md). 



**Returns:**

A reference to this object. 





        

<hr>



### function operator= 

_Move-assigns._ 
```C++
LagrangeEvaluator & LagrangeEvaluator::operator= (
    LagrangeEvaluator && x
) = default
```





**Parameters:**


* `x` An rvalue to another [**LagrangeEvaluator**](classLagrangeEvaluator.md). 



**Returns:**

A reference to this object. 





        

<hr>



### function upper\_extrapolation\_rule 

_Get the upper extrapolation rule._ 
```C++
inline upper_extrapolation_rule_type LagrangeEvaluator::upper_extrapolation_rule () const
```



Extrapolation rules are functors used to define the behaviour of the [**LagrangeEvaluator**](classLagrangeEvaluator.md) outside the domain where the break points of the Lagrange basis are defined.




**Returns:**

The upper extrapolation rule.




**See also:** [**NullExtrapolationRule**](structNullExtrapolationRule.md) ConstantExtrapolationRule PeriodicExtrapolationRule 



        

<hr>



### function ~LagrangeEvaluator 

_Destructs._ 
```C++
LagrangeEvaluator::~LagrangeEvaluator () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/lagrange_evaluator.hpp`

