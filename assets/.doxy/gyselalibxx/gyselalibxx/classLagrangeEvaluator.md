

# Class LagrangeEvaluator

**template &lt;class ExecSpace, class MemorySpace, class DataType, class LagrangeBasis, class InterpolationGrid, class LowerExtrapolationRule, class UpperExtrapolationRule&gt;**



[**ClassList**](annotated.md) **>** [**LagrangeEvaluator**](classLagrangeEvaluator.md)



_A class to evaluate, differentiate or integrate a_ [_**Lagrange**_](classLagrange.md) _function._[More...](#detailed-description)

* `#include <lagrange_evaluator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename LagrangeBasis::template Impl&lt; LagrangeBasis, MemorySpace &gt;::knot\_grid | [**basis\_domain\_type**](#typedef-basis_domain_type)  <br>_The grid on which the interpolation coefficients should be provided._  |
| typedef ddc::remove\_dims\_of\_t&lt; BatchedInterpolationGrid, InterpolationGrid &gt; | [**batch\_domain\_type**](#typedef-batch_domain_type)  <br>_The type of the batch domain (obtained by removing the dimension of interest from the whole domain)._  |
| typedef BatchedInterpolationGrid | [**batched\_evaluation\_domain\_type**](#typedef-batched_evaluation_domain_type)  <br>_The type of the whole domain representing evaluation points._  |
| typedef ddc::replace\_dim\_of\_t&lt; BatchedInterpolationDDom, InterpolationGrid, [**basis\_domain\_type**](classLagrangeEvaluator.md#typedef-basis_domain_type) &gt; | [**batched\_lagrange\_domain\_type**](#typedef-batched_lagrange_domain_type)  <br>_The type of the whole_ [_**Lagrange**_](classLagrange.md) _domain (cartesian product of 1D_[_**Lagrange**_](classLagrange.md) _domain and batch domain) preserving the order of dimensions._ |
| typedef typename LagrangeBasis::continuous\_dimension\_type | [**continuous\_dimension\_type**](#typedef-continuous_dimension_type)  <br>_The type of the evaluation continuous dimension (continuous dimension of interest) used by this class._  |
| typedef IdxRange&lt; InterpolationGrid &gt; | [**evaluation\_domain\_type**](#typedef-evaluation_domain_type)  <br>_The type of the domain for the 1D evaluation mesh used by this class._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The type of the Kokkos execution space used by this class._  |
| typedef LagrangeBasis | [**lagrange\_basis\_type**](#typedef-lagrange_basis_type)  <br>_The discrete dimension representing the_ [_**Lagrange**_](classLagrange.md) _basis._ |
| typedef IdxRange&lt; [**basis\_domain\_type**](classLagrangeEvaluator.md#typedef-basis_domain_type) &gt; | [**lagrange\_domain\_type**](#typedef-lagrange_domain_type)  <br>_The type of the 1D_ [_**Lagrange**_](classLagrange.md) _domain corresponding to the dimension of interest._ |
| typedef LowerExtrapolationRule | [**lower\_extrapolation\_rule\_type**](#typedef-lower_extrapolation_rule_type)  <br>_The type of the extrapolation rule at the lower boundary._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The type of the Kokkos memory space used by this class._  |
| typedef UpperExtrapolationRule | [**upper\_extrapolation\_rule\_type**](#typedef-upper_extrapolation_rule_type)  <br>_The type of the extrapolation rule at the upper boundary._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**LagrangeEvaluator**](#function-lagrangeevaluator-13) (LowerExtrapolationRule const & lower\_extrap\_rule, UpperExtrapolationRule const & upper\_extrap\_rule) <br>_Build a_ [_**LagrangeEvaluator**_](classLagrangeEvaluator.md) _acting on batched\_lagrange\_domain._ |
|   | [**LagrangeEvaluator**](#function-lagrangeevaluator-23) ([**LagrangeEvaluator**](classLagrangeEvaluator.md) const & x) = default<br>_Copy-constructs._  |
|   | [**LagrangeEvaluator**](#function-lagrangeevaluator-33) ([**LagrangeEvaluator**](classLagrangeEvaluator.md) && x) = default<br>_Move-constructs._  |
|  [**lower\_extrapolation\_rule\_type**](classLagrangeEvaluator.md#typedef-lower_extrapolation_rule_type) | [**lower\_extrapolation\_rule**](#function-lower_extrapolation_rule) () const<br>_Get the lower extrapolation rule._  |
|  KOKKOS\_FUNCTION DataType | [**operator()**](#function-operator) (Coord&lt; CoordsDims... &gt; const & coord\_eval, ConstField&lt; DataType, BatchedLagrangeIdxRange, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout &gt; const lagrange\_coef) const<br>_Evaluate 1D_ [_**Lagrange**_](classLagrange.md) _polynomial (described by its_[_**Lagrange**_](classLagrange.md) _coefficients) at a given coordinate._ |
|  void | [**operator()**](#function-operator_1) (Field&lt; DataType, IdxRangeBatchedInterpolation, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; const lagrange\_eval, ConstField&lt; Coord&lt; CoordsDims... &gt;, IdxRangeBatchedInterpolation, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; const coords\_eval, ConstField&lt; DataType, [**batched\_lagrange\_domain\_type**](classLagrangeEvaluator.md#typedef-batched_lagrange_domain_type)&lt; IdxRangeBatchedInterpolation &gt;, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout3 &gt; const lagrange\_coef) const<br>_Evaluate_ [_**Lagrange**_](classLagrange.md) _polynomials (described by their_[_**Lagrange**_](classLagrange.md) _coefficients) on a mesh._ |
|  void | [**operator()**](#function-operator_2) (Field&lt; DataType, BatchedInterpolationGrid, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout1 &gt; const lagrange\_eval, ConstField&lt; DataType, [**batched\_lagrange\_domain\_type**](classLagrangeEvaluator.md#typedef-batched_lagrange_domain_type)&lt; BatchedInterpolationGrid &gt;, [**memory\_space**](classLagrangeEvaluator.md#typedef-memory_space), Layout2 &gt; const lagrange\_coef) const<br>_Evaluate_ [_**Lagrange**_](classLagrange.md) _polynomials (described by their_[_**Lagrange**_](classLagrange.md) _coefficients) on a mesh._ |
|  [**LagrangeEvaluator**](classLagrangeEvaluator.md) & | [**operator=**](#function-operator_3) ([**LagrangeEvaluator**](classLagrangeEvaluator.md) const & x) = default<br>_Copy-assigns._  |
|  [**LagrangeEvaluator**](classLagrangeEvaluator.md) & | [**operator=**](#function-operator_4) ([**LagrangeEvaluator**](classLagrangeEvaluator.md) && x) = default<br>_Move-assigns._  |
|  [**upper\_extrapolation\_rule\_type**](classLagrangeEvaluator.md#typedef-upper_extrapolation_rule_type) | [**upper\_extrapolation\_rule**](#function-upper_extrapolation_rule) () const<br>_Get the upper extrapolation rule._  |
|   | [**~LagrangeEvaluator**](#function-lagrangeevaluator) () = default<br>_Destructs._  |




























## Detailed Description


A class which contains an operator () which can be used to evaluate, or differentiate a [**Lagrange**](classLagrange.md) polynomial.




**Template parameters:**


* `ExecSpace` The Kokkos execution space on which the [**Lagrange**](classLagrange.md) evaluation is performed. 
* `MemorySpace` The Kokkos memory space on which the data ([**Lagrange**](classLagrange.md) coefficients and evaluation) is stored. 
* `DataType` The data type on which calculations are made. 
* `LagrangeBasis` The discrete dimension representing the [**Lagrange**](classLagrange.md) basis. 
* `InterpolationGrid` The discrete dimension on which evaluation points are defined. 
* `LowerExtrapolationRule` The lower extrapolation rule type. 
* `UpperExtrapolationRule` The upper extrapolation rule type. 




    
## Public Types Documentation




### typedef basis\_domain\_type 

_The grid on which the interpolation coefficients should be provided._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::basis_domain_type =  typename LagrangeBasis::template Impl<LagrangeBasis, MemorySpace>::knot_grid;
```




<hr>



### typedef batch\_domain\_type 

_The type of the batch domain (obtained by removing the dimension of interest from the whole domain)._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::batch_domain_type =  ddc::remove_dims_of_t<BatchedInterpolationGrid, InterpolationGrid>;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef batched\_evaluation\_domain\_type 

_The type of the whole domain representing evaluation points._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::batched_evaluation_domain_type =  BatchedInterpolationGrid;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef batched\_lagrange\_domain\_type 

_The type of the whole_ [_**Lagrange**_](classLagrange.md) _domain (cartesian product of 1D_[_**Lagrange**_](classLagrange.md) _domain and batch domain) preserving the order of dimensions._
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::batched_lagrange_domain_type =  ddc::replace_dim_of_t<BatchedInterpolationDDom, InterpolationGrid, basis_domain_type>;
```





**Template parameters:**


* `The` batched discrete domain on which the interpolation points are defined. 




        

<hr>



### typedef continuous\_dimension\_type 

_The type of the evaluation continuous dimension (continuous dimension of interest) used by this class._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::continuous_dimension_type =  typename LagrangeBasis::continuous_dimension_type;
```




<hr>



### typedef evaluation\_domain\_type 

_The type of the domain for the 1D evaluation mesh used by this class._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::evaluation_domain_type =  IdxRange<InterpolationGrid>;
```




<hr>



### typedef exec\_space 

_The type of the Kokkos execution space used by this class._ 
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::exec_space =  ExecSpace;
```




<hr>



### typedef lagrange\_basis\_type 

_The discrete dimension representing the_ [_**Lagrange**_](classLagrange.md) _basis._
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::lagrange_basis_type =  LagrangeBasis;
```




<hr>



### typedef lagrange\_domain\_type 

_The type of the 1D_ [_**Lagrange**_](classLagrange.md) _domain corresponding to the dimension of interest._
```C++
using LagrangeEvaluator< ExecSpace, MemorySpace, DataType, LagrangeBasis, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule >::lagrange_domain_type =  IdxRange<basis_domain_type>;
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



### function lower\_extrapolation\_rule 

_Get the lower extrapolation rule._ 
```C++
inline lower_extrapolation_rule_type LagrangeEvaluator::lower_extrapolation_rule () const
```



Extrapolation rules are functors used to define the behaviour of the [**LagrangeEvaluator**](classLagrangeEvaluator.md) outside the domain where the break points of the [**Lagrange**](classLagrange.md) basis are defined.




**Returns:**

The lower extrapolation rule.




**See also:** [**NullExtrapolationRule**](structNullExtrapolationRule.md) ConstantExtrapolationRule PeriodicExtrapolationRule 



        

<hr>



### function operator() 

_Evaluate 1D_ [_**Lagrange**_](classLagrange.md) _polynomial (described by its_[_**Lagrange**_](classLagrange.md) _coefficients) at a given coordinate._
```C++
template<class Layout, class... CoordsDims, class BatchedLagrangeIdxRange>
inline KOKKOS_FUNCTION DataType LagrangeEvaluator::operator() (
    Coord< CoordsDims... > const & coord_eval,
    ConstField< DataType, BatchedLagrangeIdxRange, memory_space , Layout > const lagrange_coef
) const
```



The [**Lagrange**](classLagrange.md) coefficients describe 1D [**Lagrange**](classLagrange.md) polynomials defined on a [**Lagrange**](classLagrange.md) basis. They can be obtained using the [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md).




**Parameters:**


* `coord_eval` The coordinate where the [**Lagrange**](classLagrange.md) polynomial is evaluated. Note that only the component along the dimension of interest is used. 
* `lagrange_coef` A Field storing the 1D [**Lagrange**](classLagrange.md) coefficients.



**Returns:**

The value of the [**Lagrange**](classLagrange.md) polynomial at the desired coordinate. 





        

<hr>



### function operator() 

_Evaluate_ [_**Lagrange**_](classLagrange.md) _polynomials (described by their_[_**Lagrange**_](classLagrange.md) _coefficients) on a mesh._
```C++
template<class Layout1, class Layout2, class Layout3, class IdxRangeBatchedInterpolation, class... CoordsDims>
inline void LagrangeEvaluator::operator() (
    Field< DataType, IdxRangeBatchedInterpolation, memory_space , Layout1 > const lagrange_eval,
    ConstField< Coord< CoordsDims... >, IdxRangeBatchedInterpolation, memory_space , Layout2 > const coords_eval,
    ConstField< DataType, batched_lagrange_domain_type < IdxRangeBatchedInterpolation >, memory_space , Layout3 > const lagrange_coef
) const
```



The [**Lagrange**](classLagrange.md) coefficients describe [**Lagrange**](classLagrange.md) polynomials defined on a cartesian product of batch\_idx\_range and [**Lagrange**](classLagrange.md) basis. They can be obtained using the [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md).


This is not a multidimensional evaluation. This is a batched 1D evaluation. This means that for each slice of coordinates identified by a batch\_domain\_type::discrete\_element\_type, the evaluation is performed with the 1D set of [**Lagrange**](classLagrange.md) coefficients identified by the same batch\_domain\_type::discrete\_element\_type.




**Parameters:**


* `lagrange_eval` The values of the [**Lagrange**](classLagrange.md) polynomials at the desired coordinates. 
* `coords_eval` The coordinates where the [**Lagrange**](classLagrange.md) polynomials are evaluated. Those are stored in a Field defined on a BatchedInterpolationIdxRange. Note that the coordinates of the points represented by this index range are unused and irrelevant. 
* `lagrange_coef` A Field storing the [**Lagrange**](classLagrange.md) coefficients. 




        

<hr>



### function operator() 

_Evaluate_ [_**Lagrange**_](classLagrange.md) _polynomials (described by their_[_**Lagrange**_](classLagrange.md) _coefficients) on a mesh._
```C++
template<class Layout1, class Layout2, class BatchedInterpolationGrid>
inline void LagrangeEvaluator::operator() (
    Field< DataType, BatchedInterpolationGrid, memory_space , Layout1 > const lagrange_eval,
    ConstField< DataType, batched_lagrange_domain_type < BatchedInterpolationGrid >, memory_space , Layout2 > const lagrange_coef
) const
```



The [**Lagrange**](classLagrange.md) coefficients describe [**Lagrange**](classLagrange.md) polynomials defined on a cartesian product of batch\_idx\_range and [**Lagrange**](classLagrange.md) basis. They can be obtained using the [**IdentityInterpolationBuilder**](classIdentityInterpolationBuilder.md).


This is not a multidimensional evaluation. This is a batched 1D evaluation. This means that for each slice of coordinates identified by a batch\_domain\_type::discrete\_element\_type, the evaluation is performed with the 1D set of [**Lagrange**](classLagrange.md) coefficients identified by the same batch\_domain\_type::discrete\_element\_type.




**Parameters:**


* `lagrange_eval` The values of the [**Lagrange**](classLagrange.md) polynomials at the mesh points. 
* `lagrange_coef` A Field storing the [**Lagrange**](classLagrange.md) coefficients. 




        

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



Extrapolation rules are functors used to define the behaviour of the [**LagrangeEvaluator**](classLagrangeEvaluator.md) outside the domain where the break points of the [**Lagrange**](classLagrange.md) basis are defined.




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

