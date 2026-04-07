

# Class NDIdentityInterpolationBuilder&lt; ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange&lt; Basis... &gt; &gt;

**template &lt;class ExecSpace, class MemorySpace, class DataType, class IdxRangeInterpolation, class... Basis&gt;**



[**ClassList**](annotated.md) **>** [**NDIdentityInterpolationBuilder&lt; ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange&lt; Basis... &gt; &gt;**](classNDIdentityInterpolationBuilder_3_01ExecSpace_00_01MemorySpace_00_01DataType_00_01IdxRangeIn4bcac5b3fb4ccc82c859be6708c9c5ff.md)



_The implementation of_ [_**NDIdentityInterpolationBuilder**_](classNDIdentityInterpolationBuilder.md) _. This is separate to allow a variadic Basis._

* `#include <nd_identity_interpolation_builder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::detail::convert\_type\_seq\_to\_discrete\_domain\_t&lt; ddc::type\_seq\_replace\_t&lt; ddc::to\_type\_seq\_t&lt; BatchedInterpolationIdxRange &gt;, ddc::to\_type\_seq\_t&lt; IdxRangeInterpolation &gt;, ddc::detail::TypeSeq&lt; typename Basis::template Impl&lt; Basis, MemorySpace &gt;::knot\_grid... &gt; &gt; &gt; | [**batched\_basis\_idx\_range\_type**](#typedef-batched_basis_idx_range_type)  <br>_Batched domain with each interpolation grid replaced by its basis domain._  |
| typedef BatchedInterpolationIdxRange | [**batched\_derivs\_idx\_range\_type**](#typedef-batched_derivs_idx_range_type)  <br>_The type of the index range on which derivatives should be provided (here unused)._  |
| typedef IdxRange&lt; typename Basis::template Impl&lt; Basis, MemorySpace &gt;::knot\_grid... &gt; | [**coeff\_idx\_range\_type**](#typedef-coeff_idx_range_type)  <br>_The type of the index range for the bases over which coefficients of an ND Lagrange interpolation are defined._  |
| typedef DataType | [**data\_type**](#typedef-data_type)  <br>_The data type that the data is saved on._  |
| typedef ExecSpace | [**exec\_space**](#typedef-exec_space)  <br>_The type of the Kokkos execution space._  |
| typedef IdxRangeInterpolation | [**interpolation\_idx\_range\_type**](#typedef-interpolation_idx_range_type)  <br>_The ND index range for the interpolation mesh._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The type of the Kokkos memory space._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**NDIdentityInterpolationBuilder**](#function-ndidentityinterpolationbuilder) () = default<br> |
|  [**batched\_basis\_idx\_range\_type**](classNDIdentityInterpolationBuilder_3_01ExecSpace_00_01MemorySpace_00_01DataType_00_01IdxRangeIn4bcac5b3fb4ccc82c859be6708c9c5ff.md#typedef-batched_basis_idx_range_type)&lt; BatchedInterpolationIdxRange &gt; | [**batched\_basis\_idx\_range**](#function-batched_basis_idx_range) (BatchedInterpolationIdxRange const & batched\_interpolation\_domain) noexcept const<br>_Get the batched basis index range for a given batched interpolation domain._  |
|  void | [**operator()**](#function-operator) (Field&lt; DataType, [**batched\_basis\_idx\_range\_type**](classNDIdentityInterpolationBuilder_3_01ExecSpace_00_01MemorySpace_00_01DataType_00_01IdxRangeIn4bcac5b3fb4ccc82c859be6708c9c5ff.md#typedef-batched_basis_idx_range_type)&lt; BatchedInterpolationIdxRange &gt;, [**memory\_space**](classNDIdentityInterpolationBuilder_3_01ExecSpace_00_01MemorySpace_00_01DataType_00_01IdxRangeIn4bcac5b3fb4ccc82c859be6708c9c5ff.md#typedef-memory_space) &gt; coeffs, ConstField&lt; DataType, BatchedInterpolationIdxRange, [**memory\_space**](classNDIdentityInterpolationBuilder_3_01ExecSpace_00_01MemorySpace_00_01DataType_00_01IdxRangeIn4bcac5b3fb4ccc82c859be6708c9c5ff.md#typedef-memory_space) &gt; vals) const<br>_Compute the interpolation coefficients for a function._  |




























## Public Types Documentation




### typedef batched\_basis\_idx\_range\_type 

_Batched domain with each interpolation grid replaced by its basis domain._ 
```C++
using NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::batched_basis_idx_range_type =  ddc::detail::convert_type_seq_to_discrete_domain_t<ddc::type_seq_replace_t< ddc::to_type_seq_t<BatchedInterpolationIdxRange>, ddc::to_type_seq_t<IdxRangeInterpolation>, ddc::detail::TypeSeq< typename Basis::template Impl<Basis, MemorySpace>::knot_grid...> >>;
```



Chains ddc::replace\_dim\_of\_t for each (InterpolationGrid\_i, BasisDomain\_i) pair. 


        

<hr>



### typedef batched\_derivs\_idx\_range\_type 

_The type of the index range on which derivatives should be provided (here unused)._ 
```C++
using NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::batched_derivs_idx_range_type =  BatchedInterpolationIdxRange;
```




<hr>



### typedef coeff\_idx\_range\_type 

_The type of the index range for the bases over which coefficients of an ND Lagrange interpolation are defined._ 
```C++
using NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::coeff_idx_range_type =  IdxRange<typename Basis::template Impl<Basis, MemorySpace>::knot_grid...>;
```




<hr>



### typedef data\_type 

_The data type that the data is saved on._ 
```C++
using NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::data_type =  DataType;
```




<hr>



### typedef exec\_space 

_The type of the Kokkos execution space._ 
```C++
using NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::exec_space =  ExecSpace;
```




<hr>



### typedef interpolation\_idx\_range\_type 

_The ND index range for the interpolation mesh._ 
```C++
using NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::interpolation_idx_range_type =  IdxRangeInterpolation;
```




<hr>



### typedef memory\_space 

_The type of the Kokkos memory space._ 
```C++
using NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::memory_space =  MemorySpace;
```




<hr>
## Public Functions Documentation




### function NDIdentityInterpolationBuilder 

```C++
NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::NDIdentityInterpolationBuilder () = default
```




<hr>



### function batched\_basis\_idx\_range 

_Get the batched basis index range for a given batched interpolation domain._ 
```C++
template<class BatchedInterpolationIdxRange>
inline batched_basis_idx_range_type < BatchedInterpolationIdxRange > NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::batched_basis_idx_range (
    BatchedInterpolationIdxRange const & batched_interpolation_domain
) noexcept const
```





**Parameters:**


* `batched_interpolation_domain` The full batched interpolation domain. 



**Returns:**

The batched basis index range. 





        

<hr>



### function operator() 

_Compute the interpolation coefficients for a function._ 
```C++
template<class BatchedInterpolationIdxRange>
inline void NDIdentityInterpolationBuilder< ExecSpace, MemorySpace, DataType, IdxRangeInterpolation, IdxRange< Basis... > >::operator() (
    Field< DataType, batched_basis_idx_range_type < BatchedInterpolationIdxRange >, memory_space > coeffs,
    ConstField< DataType, BatchedInterpolationIdxRange, memory_space > vals
) const
```



Copies vals directly to coeffs. No computation is needed because the interpolation coefficients equal the function values at the grid nodes.




**Note:**

Periodic bases are not yet supported. The copy logic for periodic dimensions (wrapping the last point back to the start) is non-trivial in ND because each combination of periodic/non-periodic dimensions must be handled separately (2^k copies for k periodic dimensions). This will be implemented when needed.




**Parameters:**


* `coeffs` The coefficients of the interpolation. 
* `vals` The values of the function on the interpolation mesh. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/nd_identity_interpolation_builder.hpp`

