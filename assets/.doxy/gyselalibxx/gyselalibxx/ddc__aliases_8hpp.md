

# File ddc\_aliases.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**ddc\_aliases.hpp**](ddc__aliases_8hpp.md)

[Go to the source code of this file](ddc__aliases_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::ChunkView&lt; ElementType, IdxRange, LayoutStridedPolicy, MemorySpace &gt; | [**ConstField**](#typedef-constfield)  <br>_An alias describing the type of a constant field defined on a grid (e.g. the electric field defined on the grid_  _)_ |
| typedef ddc::Coordinate&lt; Dims... &gt; | [**Coord**](#typedef-coord)  <br>_An alias describing the type of a coordinate (e.g. a coordinate in phase-space (x, vx))._  |
| typedef ConstField&lt; double, IdxRange, MemorySpace, LayoutStridedPolicy &gt; | [**DConstField**](#typedef-dconstfield)  <br>_An alias describing the type of a constant field of doubles defined on a grid (e.g. the electric field defined on the grid_  _)_ |
| typedef Field&lt; double, IdxRange, MemorySpace, LayoutStridedPolicy &gt; | [**DField**](#typedef-dfield)  <br>_An alias describing the type of a field of doubles defined on a grid (e.g. the electric field defined on the grid_  _)_ |
| typedef FieldMem&lt; double, IdxRange, MemSpace &gt; | [**DFieldMem**](#typedef-dfieldmem)  <br>_An alias describing the type of an object which will allocate memory for a field of doubles when it is created._  |
| typedef ddc::ChunkSpan&lt; ElementType, IdxRange, LayoutStridedPolicy, MemorySpace &gt; | [**Field**](#typedef-field)  <br>_An alias describing the type of a field defined on a grid (e.g. the electric field defined on the grid_  _)_ |
| typedef ddc::Chunk&lt; ElementType, IdxRange, ddc::KokkosAllocator&lt; ElementType, MemSpace &gt; &gt; | [**FieldMem**](#typedef-fieldmem)  <br>_An alias describing the type of an object which will allocate memory for a field when it is created._  |
| typedef ddc::DiscreteElement&lt; GridTypes... &gt; | [**Idx**](#typedef-idx)  <br>_An alias describing the type of an index that is used to access the values of a field defined on a grid._  |
| typedef ddc::DiscreteDomain&lt; GridTypes... &gt; | [**IdxRange**](#typedef-idxrange)  <br>_An alias describing the type of an index range describing the subsection of a grid on which a field is defined._  |
| typedef ddc::DiscreteVector&lt; GridTypes... &gt; | [**IdxStep**](#typedef-idxstep)  <br>_An alias describing the type of a distance between two indexes._  |
| typedef ddc::NonUniformPointSampling&lt; Dim &gt; | [**NonUniformGridBase**](#typedef-nonuniformgridbase)  <br>_An alias describing the type from which a non-uniform grid must inherit._  |
| typedef ddc::UniformPointSampling&lt; Dim &gt; | [**UniformGridBase**](#typedef-uniformgridbase)  <br>_An alias describing the type from which a uniform grid must inherit._  |
| typedef typename detail::GetSplineBatchedBuilder1D&lt; ExecSpace, BSplines, InterpolationGrid, BCLBound, BCUBound, SolverType, IdxRangeType &gt;::type | [**get\_spline\_batched\_builder1d\_t**](#typedef-get_spline_batched_builder1d_t)  <br>_An alias for creating a 1D spline builder based on GetSplineBatchedBuilder1D class._  |
| typedef typename detail::GetSplineBatchedBuilder2D&lt; ExecSpace, BSplinesX1, BSplinesX2, InterpolationGridX1, InterpolationGridX2, BCLBoundX1, BCUBoundX1, BCLBoundX2, BCUBoundX2, SolverType, IdxRangeType &gt;::type | [**get\_spline\_batched\_builder2d\_t**](#typedef-get_spline_batched_builder2d_t)  <br>_An alias for creating a 2D spline builder based on GetSplineBatchedBuilder1D class._  |
| typedef typename detail::GetSplineBatchedEvaluator1D&lt; ExecSpace, BSplines, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule, IdxRangeType &gt;::type | [**get\_spline\_batched\_evaluator1d\_t**](#typedef-get_spline_batched_evaluator1d_t)  <br>_An alias for creating a 1D spline evaluator based on GetSplineBatchedEvaluator1D class._  |
| typedef typename detail::GetSplineBatchedEvaluator2D&lt; ExecSpace, BSplinesX1, BSplinesX2, InterpolationGridX1, InterpolationGridX2, LowerExtrapolationRuleX1, UpperExtrapolationRuleX1, LowerExtrapolationRuleX2, UpperExtrapolationRuleX2, IdxRangeType &gt;::type | [**get\_spline\_batched\_evaluator2d\_t**](#typedef-get_spline_batched_evaluator2d_t)  <br>_An alias for creating a 2D spline evaluator based on GetSplineBatchedEvaluator2D class._  |
















































## Public Types Documentation




### typedef ConstField 

_An alias describing the type of a constant field defined on a grid (e.g. the electric field defined on the grid_  _)_
```C++
using ConstField =  ddc::ChunkView<ElementType, IdxRange, LayoutStridedPolicy, MemorySpace>;
```




<hr>



### typedef Coord 

_An alias describing the type of a coordinate (e.g. a coordinate in phase-space (x, vx))._ 
```C++
using Coord =  ddc::Coordinate<Dims...>;
```



This file contains aliases for DDC. The documentation for DDC can be found at [https://ddc.mdls.fr/](https://ddc.mdls.fr/). The names used in the DDC project are not always intuitive for mathematicians/physicists therefore Gysela has chosen to introduce this file to provide names that are hopefully more intuitive. The documentation for these concepts can be found in the Gysela documentation at: [https://gyselax.github.io/gyselalibxx/docs\_DDC\_in\_gyselalibxx.html](https://gyselax.github.io/gyselalibxx/docs_DDC_in_gyselalibxx.html) 


        

<hr>



### typedef DConstField 

_An alias describing the type of a constant field of doubles defined on a grid (e.g. the electric field defined on the grid_  _)_
```C++
using DConstField =  ConstField<double, IdxRange, MemorySpace, LayoutStridedPolicy>;
```




<hr>



### typedef DField 

_An alias describing the type of a field of doubles defined on a grid (e.g. the electric field defined on the grid_  _)_
```C++
using DField =  Field<double, IdxRange, MemorySpace, LayoutStridedPolicy>;
```




<hr>



### typedef DFieldMem 

_An alias describing the type of an object which will allocate memory for a field of doubles when it is created._ 
```C++
using DFieldMem =  FieldMem<double, IdxRange, MemSpace>;
```




<hr>



### typedef Field 

_An alias describing the type of a field defined on a grid (e.g. the electric field defined on the grid_  _)_
```C++
using Field =  ddc::ChunkSpan<ElementType, IdxRange, LayoutStridedPolicy, MemorySpace>;
```




<hr>



### typedef FieldMem 

_An alias describing the type of an object which will allocate memory for a field when it is created._ 
```C++
using FieldMem =  ddc::Chunk<ElementType, IdxRange, ddc::KokkosAllocator<ElementType, MemSpace> >;
```




<hr>



### typedef Idx 

_An alias describing the type of an index that is used to access the values of a field defined on a grid._ 
```C++
using Idx =  ddc::DiscreteElement<GridTypes...>;
```




<hr>



### typedef IdxRange 

_An alias describing the type of an index range describing the subsection of a grid on which a field is defined._ 
```C++
using IdxRange =  ddc::DiscreteDomain<GridTypes...>;
```




<hr>



### typedef IdxStep 

_An alias describing the type of a distance between two indexes._ 
```C++
using IdxStep =  ddc::DiscreteVector<GridTypes...>;
```




<hr>



### typedef NonUniformGridBase 

_An alias describing the type from which a non-uniform grid must inherit._ 
```C++
using NonUniformGridBase =  ddc::NonUniformPointSampling<Dim>;
```




<hr>



### typedef UniformGridBase 

_An alias describing the type from which a uniform grid must inherit._ 
```C++
using UniformGridBase =  ddc::UniformPointSampling<Dim>;
```




<hr>



### typedef get\_spline\_batched\_builder1d\_t 

_An alias for creating a 1D spline builder based on GetSplineBatchedBuilder1D class._ 
```C++
using get_spline_batched_builder1d_t =  typename detail::GetSplineBatchedBuilder1D< ExecSpace, BSplines, InterpolationGrid, BCLBound, BCUBound, SolverType, IdxRangeType>::type;
```




<hr>



### typedef get\_spline\_batched\_builder2d\_t 

_An alias for creating a 2D spline builder based on GetSplineBatchedBuilder1D class._ 
```C++
using get_spline_batched_builder2d_t =  typename detail::GetSplineBatchedBuilder2D< ExecSpace, BSplinesX1, BSplinesX2, InterpolationGridX1, InterpolationGridX2, BCLBoundX1, BCUBoundX1, BCLBoundX2, BCUBoundX2, SolverType, IdxRangeType>::type;
```




<hr>



### typedef get\_spline\_batched\_evaluator1d\_t 

_An alias for creating a 1D spline evaluator based on GetSplineBatchedEvaluator1D class._ 
```C++
using get_spline_batched_evaluator1d_t =  typename detail::GetSplineBatchedEvaluator1D< ExecSpace, BSplines, InterpolationGrid, LowerExtrapolationRule, UpperExtrapolationRule, IdxRangeType>::type;
```




<hr>



### typedef get\_spline\_batched\_evaluator2d\_t 

_An alias for creating a 2D spline evaluator based on GetSplineBatchedEvaluator2D class._ 
```C++
using get_spline_batched_evaluator2d_t =  typename detail::GetSplineBatchedEvaluator2D< ExecSpace, BSplinesX1, BSplinesX2, InterpolationGridX1, InterpolationGridX2, LowerExtrapolationRuleX1, UpperExtrapolationRuleX1, LowerExtrapolationRuleX2, UpperExtrapolationRuleX2, IdxRangeType>::type;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/ddc_aliases.hpp`

