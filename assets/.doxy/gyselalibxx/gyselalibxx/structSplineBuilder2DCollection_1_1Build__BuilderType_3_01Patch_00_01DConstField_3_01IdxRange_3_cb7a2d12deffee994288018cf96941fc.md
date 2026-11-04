

# Struct SplineBuilder2DCollection::Build\_BuilderType&lt; Patch, DConstField&lt; IdxRange&lt; Grid1D... &gt;, MemorySpace &gt; &gt;

**template &lt;class [**Patch**](structPatch.md), class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**Build\_BuilderType&lt; Patch, DConstField&lt; IdxRange&lt; Grid1D... &gt;, MemorySpace &gt; &gt;**](structSplineBuilder2DCollection_1_1Build__BuilderType_3_01Patch_00_01DConstField_3_01IdxRange_3_cb7a2d12deffee994288018cf96941fc.md)






















## Public Types

| Type | Name |
| ---: | :--- |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), Grid1OnPatch&lt; [**Patch**](structPatch.md) &gt;, FRONT &gt;, typename Connectivity::interface\_collection &gt; | [**lower\_matching\_edge1**](#typedef-lower_matching_edge1)  <br> |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), Grid2OnPatch&lt; [**Patch**](structPatch.md) &gt;, FRONT &gt;, typename Connectivity::interface\_collection &gt; | [**lower\_matching\_edge2**](#typedef-lower_matching_edge2)  <br> |
| typedef ddc::SplineBuilder2D&lt; ExecSpace, MemorySpace, BSpline1OnPatch&lt; [**Patch**](structPatch.md) &gt;, BSpline2OnPatch&lt; [**Patch**](structPatch.md) &gt;, Grid1OnPatch&lt; [**Patch**](structPatch.md) &gt;, Grid2OnPatch&lt; [**Patch**](structPatch.md) &gt;, std::is\_same\_v&lt; lower\_matching\_edge1, [**OutsideEdge**](structOutsideEdge.md) &gt; ? SBCLower1 :BcTransition, std::is\_same\_v&lt; upper\_matching\_edge1, [**OutsideEdge**](structOutsideEdge.md) &gt; ? SBCUpper1 :BcTransition, std::is\_same\_v&lt; lower\_matching\_edge2, [**OutsideEdge**](structOutsideEdge.md) &gt; ? SBCLower2 :BcTransition, std::is\_same\_v&lt; upper\_matching\_edge2, [**OutsideEdge**](structOutsideEdge.md) &gt; ? SBCUpper2 :BcTransition, Solver &gt; | [**type**](#typedef-type)  <br> |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), Grid1OnPatch&lt; [**Patch**](structPatch.md) &gt;, BACK &gt;, typename Connectivity::interface\_collection &gt; | [**upper\_matching\_edge1**](#typedef-upper_matching_edge1)  <br> |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), Grid2OnPatch&lt; [**Patch**](structPatch.md) &gt;, BACK &gt;, typename Connectivity::interface\_collection &gt; | [**upper\_matching\_edge2**](#typedef-upper_matching_edge2)  <br> |
















































## Public Types Documentation




### typedef lower\_matching\_edge1 

```C++
using SplineBuilder2DCollection< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::lower_matching_edge1 =  equivalent_edge_t<
                Edge<Patch, Grid1OnPatch<Patch>, FRONT>,
                typename Connectivity::interface_collection>;
```




<hr>



### typedef lower\_matching\_edge2 

```C++
using SplineBuilder2DCollection< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::lower_matching_edge2 =  equivalent_edge_t<
                Edge<Patch, Grid2OnPatch<Patch>, FRONT>,
                typename Connectivity::interface_collection>;
```




<hr>



### typedef type 

```C++
using SplineBuilder2DCollection< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::type =  ddc::SplineBuilder2D<
                ExecSpace,
                MemorySpace,
                BSpline1OnPatch<Patch>,
                BSpline2OnPatch<Patch>,
                Grid1OnPatch<Patch>,
                Grid2OnPatch<Patch>,
                std::is_same_v<lower_matching_edge1, OutsideEdge> ? SBCLower1 : BcTransition,
                std::is_same_v<upper_matching_edge1, OutsideEdge> ? SBCUpper1 : BcTransition,
                std::is_same_v<lower_matching_edge2, OutsideEdge> ? SBCLower2 : BcTransition,
                std::is_same_v<upper_matching_edge2, OutsideEdge> ? SBCUpper2 : BcTransition,
                Solver>;
```




<hr>



### typedef upper\_matching\_edge1 

```C++
using SplineBuilder2DCollection< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::upper_matching_edge1 =  equivalent_edge_t<
                Edge<Patch, Grid1OnPatch<Patch>, BACK>,
                typename Connectivity::interface_collection>;
```




<hr>



### typedef upper\_matching\_edge2 

```C++
using SplineBuilder2DCollection< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, SBCLower1, SBCUpper1, SBCLower2, SBCUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::upper_matching_edge2 =  equivalent_edge_t<
                Edge<Patch, Grid2OnPatch<Patch>, BACK>,
                typename Connectivity::interface_collection>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/spline_builder_2d_collection.hpp`

