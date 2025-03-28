

# Struct MultipatchSplineBuilder2D::Build\_BuilderType&lt; Patch, DConstField&lt; IdxRange&lt; Grid1D... &gt;, MemorySpace &gt; &gt;

**template &lt;class [**Patch**](structPatch.md), class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**Build\_BuilderType&lt; Patch, DConstField&lt; IdxRange&lt; Grid1D... &gt;, MemorySpace &gt; &gt;**](structMultipatchSplineBuilder2D_1_1Build__BuilderType_3_01Patch_00_01DConstField_3_01IdxRange_3_388990a8744187d12e0f612652c86727.md)






















## Public Types

| Type | Name |
| ---: | :--- |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), Grid1OnPatch&lt; [**Patch**](structPatch.md) &gt;, FRONT &gt;, typename Connectivity::interface\_collection &gt; | [**lower\_matching\_edge1**](#typedef-lower_matching_edge1)  <br> |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), Grid2OnPatch&lt; [**Patch**](structPatch.md) &gt;, FRONT &gt;, typename Connectivity::interface\_collection &gt; | [**lower\_matching\_edge2**](#typedef-lower_matching_edge2)  <br> |
| typedef ddc::SplineBuilder2D&lt; ExecSpace, MemorySpace, BSpline1OnPatch&lt; [**Patch**](structPatch.md) &gt;, BSpline2OnPatch&lt; [**Patch**](structPatch.md) &gt;, Grid1OnPatch&lt; [**Patch**](structPatch.md) &gt;, Grid2OnPatch&lt; [**Patch**](structPatch.md) &gt;, std::is\_same\_v&lt; lower\_matching\_edge1, [**OutsideEdge**](structOutsideEdge.md) &gt; ? BcLower1 :BcTransition, std::is\_same\_v&lt; upper\_matching\_edge1, [**OutsideEdge**](structOutsideEdge.md) &gt; ? BcUpper1 :BcTransition, std::is\_same\_v&lt; lower\_matching\_edge2, [**OutsideEdge**](structOutsideEdge.md) &gt; ? BcLower2 :BcTransition, std::is\_same\_v&lt; upper\_matching\_edge2, [**OutsideEdge**](structOutsideEdge.md) &gt; ? BcUpper2 :BcTransition, Solver, Grid1D... &gt; | [**type**](#typedef-type)  <br> |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), Grid1OnPatch&lt; [**Patch**](structPatch.md) &gt;, BACK &gt;, typename Connectivity::interface\_collection &gt; | [**upper\_matching\_edge1**](#typedef-upper_matching_edge1)  <br> |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), Grid2OnPatch&lt; [**Patch**](structPatch.md) &gt;, BACK &gt;, typename Connectivity::interface\_collection &gt; | [**upper\_matching\_edge2**](#typedef-upper_matching_edge2)  <br> |
















































## Public Types Documentation




### typedef lower\_matching\_edge1 

```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, BcLower1, BcUpper1, BcLower2, BcUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::lower_matching_edge1 =  equivalent_edge_t< Edge<Patch, Grid1OnPatch<Patch>, FRONT>, typename Connectivity::interface_collection>;
```




<hr>



### typedef lower\_matching\_edge2 

```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, BcLower1, BcUpper1, BcLower2, BcUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::lower_matching_edge2 =  equivalent_edge_t< Edge<Patch, Grid2OnPatch<Patch>, FRONT>, typename Connectivity::interface_collection>;
```




<hr>



### typedef type 

```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, BcLower1, BcUpper1, BcLower2, BcUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::type =  ddc::SplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch<Patch>, BSpline2OnPatch<Patch>, Grid1OnPatch<Patch>, Grid2OnPatch<Patch>, std::is_same_v<lower_matching_edge1, OutsideEdge> ? BcLower1 : BcTransition, std::is_same_v<upper_matching_edge1, OutsideEdge> ? BcUpper1 : BcTransition, std::is_same_v<lower_matching_edge2, OutsideEdge> ? BcLower2 : BcTransition, std::is_same_v<upper_matching_edge2, OutsideEdge> ? BcUpper2 : BcTransition, Solver, Grid1D...>;
```




<hr>



### typedef upper\_matching\_edge1 

```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, BcLower1, BcUpper1, BcLower2, BcUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::upper_matching_edge1 =  equivalent_edge_t< Edge<Patch, Grid1OnPatch<Patch>, BACK>, typename Connectivity::interface_collection>;
```




<hr>



### typedef upper\_matching\_edge2 

```C++
using MultipatchSplineBuilder2D< ExecSpace, MemorySpace, BSpline1OnPatch, BSpline2OnPatch, Grid1OnPatch, Grid2OnPatch, BcLower1, BcUpper1, BcLower2, BcUpper2, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::upper_matching_edge2 =  equivalent_edge_t< Edge<Patch, Grid2OnPatch<Patch>, BACK>, typename Connectivity::interface_collection>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/multipatch_spline_builder_2d.hpp`

