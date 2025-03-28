

# Struct MultipatchSplineBuilder::Build\_BuilderType&lt; Patch, DConstField&lt; IdxRange&lt; Grid1D... &gt;, MemorySpace &gt; &gt;

**template &lt;class [**Patch**](structPatch.md), class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**Build\_BuilderType&lt; Patch, DConstField&lt; IdxRange&lt; Grid1D... &gt;, MemorySpace &gt; &gt;**](structMultipatchSplineBuilder_1_1Build__BuilderType_3_01Patch_00_01DConstField_3_01IdxRange_3_0143963206948829b7e1024cedaaccd580.md)






















## Public Types

| Type | Name |
| ---: | :--- |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), GridOnPatch&lt; [**Patch**](structPatch.md) &gt;, FRONT &gt;, typename Connectivity::interface\_collection &gt; | [**lower\_matching\_edge**](#typedef-lower_matching_edge)  <br> |
| typedef ddc::SplineBuilder&lt; ExecSpace, MemorySpace, BSplineOnPatch&lt; [**Patch**](structPatch.md) &gt;, GridOnPatch&lt; [**Patch**](structPatch.md) &gt;, std::is\_same\_v&lt; lower\_matching\_edge, [**OutsideEdge**](structOutsideEdge.md) &gt; ? BcLower :BcTransition, std::is\_same\_v&lt; upper\_matching\_edge, [**OutsideEdge**](structOutsideEdge.md) &gt; ? BcUpper :BcTransition, Solver, Grid1D... &gt; | [**type**](#typedef-type)  <br> |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), GridOnPatch&lt; [**Patch**](structPatch.md) &gt;, BACK &gt;, typename Connectivity::interface\_collection &gt; | [**upper\_matching\_edge**](#typedef-upper_matching_edge)  <br> |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**patch\_id**](#variable-patch_id)   = `ddc::type\_seq\_rank\_v&lt;[**Patch**](structPatch.md), PatchOrdering&gt;`<br> |










































## Public Types Documentation




### typedef lower\_matching\_edge 

```C++
using MultipatchSplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch, GridOnPatch, BcLower, BcUpper, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::lower_matching_edge =  equivalent_edge_t< Edge<Patch, GridOnPatch<Patch>, FRONT>, typename Connectivity::interface_collection>;
```




<hr>



### typedef type 

```C++
using MultipatchSplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch, GridOnPatch, BcLower, BcUpper, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::type =  ddc::SplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch<Patch>, GridOnPatch<Patch>, std::is_same_v<lower_matching_edge, OutsideEdge> ? BcLower : BcTransition, std::is_same_v<upper_matching_edge, OutsideEdge> ? BcUpper : BcTransition, Solver, Grid1D...>;
```




<hr>



### typedef upper\_matching\_edge 

```C++
using MultipatchSplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch, GridOnPatch, BcLower, BcUpper, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::upper_matching_edge =  equivalent_edge_t< Edge<Patch, GridOnPatch<Patch>, BACK>, typename Connectivity::interface_collection>;
```




<hr>
## Public Static Attributes Documentation




### variable patch\_id 

```C++
constexpr std::size_t MultipatchSplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch, GridOnPatch, BcLower, BcUpper, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch, DConstField< IdxRange< Grid1D... >, MemorySpace > >::patch_id;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/multipatch_spline_builder.hpp`

