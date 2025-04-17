

# Struct MultipatchSplineBuilder::Build\_BuilderType

**template &lt;class [**Patch**](structPatch.md)&gt;**



[**ClassList**](annotated.md) **>** [**Build\_BuilderType**](structMultipatchSplineBuilder_1_1Build__BuilderType.md)



[More...](#detailed-description)


















## Public Types

| Type | Name |
| ---: | :--- |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), GridOnPatch&lt; [**Patch**](structPatch.md) &gt;, FRONT &gt;, typename Connectivity::interface\_collection &gt; | [**lower\_matching\_edge**](#typedef-lower_matching_edge)  <br> |
| typedef ddc::SplineBuilder&lt; ExecSpace, MemorySpace, BSplineOnPatch&lt; [**Patch**](structPatch.md) &gt;, GridOnPatch&lt; [**Patch**](structPatch.md) &gt;, std::is\_same\_v&lt; lower\_matching\_edge, [**OutsideEdge**](structOutsideEdge.md) &gt; ? BcLower :BcTransition, std::is\_same\_v&lt; upper\_matching\_edge, [**OutsideEdge**](structOutsideEdge.md) &gt; ? BcUpper :BcTransition, Solver &gt; | [**type**](#typedef-type)  <br> |
| typedef equivalent\_edge\_t&lt; [**Edge**](structEdge.md)&lt; [**Patch**](structPatch.md), GridOnPatch&lt; [**Patch**](structPatch.md) &gt;, BACK &gt;, typename Connectivity::interface\_collection &gt; | [**upper\_matching\_edge**](#typedef-upper_matching_edge)  <br> |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr std::size\_t | [**patch\_id**](#variable-patch_id)   = `ddc::type\_seq\_rank\_v&lt;[**Patch**](structPatch.md), PatchOrdering&gt;`<br> |










































## Detailed Description


A small structure allowing the multiple grids to be unpacked from a field and repacked into a SplineBuilder type. 


    
## Public Types Documentation




### typedef lower\_matching\_edge 

```C++
using MultipatchSplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch, GridOnPatch, BcLower, BcUpper, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch >::lower_matching_edge =  equivalent_edge_t< Edge<Patch, GridOnPatch<Patch>, FRONT>, typename Connectivity::interface_collection>;
```




<hr>



### typedef type 

```C++
using MultipatchSplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch, GridOnPatch, BcLower, BcUpper, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch >::type =  ddc::SplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch<Patch>, GridOnPatch<Patch>, std::is_same_v<lower_matching_edge, OutsideEdge> ? BcLower : BcTransition, std::is_same_v<upper_matching_edge, OutsideEdge> ? BcUpper : BcTransition, Solver>;
```




<hr>



### typedef upper\_matching\_edge 

```C++
using MultipatchSplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch, GridOnPatch, BcLower, BcUpper, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch >::upper_matching_edge =  equivalent_edge_t< Edge<Patch, GridOnPatch<Patch>, BACK>, typename Connectivity::interface_collection>;
```




<hr>
## Public Static Attributes Documentation




### variable patch\_id 

```C++
constexpr std::size_t MultipatchSplineBuilder< ExecSpace, MemorySpace, BSplineOnPatch, GridOnPatch, BcLower, BcUpper, BcTransition, Connectivity, Solver, ValuesOnPatch, Patches >::Build_BuilderType< Patch >::patch_id;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/spline/multipatch_spline_builder.hpp`

