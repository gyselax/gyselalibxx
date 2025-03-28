

# Struct connectivity\_details::CollectAllGridsOnDim

**template &lt;class StartPatch, class Grid1D, class InterfaceTypeSeq&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**CollectAllGridsOnDim**](structconnectivity__details_1_1CollectAllGridsOnDim.md)



_A class which collects all grids along a given dimension in both directions._ [More...](#detailed-description)

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**CollectGridsAlongDim**](structconnectivity__details_1_1CollectGridsAlongDim.md)&lt; [**Edge**](structEdge.md)&lt; StartPatch, Grid1D, FRONT &gt;, InterfaceTypeSeq, BackInsert &gt;[**::type**](structconnectivity__details_1_1CollectAllGridsOnDim.md#typedef-type) | [**BackwardTypeSeq**](#typedef-backwardtypeseq)  <br>_The type sequence describing all grids found by iterating along this dimension in the backwards direction._  |
| typedef ddc::type\_seq\_merge\_t&lt; [**BackwardTypeSeq**](structconnectivity__details_1_1CollectAllGridsOnDim.md#typedef-backwardtypeseq), typename [**CollectGridsAlongDim**](structconnectivity__details_1_1CollectGridsAlongDim.md)&lt; [**Edge**](structEdge.md)&lt; StartPatch, Grid1D, BACK &gt;, InterfaceTypeSeq, FrontInsert, [**BackwardTypeSeq**](structconnectivity__details_1_1CollectAllGridsOnDim.md#typedef-backwardtypeseq) &gt;::type &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Detailed Description




**Template parameters:**


* `StartPatch` The patch from which the collection should begin. 
* `Grid1D` The first grid to be included (this describes the dimension along which grids are collected). 
* `InterfaceTypeSeq` A DDC type sequence containing all the possible Interfaces. 




    
## Public Types Documentation




### typedef BackwardTypeSeq 

_The type sequence describing all grids found by iterating along this dimension in the backwards direction._ 
```C++
using connectivity_details::CollectAllGridsOnDim< StartPatch, Grid1D, InterfaceTypeSeq >::BackwardTypeSeq =  typename CollectGridsAlongDim< Edge<StartPatch, Grid1D, FRONT>, InterfaceTypeSeq, BackInsert>::type;
```



This is found by working backward from front (start) of grid inserting each new grid at the start of the sequence. 


        

<hr>



### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::CollectAllGridsOnDim< StartPatch, Grid1D, InterfaceTypeSeq >::type =  ddc::type_seq_merge_t< BackwardTypeSeq, typename CollectGridsAlongDim< Edge<StartPatch, Grid1D, BACK>, InterfaceTypeSeq, FrontInsert, BackwardTypeSeq>::type>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

