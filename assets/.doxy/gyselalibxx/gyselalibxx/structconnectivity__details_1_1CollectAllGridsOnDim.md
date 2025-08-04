

# Struct connectivity\_details::CollectAllGridsOnDim

**template &lt;class StartPatch, class Grid1D, class InterfaceTypeSeq&gt;**



[**ClassList**](annotated.md) **>** [**connectivity\_details**](namespaceconnectivity__details.md) **>** [**CollectAllGridsOnDim**](structconnectivity__details_1_1CollectAllGridsOnDim.md)



_A class which collects all grids along a given dimension in both directions._ [More...](#detailed-description)

* `#include <connectivity_details.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::type\_seq\_merge\_t&lt; BackwardTypeSeq, ForwardTypeSeq &gt; | [**type**](#typedef-type)  <br>_The type found by the class._  |
















































## Detailed Description




**Template parameters:**


* `StartPatch` The patch from which the collection should begin. 
* `Grid1D` The first grid to be included (this describes the dimension along which grids are collected). 
* `InterfaceTypeSeq` A DDC type sequence containing all the possible Interfaces. 




    
## Public Types Documentation




### typedef type 

_The type found by the class._ 
```C++
using connectivity_details::CollectAllGridsOnDim< StartPatch, Grid1D, InterfaceTypeSeq >::type =  ddc::type_seq_merge_t<BackwardTypeSeq, ForwardTypeSeq>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity_details.hpp`

