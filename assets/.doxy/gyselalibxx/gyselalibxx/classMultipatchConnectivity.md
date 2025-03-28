

# Class MultipatchConnectivity

**template &lt;class... Interfaces&gt;**



[**ClassList**](annotated.md) **>** [**MultipatchConnectivity**](classMultipatchConnectivity.md)



_A helper class which provides functionalities to recognise how different patches are connected._ [More...](#detailed-description)

* `#include <connectivity.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef extract\_patches\_t&lt; [**inner\_edges**](classMultipatchConnectivity.md#typedef-inner_edges) &gt; | [**all\_patches**](#typedef-all_patches)  <br>_A type sequence of all the patches handled by this class._  |
| typedef ddcHelper::type\_seq\_intersection\_t&lt; [**get\_type\_seq\_connections\_t**](classMultipatchConnectivity.md#typedef-get_type_seq_connections_t)&lt; QueryPatch1 &gt;, [**get\_type\_seq\_connections\_t**](classMultipatchConnectivity.md#typedef-get_type_seq_connections_t)&lt; QueryPatch2 &gt; &gt; | [**find\_connections\_t**](#typedef-find_connections_t)  <br>_A tool to find all interfaces linking 2 patches._  |
| typedef collect\_interfaces\_on\_dim\_t&lt; find\_patch\_t&lt; Grid1D, [**all\_patches**](classMultipatchConnectivity.md#typedef-all_patches) &gt;, Grid1D, [**interface\_collection**](classMultipatchConnectivity.md#typedef-interface_collection) &gt; | [**get\_all\_interfaces\_along\_direction\_t**](#typedef-get_all_interfaces_along_direction_t)  <br>_A tool to find all interfaces along a line which passes through the requested grid._  |
| typedef to\_tuple\_t&lt; [**get\_type\_seq\_connections\_t**](classMultipatchConnectivity.md#typedef-get_type_seq_connections_t)&lt; QueryPatch &gt; &gt; | [**get\_connections\_t**](#typedef-get_connections_t)  <br>_A tool to find a tuple of all interfaces of a patch._  |
| typedef interfaces\_of\_patch\_t&lt; QueryPatch, [**interface\_collection**](classMultipatchConnectivity.md#typedef-interface_collection) &gt; | [**get\_type\_seq\_connections\_t**](#typedef-get_type_seq_connections_t)  <br>_A tool to find a type sequence of all interfaces of a patch._  |
| typedef strip\_outside\_edges\_t&lt; typename Interfaces::Edge1..., typename Interfaces::Edge2... &gt; | [**inner\_edges**](#typedef-inner_edges)  <br>_A type sequence of all the edges handled by this class except the OuterEdge types._  |
| typedef ddc::detail::TypeSeq&lt; Interfaces... &gt; | [**interface\_collection**](#typedef-interface_collection)  <br>_A type sequence of all the interfaces handled by this class._  |






















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  auto | [**get\_all\_idx\_ranges\_along\_direction**](#function-get_all_idx_ranges_along_direction-12) (std::tuple&lt; IdxRanges... &gt; all\_idx\_ranges) <br>_A function to return all index ranges which can be used to obtain coordinates along a line which passes through the requested grid._  |
|  auto | [**get\_all\_idx\_ranges\_along\_direction**](#function-get_all_idx_ranges_along_direction-22) ([**MultipatchType**](classMultipatchType.md)&lt; [**T**](structT.md), Patches... &gt; all\_idx\_ranges) <br>_A function to return all index ranges which can be used to obtain coordinates along a line which passes through the requested grid._  |


























## Detailed Description




**Template parameters:**


* `Interfaces` The interfaces which describe the connectivity in this geometry. 




    
## Public Types Documentation




### typedef all\_patches 

_A type sequence of all the patches handled by this class._ 
```C++
using MultipatchConnectivity< Interfaces >::all_patches =  extract_patches_t<inner_edges>;
```




<hr>



### typedef find\_connections\_t 

_A tool to find all interfaces linking 2 patches._ 
```C++
using MultipatchConnectivity< Interfaces >::find_connections_t =  ddcHelper::type_seq_intersection_t< get_type_seq_connections_t<QueryPatch1>, get_type_seq_connections_t<QueryPatch2> >;
```





**Template parameters:**


* `QueryPatch1` The first patch. 
* `QueryPatch2` The second patch. 




        

<hr>



### typedef get\_all\_interfaces\_along\_direction\_t 

_A tool to find all interfaces along a line which passes through the requested grid._ 
```C++
using MultipatchConnectivity< Interfaces >::get_all_interfaces_along_direction_t =  collect_interfaces_on_dim_t< find_patch_t<Grid1D, all_patches>, Grid1D, interface_collection>;
```





**Template parameters:**


* `Grid1D` The grid indicating the direction of interest. 




        

<hr>



### typedef get\_connections\_t 

_A tool to find a tuple of all interfaces of a patch._ 
```C++
using MultipatchConnectivity< Interfaces >::get_connections_t =  to_tuple_t<get_type_seq_connections_t<QueryPatch> >;
```





**Template parameters:**


* `QueryPatch` The patch whose connections we are interested in. 




        

<hr>



### typedef get\_type\_seq\_connections\_t 

_A tool to find a type sequence of all interfaces of a patch._ 
```C++
using MultipatchConnectivity< Interfaces >::get_type_seq_connections_t =  interfaces_of_patch_t<QueryPatch, interface_collection>;
```





**Template parameters:**


* `QueryPatch` The patch whose connections we are interested in. 




        

<hr>



### typedef inner\_edges 

_A type sequence of all the edges handled by this class except the OuterEdge types._ 
```C++
using MultipatchConnectivity< Interfaces >::inner_edges =  strip_outside_edges_t<typename Interfaces::Edge1..., typename Interfaces::Edge2...>;
```




<hr>



### typedef interface\_collection 

_A type sequence of all the interfaces handled by this class._ 
```C++
using MultipatchConnectivity< Interfaces >::interface_collection =  ddc::detail::TypeSeq<Interfaces...>;
```




<hr>
## Public Static Functions Documentation




### function get\_all\_idx\_ranges\_along\_direction [1/2]

_A function to return all index ranges which can be used to obtain coordinates along a line which passes through the requested grid._ 
```C++
template<class Grid1D, class... IdxRanges>
static inline auto MultipatchConnectivity::get_all_idx_ranges_along_direction (
    std::tuple< IdxRanges... > all_idx_ranges
) 
```





**Template parameters:**


* `Grid1D` The grid indicating the direction of interest. 



**Parameters:**


* `all_idx_ranges` A tuple containing all available index ranges.



**Returns:**

A tuple of index ranges along the line of interest. 





        

<hr>



### function get\_all\_idx\_ranges\_along\_direction [2/2]

_A function to return all index ranges which can be used to obtain coordinates along a line which passes through the requested grid._ 
```C++
template<class Grid1D, template< typename P > typename T, class... Patches>
static inline auto MultipatchConnectivity::get_all_idx_ranges_along_direction (
    MultipatchType < T , Patches... > all_idx_ranges
) 
```





**Template parameters:**


* `Grid1D` The grid indicating the direction of interest. 



**Parameters:**


* `all_idx_ranges` A [**MultipatchType**](classMultipatchType.md) containing all available index ranges.



**Returns:**

A tuple of index ranges along the line of interest. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/connectivity.hpp`

