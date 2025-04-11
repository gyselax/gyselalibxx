

# Class EdgeTransformation

**template &lt;class [**Interface**](structInterface.md)&gt;**



[**ClassList**](annotated.md) **>** [**EdgeTransformation**](classEdgeTransformation.md)



_Transform a coordinate or an index from one edge to the one on the other edge._ [More...](#detailed-description)

* `#include <edge_transformation.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**EdgeTransformation**](#function-edgetransformation) (IdxRangeEdge1 const & idx\_range\_patch\_1, IdxRangeEdge2 const & idx\_range\_patch\_2) <br>_Instantiate an_ [_**EdgeTransformation**_](classEdgeTransformation.md) _._ |
|  bool | [**is\_match\_available**](#function-is_match_available) (CurrentIdx const & current\_idx) const<br>_Check if a given index has an equivalent index on the other patch of an interface._  |
|  Coord&lt; std::conditional\_t&lt; std::is\_same\_v&lt; CurrentDim, EdgeDim1 &gt;, EdgeDim2, EdgeDim1 &gt; &gt; | [**operator()**](#function-operator) (Coord&lt; CurrentDim &gt; const & current\_coord) const<br>_Transform a coordinate on the edge in the dimension of the current patch to the analogous coordinate on the target patch._  |
|  auto | [**operator()**](#function-operator_1) (CurrentIdx const & current\_idx) const<br>_Transform an index on the edge in the dimension of the current patch to the analogous index on the target patch._  |
|  bool | [**search\_for\_match**](#function-search_for_match) (Idx&lt; TargetGrid &gt; & target\_idx, Idx&lt; CurrentGrid &gt; current\_idx) const<br>_Check if a given index has an equivalent index and transform an index on the edge in the dimension of the current patch to the analogous index on the target patch._  |
|  Coord&lt; std::conditional\_t&lt; std::is\_same\_v&lt; CurrentPatch, Patch1 &gt;, EdgeDim2, EdgeDim1 &gt; &gt; | [**transform\_edge\_coord**](#function-transform_edge_coord) (Coord&lt; std::conditional\_t&lt; std::is\_same\_v&lt; CurrentPatch, Patch1 &gt;, EdgeDim1, EdgeDim2 &gt; &gt; const & current\_coord) const<br>_Transform a coordinate on the edge in the dimension of the current patch to the analogous coordinate on the target patch._  |
|   | [**~EdgeTransformation**](#function-edgetransformation) () = default<br> |




























## Detailed Description


According to the orientation of the interface, we compute the equivalent coordinate
* if True, 
* if False, 




where  and  are the minimum and maximum coordinates of the edge .


For the indices, we look for an equivalent index corresponding to a coordinate equivalent to the coordinate of the initial index.




**Template parameters:**


* [**Interface**](structInterface.md) The [**Interface**](structInterface.md) type where we want to compute the transformation. 




    
## Public Functions Documentation




### function EdgeTransformation 

_Instantiate an_ [_**EdgeTransformation**_](classEdgeTransformation.md) _._
```C++
inline EdgeTransformation::EdgeTransformation (
    IdxRangeEdge1 const & idx_range_patch_1,
    IdxRangeEdge2 const & idx_range_patch_2
) 
```





**Parameters:**


* `idx_range_patch_1` 1D index range on the patch 1 of the interface. 
* `idx_range_patch_2` 1D index range on the patch 2 of the interface. 




        

<hr>



### function is\_match\_available 

_Check if a given index has an equivalent index on the other patch of an interface._ 
```C++
template<class CurrentIdx>
inline bool EdgeTransformation::is_match_available (
    CurrentIdx const & current_idx
) const
```



If the grids are uniform, we can simplify the algorithm by using modulo. Otherwise, we need to check all the indices of the target grid. We suppose the coordinate transformation bijective, so we can use a dichotomy method. 



This method mainly calls search\_for\_match.




**Parameters:**


* `current_idx` A index on the edge of the current patch. 



**Template parameters:**


* `CurrentIdx` The current index type of the given coordinate index.



**Returns:**

Boolean stating if there is an equivalent index. 





        

<hr>



### function operator() 

_Transform a coordinate on the edge in the dimension of the current patch to the analogous coordinate on the target patch._ 
```C++
template<class CurrentDim>
inline Coord< std::conditional_t< std::is_same_v< CurrentDim, EdgeDim1 >, EdgeDim2, EdgeDim1 > > EdgeTransformation::operator() (
    Coord< CurrentDim > const & current_coord
) const
```





**Parameters:**


* `current_coord` A coordinate on the edge of the current patch.



**Template parameters:**


* `CurrentDim` The current continuous dimension of the given coordinate coord.



**Returns:**

The analogous coordinate on the target patch. 




**Warning:**

This operator is ill-defined when the two patches have the same continuous dimension. 





        

<hr>



### function operator() 

_Transform an index on the edge in the dimension of the current patch to the analogous index on the target patch._ 
```C++
template<class CurrentIdx>
inline auto EdgeTransformation::operator() (
    CurrentIdx const & current_idx
) const
```



If the grids are uniform, we can simplify the algorithm by using modulo. Otherwise, we need to check all the indices of the target grid. We suppose the coordinate transformation bijective, so we can use a dichotomy method.


This method mainly calls search\_for\_match.




**Parameters:**


* `current_idx` A index on the edge of the current patch. 



**Template parameters:**


* `CurrentIdx` The current index type of the given coordinate index.



**Returns:**

The analogous index on the target patch. 





        

<hr>



### function search\_for\_match 

_Check if a given index has an equivalent index and transform an index on the edge in the dimension of the current patch to the analogous index on the target patch._ 
```C++
template<class CurrentGrid, class TargetGrid>
inline bool EdgeTransformation::search_for_match (
    Idx< TargetGrid > & target_idx,
    Idx< CurrentGrid > current_idx
) const
```



If the grids are uniform, we can simplify the algorithm by using modulo. Otherwise, we need to check all the indices of the target grid. We suppose the coordinate transformation bijective, so we can use a dichotomy method.




**Warning:**

target\_idx is always replaced by the suspected index. If there is not equivalent index, the returned index is wrong but the closest that the algorithm found.




**Template parameters:**


* `CurrentGrid` The grid where the input index is defined. 
* `TargetGrid` The grid where the output index is defined.



**Parameters:**


* `target_idx` A index on the edge of the target patch. 
* `current_idx` A index on the edge of the current patch. 



**Template parameters:**


* `CurrentIdx` The current index type of the given coordinate index.



**Returns:**

Boolean stating if there is an equivalent index. 





        

<hr>



### function transform\_edge\_coord 

_Transform a coordinate on the edge in the dimension of the current patch to the analogous coordinate on the target patch._ 
```C++
template<class CurrentPatch>
inline Coord< std::conditional_t< std::is_same_v< CurrentPatch, Patch1 >, EdgeDim2, EdgeDim1 > > EdgeTransformation::transform_edge_coord (
    Coord< std::conditional_t< std::is_same_v< CurrentPatch, Patch1 >, EdgeDim1, EdgeDim2 > > const & current_coord
) const
```





**Parameters:**


* `current_coord` A coordinate on the edge of the current patch.



**Template parameters:**


* `CurrentPatch` The current patch of the given coordinate coord.



**Returns:**

The analogous coordinate on the target patch. 





        

<hr>



### function ~EdgeTransformation 

```C++
EdgeTransformation::~EdgeTransformation () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/edge_transformation.hpp`

