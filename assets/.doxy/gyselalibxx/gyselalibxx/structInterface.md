

# Struct Interface

**template &lt;class Edge1Type, class Edge2Type, bool orientations\_agree\_bool&gt;**



[**ClassList**](annotated.md) **>** [**Interface**](structInterface.md)



_Represent a simple sticking of two edges._ [More...](#detailed-description)

* `#include <interface.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Edge1Type | [**Edge1**](#typedef-edge1)  <br>[_**Edge**_](structEdge.md) _type of the first patch._ |
| typedef Edge2Type | [**Edge2**](#typedef-edge2)  <br>[_**Edge**_](structEdge.md) _type of the second patch._ |
| typedef std::conditional\_t&lt; std::is\_same\_v&lt; [**Edge**](structEdge.md), [**Edge1**](structInterface.md#typedef-edge1) &gt;, [**Edge2**](structInterface.md#typedef-edge2), [**Edge1**](structInterface.md#typedef-edge1) &gt; | [**OtherEdge**](#typedef-otheredge)  <br>_A tool to get the other edge of an interface. I.e. to get Edge1 given Edge2 or Edge2 given Edge1._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**orientations\_agree**](#variable-orientations_agree)   = `orientations\_agree\_bool`<br> |
















## Public Static Functions

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**connected\_to\_patch**](#function-connected_to_patch) () <br>_A compile-time function to check if an interface is on the edge of a given patch._  |


























## Detailed Description




**Template parameters:**


* `Edge1Type` [**Edge**](structEdge.md) type defined on first patch. 
* `Edge2Type` [**Edge**](structEdge.md) type defined on second patch. 
* `orientations_agree_bool` Boolean, if true, the parametrisations of the physical edge have the same orientation; else, the parametrisations of the physical edge have the opposite orientation. 



**See also:** [**EdgeTransformation**](classEdgeTransformation.md) 



    
## Public Types Documentation




### typedef Edge1 

[_**Edge**_](structEdge.md) _type of the first patch._
```C++
using Interface< Edge1Type, Edge2Type, orientations_agree_bool >::Edge1 =  Edge1Type;
```




<hr>



### typedef Edge2 

[_**Edge**_](structEdge.md) _type of the second patch._
```C++
using Interface< Edge1Type, Edge2Type, orientations_agree_bool >::Edge2 =  Edge2Type;
```




<hr>



### typedef OtherEdge 

_A tool to get the other edge of an interface. I.e. to get Edge1 given Edge2 or Edge2 given Edge1._ 
```C++
using Interface< Edge1Type, Edge2Type, orientations_agree_bool >::OtherEdge =  std::conditional_t<std::is_same_v<Edge, Edge1>, Edge2, Edge1>;
```




<hr>
## Public Static Attributes Documentation




### variable orientations\_agree 

```C++
constexpr bool Interface< Edge1Type, Edge2Type, orientations_agree_bool >::orientations_agree;
```



If True, the parametrisations of the physical edge have the same orientation. If False, the parametrisations of the physical edge have the opposite orientation. (See [**EdgeTransformation**](classEdgeTransformation.md)). 


        

<hr>
## Public Static Functions Documentation




### function connected\_to\_patch 

_A compile-time function to check if an interface is on the edge of a given patch._ 
```C++
template<class Patch>
static inline constexpr bool Interface::connected_to_patch () 
```





**Template parameters:**


* `The` patch which may be connected 



**Returns:**

True if the patch is connected, False otherwise. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/interface.hpp`

