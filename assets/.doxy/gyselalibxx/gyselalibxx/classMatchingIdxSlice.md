

# Class MatchingIdxSlice

**template &lt;class [**Interface**](structInterface.md)&gt;**



[**ClassList**](annotated.md) **>** [**MatchingIdxSlice**](classMatchingIdxSlice.md)



_Store the conforming indexes of each patch of a given interface._ [More...](#detailed-description)

* `#include <matching_idx_slice.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MatchingIdxSlice**](#function-matchingidxslice-12) (IdxRange1D\_1 const & idx\_range\_1, IdxRange1D\_2 const & idx\_range\_2) <br>_Instantiate the class from 1D index ranges._  |
|   | [**MatchingIdxSlice**](#function-matchingidxslice-22) (IdxRange2D\_1 const & idx\_range\_1, IdxRange2D\_2 const & idx\_range\_2) <br>_Instantiate the class from 2D index ranges._  |
|  [**IdxRangeSlice**](classIdxRangeSlice.md)&lt; ParallelGrid &gt; | [**get**](#function-get) () const<br>_Get the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _containing the conforming indexes._ |
|  auto | [**get\_from\_perp**](#function-get_from_perp) () const<br>_Get the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _containing the conforming indexes._ |
|   | [**~MatchingIdxSlice**](#function-matchingidxslice) () = default<br> |




























## Detailed Description


The conforming indexes are the indexes with an equivalent index on the parallel grid of the other edge of a given interface. The conforming indexes are stored in an [**IdxRangeSlice**](classIdxRangeSlice.md). The index step between two indexes in the [**IdxRangeSlice**](classIdxRangeSlice.md) are supposed to be uniform. If they are not, the instantiation of the class fails.


If the grids are uniform, the index steps between the conforming indexes are uniform. The uniform index steps can be deduced from the greatest common divisor.


E.g. for a first grid of N cells and a second grid of M cells, gcd = gcd(M, N) idx\_step on the first grid: N / gcd idx\_step on the second grid: M / gcd
* 






**Template parameters:**


* [**Interface**](structInterface.md) [**Interface**](structInterface.md) type between two edges of patches ([**Interface**](structInterface.md) with [**OutsideEdge**](structOutsideEdge.md) not allowed). 




    
## Public Functions Documentation




### function MatchingIdxSlice [1/2]

_Instantiate the class from 1D index ranges._ 
```C++
inline MatchingIdxSlice::MatchingIdxSlice (
    IdxRange1D_1 const & idx_range_1,
    IdxRange1D_2 const & idx_range_2
) 
```



To define the IdxRangeSlices containing the conforming indexes, we first check that the index steps between two conforming indexes are uniform, for the 1D grid of each edge of the interface.


If true, the index step of the patch is applied to instantiate the associated [**IdxRangeSlice**](classIdxRangeSlice.md).


If the grids are uniform, it is true and we can use the greatest common divisor between the two number of cells to compute the index steps of each [**IdxRangeSlice**](classIdxRangeSlice.md).




**Parameters:**


* `idx_range_1` 1D IdxRange of the first [**Edge**](structEdge.md) of the [**Interface**](structInterface.md). 
* `idx_range_2` 1D IdxRange of the second [**Edge**](structEdge.md) of the [**Interface**](structInterface.md). 




        

<hr>



### function MatchingIdxSlice [2/2]

_Instantiate the class from 2D index ranges._ 
```C++
inline MatchingIdxSlice::MatchingIdxSlice (
    IdxRange2D_1 const & idx_range_1,
    IdxRange2D_2 const & idx_range_2
) 
```



To define the IdxRangeSlices containing the conforming indexes, we first check that the index steps between two conforming indexes are uniform, for the 1D grid of each edge of the interface.


If true, the index step of the patch is applied to instantiate the associated [**IdxRangeSlice**](classIdxRangeSlice.md).


If the grids are uniform, it is true and we can use the greatest common divisor between the two number of cells to compute the index steps of each [**IdxRangeSlice**](classIdxRangeSlice.md).




**Parameters:**


* `idx_range_1` 2D IdxRange of the first [**Edge**](structEdge.md) of the [**Interface**](structInterface.md). 
* `idx_range_2` 2D IdxRange of the second [**Edge**](structEdge.md) of the [**Interface**](structInterface.md). 




        

<hr>



### function get 

_Get the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _containing the conforming indexes._
```C++
template<class ParallelGrid, std::enable_if_t<(std::is_same_v< ParallelGrid, EdgeGrid1 >)||(std::is_same_v< ParallelGrid, EdgeGrid2 >), bool >>
inline IdxRangeSlice < ParallelGrid > MatchingIdxSlice::get () const
```





**Template parameters:**


* `ParallelGrid` The parallel grid to the edge. 



**Returns:**

[**IdxRangeSlice**](classIdxRangeSlice.md) of conforming indexes on the given ParallelGrid. 





        

<hr>



### function get\_from\_perp 

_Get the_ [_**IdxRangeSlice**_](classIdxRangeSlice.md) _containing the conforming indexes._
```C++
template<class PerpendicularGrid, std::enable_if_t<(std::is_same_v< PerpendicularGrid, PerpEdgeGrid1 >)||(std::is_same_v< PerpendicularGrid, PerpEdgeGrid2 >), bool >>
inline auto MatchingIdxSlice::get_from_perp () const
```





**Template parameters:**


* `PerpendicularGrid` The perpendicular grid to the edge. 



**Returns:**

[**IdxRangeSlice**](classIdxRangeSlice.md) of conforming indexes on the perpendicular grid to the given PerpendicularGrid. 





        

<hr>



### function ~MatchingIdxSlice 

```C++
MatchingIdxSlice::~MatchingIdxSlice () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/connectivity/matching_idx_slice.hpp`

