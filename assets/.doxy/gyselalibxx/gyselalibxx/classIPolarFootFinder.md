

# Class IPolarFootFinder

**template &lt;class GridRadial, class GridPoloidal, class VectorIndexSetAdvDims, class IdxRangeBatched, class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**IPolarFootFinder**](classIPolarFootFinder.md)



_Define a base class for all the time integration methods used to find the foot of a characteristic on a polar domain (a polar domain is a domain defined on the_  _plane)._[More...](#detailed-description)

* `#include <ipolar_foot_finder.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef IdxRangeBatched | [**IdxRangeOperator**](#typedef-idxrangeoperator)  <br>_The type of the index range over which the operator works._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](#function-operator) (Field&lt; Coord&lt; [**R**](classIPolarFootFinder.md#typedef-r), [**Theta**](classIPolarFootFinder.md#typedef-theta) &gt;, [**IdxRangeOperator**](classIPolarFootFinder.md#typedef-idxrangeoperator), [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space) &gt; feet, [**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeOperator**](classIPolarFootFinder.md#typedef-idxrangeoperator), [**VectorIndexSetAdvectionDims**](classIPolarFootFinder.md#typedef-vectorindexsetadvectiondims), [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space) &gt; advection\_field, double dt) const = 0<br>_Advect the feet over_  _._ |
| virtual  | [**~IPolarFootFinder**](#function-ipolarfootfinder) () = default<br> |




## Protected Types

| Type | Name |
| ---: | :--- |
| typedef GridRadial | [**GridR**](#typedef-gridr)  <br>_The continuous radial dimension._  |
| typedef GridPoloidal | [**GridTheta**](#typedef-gridtheta)  <br>_The continuous poloidal dimension._  |
| typedef typename GridR::continuous\_dimension\_type | [**R**](#typedef-r)  <br>_The continuous radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](#typedef-theta)  <br>_The continuous poloidal dimension._  |
| typedef VectorIndexSetAdvDims | [**VectorIndexSetAdvectionDims**](#typedef-vectorindexsetadvectiondims)  <br>_The continuous radial dimension._  |
























## Detailed Description


)




**Template parameters:**


* `GridRadial` The radial grid on which the distribution function is defined. 
* `GridPoloidal` The poloidial grid on which the distribution function is defined. 
* `VectorIndexSetAdvDims` A vector index set containing the set of dimensions that can be used to index the advection dimensions. 
* `IdxRangeBatched` The index range on which this batched operator operates. I.e. the index range of the distribution function. 
* `MemorySpace` The memory space where the data is saved (CPU/GPU). 




    
## Public Types Documentation




### typedef IdxRangeOperator 

_The type of the index range over which the operator works._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::IdxRangeOperator =  IdxRangeBatched;
```




<hr>



### typedef memory\_space 

_The type of the memory space where the field is saved (CPU vs GPU)._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::memory_space =  MemorySpace;
```




<hr>
## Public Functions Documentation




### function operator() 

_Advect the feet over_  _._
```C++
virtual void IPolarFootFinder::operator() (
    Field< Coord< R , Theta >, IdxRangeOperator , memory_space > feet,
    DVectorConstField < IdxRangeOperator , VectorIndexSetAdvectionDims , memory_space > advection_field,
    double dt
) const = 0
```





**Parameters:**


* `feet` On input: the mesh points. On output: the characteristic feet. 
* `advection_field` The advection field in the physical domain. 
* `dt` The time step. 




        

<hr>



### function ~IPolarFootFinder 

```C++
virtual IPolarFootFinder::~IPolarFootFinder () = default
```




<hr>
## Protected Types Documentation




### typedef GridR 

_The continuous radial dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::GridR =  GridRadial;
```




<hr>



### typedef GridTheta 

_The continuous poloidal dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::GridTheta =  GridPoloidal;
```




<hr>



### typedef R 

_The continuous radial dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef Theta 

_The continuous poloidal dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>



### typedef VectorIndexSetAdvectionDims 

_The continuous radial dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, VectorIndexSetAdvDims, IdxRangeBatched, MemorySpace >::VectorIndexSetAdvectionDims =  VectorIndexSetAdvDims;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/ipolar_foot_finder.hpp`

