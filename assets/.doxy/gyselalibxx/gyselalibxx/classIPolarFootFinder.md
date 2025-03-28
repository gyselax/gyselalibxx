

# Class IPolarFootFinder

**template &lt;class GridRadial, class GridPoloidal, class AdvectionDim1, class AdvectionDim2, class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**IPolarFootFinder**](classIPolarFootFinder.md)



_Define a base class for all the time integration methods used to find the foot of a characteristic on a polar domain (a polar domain is a domain defined on the_  _plane)._[More...](#detailed-description)

* `#include <ipolar_foot_finder.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](#function-operator) (Field&lt; Coord&lt; [**R**](classIPolarFootFinder.md#typedef-r), [**Theta**](classIPolarFootFinder.md#typedef-theta) &gt;, [**IdxRangeRTheta**](classIPolarFootFinder.md#typedef-idxrangertheta), [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space) &gt; feet, [**DVectorConstField**](classVectorField.md)&lt; [**IdxRangeRTheta**](classIPolarFootFinder.md#typedef-idxrangertheta), VectorIndexSet&lt; [**X**](classIPolarFootFinder.md#typedef-x), [**Y**](classIPolarFootFinder.md#typedef-y) &gt;, [**memory\_space**](classIPolarFootFinder.md#typedef-memory_space) &gt; advection\_field, double dt) const = 0<br>_Advect the feet over_  _._ |
| virtual  | [**~IPolarFootFinder**](#function-ipolarfootfinder) () = default<br> |




## Protected Types

| Type | Name |
| ---: | :--- |
| typedef GridRadial | [**GridR**](#typedef-gridr)  <br>_The continuous radial dimension._  |
| typedef GridPoloidal | [**GridTheta**](#typedef-gridtheta)  <br>_The continuous poloidal dimension._  |
| typedef IdxRange&lt; [**GridR**](classIPolarFootFinder.md#typedef-gridr), [**GridTheta**](classIPolarFootFinder.md#typedef-gridtheta) &gt; | [**IdxRangeRTheta**](#typedef-idxrangertheta)  <br>_The type of the index range over which the operator works._  |
| typedef typename GridR::continuous\_dimension\_type | [**R**](#typedef-r)  <br>_The continuous radial dimension._  |
| typedef typename GridTheta::continuous\_dimension\_type | [**Theta**](#typedef-theta)  <br>_The continuous poloidal dimension._  |
| typedef AdvectionDim1 | [**X**](#typedef-x)  <br>_The continuous radial dimension._  |
| typedef AdvectionDim2 | [**Y**](#typedef-y)  <br>_The continuous poloidal dimension._  |
| typedef MemorySpace | [**memory\_space**](#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |
























## Detailed Description


)




**Template parameters:**


* `GridRadial` The radial grid on which the distribution function is defined. 
* `GridPoloidal` The poloidial grid on which the distribution function is defined. 
* `AdvectionDim1` The first dimension of the advection field vector. 
* `AdvectionDim2` The second dimension of the advection field vector. 
* `MemorySpace` The memory space where the data is saved (CPU/GPU). 




    
## Public Functions Documentation




### function operator() 

_Advect the feet over_  _._
```C++
virtual void IPolarFootFinder::operator() (
    Field< Coord< R , Theta >, IdxRangeRTheta , memory_space > feet,
    DVectorConstField < IdxRangeRTheta , VectorIndexSet< X , Y >, memory_space > advection_field,
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
using IPolarFootFinder< GridRadial, GridPoloidal, AdvectionDim1, AdvectionDim2, MemorySpace >::GridR =  GridRadial;
```




<hr>



### typedef GridTheta 

_The continuous poloidal dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, AdvectionDim1, AdvectionDim2, MemorySpace >::GridTheta =  GridPoloidal;
```




<hr>



### typedef IdxRangeRTheta 

_The type of the index range over which the operator works._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, AdvectionDim1, AdvectionDim2, MemorySpace >::IdxRangeRTheta =  IdxRange<GridR, GridTheta>;
```




<hr>



### typedef R 

_The continuous radial dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, AdvectionDim1, AdvectionDim2, MemorySpace >::R =  typename GridR::continuous_dimension_type;
```




<hr>



### typedef Theta 

_The continuous poloidal dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, AdvectionDim1, AdvectionDim2, MemorySpace >::Theta =  typename GridTheta::continuous_dimension_type;
```




<hr>



### typedef X 

_The continuous radial dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, AdvectionDim1, AdvectionDim2, MemorySpace >::X =  AdvectionDim1;
```




<hr>



### typedef Y 

_The continuous poloidal dimension._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, AdvectionDim1, AdvectionDim2, MemorySpace >::Y =  AdvectionDim2;
```




<hr>



### typedef memory\_space 

_The type of the memory space where the field is saved (CPU vs GPU)._ 
```C++
using IPolarFootFinder< GridRadial, GridPoloidal, AdvectionDim1, AdvectionDim2, MemorySpace >::memory_space =  MemorySpace;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/ipolar_foot_finder.hpp`

