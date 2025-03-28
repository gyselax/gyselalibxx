

# Class VectorMapper&lt; VectorIndexSet&lt; XIn, YIn &gt;, VectorIndexSet&lt; XOut, YOut &gt;, Mapping, ExecSpace &gt;

**template &lt;class XIn, class YIn, class XOut, class YOut, class Mapping, class ExecSpace&gt;**



[**ClassList**](annotated.md) **>** [**VectorMapper&lt; VectorIndexSet&lt; XIn, YIn &gt;, VectorIndexSet&lt; XOut, YOut &gt;, Mapping, ExecSpace &gt;**](classVectorMapper_3_01VectorIndexSet_3_01XIn_00_01YIn_01_4_00_01VectorIndexSet_3_01XOut_00_01YOu77c12468788509067d2c0ef34f5e389c.md)



_A class to map vector fields from one coordinate system to another._ [More...](#detailed-description)

* `#include <vector_mapper.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename ExecSpace::memory\_space | [**memory\_space**](#typedef-memory_space)  <br>_The type of the memory space where the field is saved (CPU vs GPU)._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**VectorMapper**](#function-vectormapper) (Mapping mapping) <br>_A constructor for the_ [_**VectorMapper**_](classVectorMapper.md) _._ |
|  void | [**operator()**](#function-operator) (ExecSpace exec\_space, [**VectorField**](classVectorField.md)&lt; double, IdxRangeType, VectorIndexSet&lt; XOut, YOut &gt;, [**memory\_space**](classVectorMapper_3_01VectorIndexSet_3_01XIn_00_01YIn_01_4_00_01VectorIndexSet_3_01XOut_00_01YOu77c12468788509067d2c0ef34f5e389c.md#typedef-memory_space), LayoutStridedPolicy1 &gt; vector\_field\_output, [**VectorConstField**](classVectorField.md)&lt; double, IdxRangeType, VectorIndexSet&lt; XIn, YIn &gt;, [**memory\_space**](classVectorMapper_3_01VectorIndexSet_3_01XIn_00_01YIn_01_4_00_01VectorIndexSet_3_01XOut_00_01YOu77c12468788509067d2c0ef34f5e389c.md#typedef-memory_space), LayoutStridedPolicy2 &gt; vector\_field\_input) <br>_Convert vectors defined in the input coordinate system to equivalent vectors in the output coordinate system._  |




























## Detailed Description


 

**Template parameters:**


* `InVectorSpace` A VectorIndexSet&lt;XIn, YIn&gt; describing the dimensions of the coordinate system taken as input. 
* `OutVectorSpace` A VectorIndexSet&lt;XOut, YOut&gt; describing the dimensions of the coordinate system returned as output. 
* `Mapping` A class describing a mapping system. 
* `ExecSpace` The space (CPU/GPU) where the calculations are carried out. 




    
## Public Types Documentation




### typedef memory\_space 

_The type of the memory space where the field is saved (CPU vs GPU)._ 
```C++
using VectorMapper< VectorIndexSet< XIn, YIn >, VectorIndexSet< XOut, YOut >, Mapping, ExecSpace >::memory_space =  typename ExecSpace::memory_space;
```




<hr>
## Public Functions Documentation




### function VectorMapper 

_A constructor for the_ [_**VectorMapper**_](classVectorMapper.md) _._
```C++
inline explicit VectorMapper< VectorIndexSet< XIn, YIn >, VectorIndexSet< XOut, YOut >, Mapping, ExecSpace >::VectorMapper (
    Mapping mapping
) 
```





**Parameters:**


* `mapping` The mapping description. 




        

<hr>



### function operator() 

_Convert vectors defined in the input coordinate system to equivalent vectors in the output coordinate system._ 
```C++
template<class IdxRangeType, class LayoutStridedPolicy1, class LayoutStridedPolicy2>
inline void VectorMapper< VectorIndexSet< XIn, YIn >, VectorIndexSet< XOut, YOut >, Mapping, ExecSpace >::operator() (
    ExecSpace exec_space,
    VectorField < double, IdxRangeType, VectorIndexSet< XOut, YOut >, memory_space , LayoutStridedPolicy1 > vector_field_output,
    VectorConstField < double, IdxRangeType, VectorIndexSet< XIn, YIn >, memory_space , LayoutStridedPolicy2 > vector_field_input
) 
```





**Parameters:**


* `exec_space` The space on which the function is executed (CPU/GPU). 
* `vector_field_output` The vector field containing the vectors in the output coordinate system. 
* `vector_field_input` The vector field containing the vectors in the input coordinate system. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/vector_mapper.hpp`

