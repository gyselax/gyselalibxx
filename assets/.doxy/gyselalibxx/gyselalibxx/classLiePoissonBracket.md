

# Class LiePoissonBracket

**template &lt;class Mapping3D, class MappingCoord&gt;**



[**ClassList**](annotated.md) **>** [**LiePoissonBracket**](classLiePoissonBracket.md)



_A class which implements a gyrokinetic Poisson bracket operator. The implemented equation is:_  _with_ _the unitary magnetic field, i.e:_ _with_ _the jacobian of the system,_ _the covariant components of b and_ _the Levi-Civita symbol._[More...](#detailed-description)

* `#include <lie_poisson_bracket.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**LiePoissonBracket**](#function-liepoissonbracket) (Mapping3D const & mapping) <br>_Build a_ [_**LiePoissonBracket**_](classLiePoissonBracket.md) _operator._ |
|  KOKKOS\_INLINE\_FUNCTION double | [**operator()**](#function-operator) ([**DTensor**](classTensor.md)&lt; CovBasisSpatial &gt; const & partial\_derivatives\_f, [**DTensor**](classTensor.md)&lt; CovBasisSpatial &gt; const & partial\_derivatives\_g, [**DTensor**](classTensor.md)&lt; BasisSpatial &gt; const & B, MappingCoord const & coord) const<br>_Compute the gyrokinetic Poisson bracket at a given coordinate, from the partial derivatives of the two fields, and the magnetic field._  |
|  void | [**operator()**](#function-operator_1) (ExecSpace exec\_space, DField&lt; IdxRange, MemorySpace &gt; poisson\_bracket, [**DVectorConstField**](classVectorField.md)&lt; IdxRange, CovBasisSpatial, MemorySpace &gt; const partial\_derivatives\_f, [**DVectorConstField**](classVectorField.md)&lt; IdxRange, CovBasisSpatial, MemorySpace &gt; const partial\_derivatives\_g, [**DVectorConstField**](classVectorField.md)&lt; IdxRange, BasisSpatial, MemorySpace &gt; const B) <br>_Compute the gyrokinetic Poisson bracket at every point on a grid from the partial derivatives of the fields f and g, and the magnetic field._  |




























## Detailed Description




**Template parameters:**


* `Mapping3D` A type representing a mapping in 3 dimensions. 
* `MappingCoord` The type of the coordinate that will be used to evaluate the mapping. This coordinate is used to calculate the determinant of the Jacobian and the metric tensor. It is almost always equal to the argument type of the 3D mapping but axi-symmetry may make it useful to evaluate the mapping on a coordinate with a reduced dimensionality. 




    
## Public Functions Documentation




### function LiePoissonBracket 

_Build a_ [_**LiePoissonBracket**_](classLiePoissonBracket.md) _operator._
```C++
inline explicit LiePoissonBracket::LiePoissonBracket (
    Mapping3D const & mapping
) 
```





**Parameters:**


* `mapping` The mapping describing the system of coordinates on which the expression is calculated. 




        

<hr>



### function operator() 

_Compute the gyrokinetic Poisson bracket at a given coordinate, from the partial derivatives of the two fields, and the magnetic field._ 
```C++
inline KOKKOS_INLINE_FUNCTION double LiePoissonBracket::operator() (
    DTensor < CovBasisSpatial > const & partial_derivatives_f,
    DTensor < CovBasisSpatial > const & partial_derivatives_g,
    DTensor < BasisSpatial > const & B,
    MappingCoord const & coord
) const
```





**Parameters:**


* `partial_derivatives_f` A vector containing the partial derivatives of the scalar field f expressed at the given coordinate. 
* `partial_derivatives_g` A vector containing the partial derivatives of the scalar field g expressed at the given coordinate. 
* `B` A vector containing the magnetic field at the given coordinate. 
* `coord` The coordinate where the calculation is carried out. 




        

<hr>



### function operator() 

_Compute the gyrokinetic Poisson bracket at every point on a grid from the partial derivatives of the fields f and g, and the magnetic field._ 
```C++
template<class ExecSpace, class IdxRange, class MemorySpace>
inline void LiePoissonBracket::operator() (
    ExecSpace exec_space,
    DField< IdxRange, MemorySpace > poisson_bracket,
    DVectorConstField < IdxRange, CovBasisSpatial, MemorySpace > const partial_derivatives_f,
    DVectorConstField < IdxRange, CovBasisSpatial, MemorySpace > const partial_derivatives_g,
    DVectorConstField < IdxRange, BasisSpatial, MemorySpace > const B
) 
```





**Parameters:**


* `exec_space` The space (CPU/GPU) where the calculation should be executed. 
* `poisson_bracket` The result of the calculation of the gyrokinetic Poisson bracket at every point on a grid. 
* `partial_derivatives_f` A vector field containing the partial derivatives of f at each point on the grid. 
* `partial_derivatives_g` A vector field containing the partial derivatives of f at each point on the grid. 
* `B` A vector field describing the magnetic field. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/lie_poisson_bracket.hpp`

