

# Class Curl

**template &lt;concepts::Mapping Mapping3D, class MappingCoord&gt;**



[**ClassList**](annotated.md) **>** [**Curl**](classCurl.md)



_A class which implements a curl operator The implemented equation is:_ \(\nabla \cross \mathbf{F}\) __\(\nabla \cross \mathbf{F} = {\cal J}_{\rm x}^{-1}\epsilon^{klm}\partial_{x^l} F_m\) _with_\({\cal J}_{\rm x}\) _the jacobian of the system,_\(F_m\) _the covariant components of F and_\(\epsilon^{klm}\) _the Levi-Civita symbol._[More...](#detailed-description)

* `#include <curl.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Curl**](#function-curl) (Mapping3D const & mapping) <br>_Build a curl operator._  |
|  KOKKOS\_INLINE\_FUNCTION [**Tensor**](classTensor.md)&lt; T, get\_contravariant\_dims\_t&lt; CovBasisF &gt; &gt; | [**operator()**](#function-operator) ([**Tensor**](classTensor.md)&lt; T, CovBasisF, CovBasisSpatial &gt; const & partial\_derivatives\_f, MappingCoord const & coord) const<br>_Compute a curl at a given coordinate, from the partial derivatives of the vector field f._  |




























## Detailed Description




**Template parameters:**


* `Mapping3D` A type representing a mapping in 3 dimensions. 
* `MappingCoord` The type of the coordinate that will be used to evaluate the mapping. This coordinate is used to calculate the determinant of the Jacobian and the metric tensor. It is almost always equal to the argument type of the 3D mapping but axi-symmetry may make it useful to evaluate the mapping on a coordinate with a reduced dimensionality. 




    
## Public Functions Documentation




### function Curl 

_Build a curl operator._ 
```C++
inline explicit Curl::Curl (
    Mapping3D const & mapping
) 
```





**Parameters:**


* `mapping` The mapping describing the system of coordinates on which the expression is calculated. 




        

<hr>



### function operator() 

_Compute a curl at a given coordinate, from the partial derivatives of the vector field f._ 
```C++
template<class T, class CovBasisF>
inline KOKKOS_INLINE_FUNCTION Tensor < T, get_contravariant_dims_t< CovBasisF > > Curl::operator() (
    Tensor < T, CovBasisF, CovBasisSpatial > const & partial_derivatives_f,
    MappingCoord const & coord
) const
```





**Parameters:**


* `partial_derivatives_f` A tensor containing the partial derivatives of the vector field f expressed at the given coordinate. 
* `coord` The coordinate where the calculation is carried out. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/curl.hpp`

