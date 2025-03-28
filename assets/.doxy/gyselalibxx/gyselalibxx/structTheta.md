

# Struct Theta



[**ClassList**](annotated.md) **>** [**Theta**](structTheta.md)



_Define periodic real contravariant_ [_**Theta**_](structTheta.md) _dimension._

* `#include <geometry.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Theta\_cov**](structTheta__cov.md) | [**Dual**](#typedef-dual)  <br>_A type-alias mapping to the covariant counterpart._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  bool constexpr | [**IS\_CONTRAVARIANT**](#variable-is_contravariant)   = `true`<br>_A boolean indicating if dimension describes a contravariant coordinate._  |
|  bool constexpr | [**IS\_COVARIANT**](#variable-is_covariant)   = `false`<br>_A boolean indicating if dimension describes a covariant coordinate._  |
|  bool constexpr | [**PERIODIC**](#variable-periodic)   = `true`<br>_Define periodicity of the dimension. Here, periodic._  |










































## Public Types Documentation




### typedef Dual 

_A type-alias mapping to the covariant counterpart._ 
```C++
using Theta::Dual =  Theta_cov;
```




<hr>
## Public Static Attributes Documentation




### variable IS\_CONTRAVARIANT 

_A boolean indicating if dimension describes a contravariant coordinate._ 
```C++
bool constexpr Theta::IS_CONTRAVARIANT;
```




<hr>



### variable IS\_COVARIANT 

_A boolean indicating if dimension describes a covariant coordinate._ 
```C++
bool constexpr Theta::IS_COVARIANT;
```




<hr>



### variable PERIODIC 

_Define periodicity of the dimension. Here, periodic._ 
```C++
bool constexpr Theta::PERIODIC;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/geometry/geometry.hpp`

