

# Struct R\_cov



[**ClassList**](annotated.md) **>** [**R\_cov**](structR__cov.md)



_Define non periodic real covariant_ [_**R**_](structR.md) _dimension._

* `#include <geometry.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**R**](structR.md) | [**Dual**](#typedef-dual)  <br>_A type-alias mapping to the covariant counterpart._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  bool constexpr | [**IS\_CONTRAVARIANT**](#variable-is_contravariant)   = `false`<br>_A boolean indicating if dimension describes a contravariant coordinate._  |
|  bool constexpr | [**IS\_COVARIANT**](#variable-is_covariant)   = `true`<br>_A boolean indicating if dimension describes a covariant coordinate._  |
|  bool constexpr | [**PERIODIC**](#variable-periodic)   = `false`<br>_Define periodicity of the dimension. Here, not periodic._  |










































## Public Types Documentation




### typedef Dual 

_A type-alias mapping to the covariant counterpart._ 
```C++
using R_cov::Dual =  R;
```




<hr>
## Public Static Attributes Documentation




### variable IS\_CONTRAVARIANT 

_A boolean indicating if dimension describes a contravariant coordinate._ 
```C++
bool constexpr R_cov::IS_CONTRAVARIANT;
```




<hr>



### variable IS\_COVARIANT 

_A boolean indicating if dimension describes a covariant coordinate._ 
```C++
bool constexpr R_cov::IS_COVARIANT;
```




<hr>



### variable PERIODIC 

_Define periodicity of the dimension. Here, not periodic._ 
```C++
bool constexpr R_cov::PERIODIC;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/geometry/geometry.hpp`

