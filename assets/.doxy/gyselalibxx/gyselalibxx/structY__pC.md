

# Struct Y\_pC



[**ClassList**](annotated.md) **>** [**Y\_pC**](structY__pC.md)



_Tag the second non periodic dimension in the pseudo\_Cartesian index range._ 

* `#include <geometry_pseudo_cartesian.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Y\_pC**](structY__pC.md) | [**Dual**](#typedef-dual)  <br>_A type-alias mapping to the co/contra-variant counterpart._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  bool constexpr | [**IS\_CONTRAVARIANT**](#variable-is_contravariant)   = `true`<br>_A boolean indicating if dimension describes a contravariant coordinate._  |
|  bool constexpr | [**IS\_COVARIANT**](#variable-is_covariant)   = `true`<br>_A boolean indicating if dimension describes a covariant coordinate._  |
|  bool constexpr | [**PERIODIC**](#variable-periodic)   = `false`<br>_Define periodicity of the dimension. Here, not periodic._  |










































## Public Types Documentation




### typedef Dual 

_A type-alias mapping to the co/contra-variant counterpart._ 
```C++
using Y_pC::Dual =  Y_pC;
```




<hr>
## Public Static Attributes Documentation




### variable IS\_CONTRAVARIANT 

_A boolean indicating if dimension describes a contravariant coordinate._ 
```C++
bool constexpr Y_pC::IS_CONTRAVARIANT;
```




<hr>



### variable IS\_COVARIANT 

_A boolean indicating if dimension describes a covariant coordinate._ 
```C++
bool constexpr Y_pC::IS_COVARIANT;
```




<hr>



### variable PERIODIC 

_Define periodicity of the dimension. Here, not periodic._ 
```C++
bool constexpr Y_pC::PERIODIC;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mapping/geometry_pseudo_cartesian.hpp`

