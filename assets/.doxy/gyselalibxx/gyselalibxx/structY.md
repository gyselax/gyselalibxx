

# Struct Y



[**ClassList**](annotated.md) **>** [**Y**](structY.md)



_Define non periodic real_ [_**Y**_](structY.md) _dimension._[More...](#detailed-description)

* `#include <geometry.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Y**](structY.md) | [**Dual**](#typedef-dual-13)  <br>_A type-alias mapping to the covariant counterpart._  |
| typedef [**Y**](structY.md) | [**Dual**](#typedef-dual-13)  <br>_A type-alias mapping to the co/contra-variant counterpart._  |
| typedef [**Y**](structY.md) | [**Dual**](#typedef-dual-13)  <br>_A type-alias mapping to the co/contra-variant counterpart._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  bool constexpr | [**IS\_CONTRAVARIANT**](#variable-is_contravariant)   = `true`<br>_A boolean indicating if dimension describes a contravariant coordinate._  |
|  bool constexpr | [**IS\_COVARIANT**](#variable-is_covariant)   = `true`<br>_A boolean indicating if dimension describes a covariant coordinate._  |
|  bool constexpr | [**PERIODIC**](#variable-periodic)   = `false`<br>_Define periodicity of the dimension. Here, not periodic._  |










































## Detailed Description


A class which describes the real space in the second spatial direction [**Y**](structY.md). 


    
## Public Types Documentation




### typedef Dual [1/3]

_A type-alias mapping to the covariant counterpart._ 
```C++
using Y::Dual =  Y;
```




<hr>



### typedef Dual [1/3]

_A type-alias mapping to the co/contra-variant counterpart._ 
```C++
using Y::Dual =  Y;
```




<hr>



### typedef Dual [1/3]

_A type-alias mapping to the co/contra-variant counterpart._ 
```C++
using Y::Dual =  Y;
```




<hr>
## Public Static Attributes Documentation




### variable IS\_CONTRAVARIANT 

_A boolean indicating if dimension describes a contravariant coordinate._ 
```C++
static bool constexpr Y::IS_CONTRAVARIANT;
```




<hr>



### variable IS\_COVARIANT 

_A boolean indicating if dimension describes a covariant coordinate._ 
```C++
static bool constexpr Y::IS_COVARIANT;
```




<hr>



### variable PERIODIC 

_Define periodicity of the dimension. Here, not periodic._ 
```C++
static bool constexpr Y::PERIODIC;
```



A boolean indicating if the dimension is periodic. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/geometry/geometry.hpp`

