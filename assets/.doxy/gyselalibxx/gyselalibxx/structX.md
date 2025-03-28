

# Struct X



[**ClassList**](annotated.md) **>** [**X**](structX.md)



_Define non periodic real_ [_**X**_](structX.md) _dimension._[More...](#detailed-description)

* `#include <geometry.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**X**](structX.md) | [**Dual**](#typedef-dual-14)  <br>_A type-alias mapping to the covariant counterpart._  |
| typedef [**X**](structX.md) | [**Dual**](#typedef-dual-14)  <br>_A type-alias mapping to the co/contra-variant counterpart._  |
| typedef [**X**](structX.md) | [**Dual**](#typedef-dual-14)  <br>_A type-alias mapping to the co/contra-variant counterpart._  |
| typedef [**X**](structX.md) | [**Dual**](#typedef-dual-14)  <br>_A type-alias mapping to the co/contra-variant counterpart._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  bool constexpr | [**IS\_CONTRAVARIANT**](#variable-is_contravariant)   = `true`<br>_A boolean indicating if dimension describes a contravariant coordinate._  |
|  bool constexpr | [**IS\_COVARIANT**](#variable-is_covariant)   = `true`<br>_A boolean indicating if dimension describes a covariant coordinate._  |
|  bool constexpr | [**PERIODIC**](#variable-periodic)   = `false`<br>_Define periodicity of the dimension. Here, not periodic._  |










































## Detailed Description


A class which describes the real space in the first spatial direction [**X**](structX.md).


A class which describes the real space in the spatial [**X**](structX.md) direction. 


    
## Public Types Documentation




### typedef Dual [1/4]

_A type-alias mapping to the covariant counterpart._ 
```C++
using X::Dual =  X;
```




<hr>



### typedef Dual [1/4]

_A type-alias mapping to the co/contra-variant counterpart._ 
```C++
using X::Dual =  X;
```




<hr>



### typedef Dual [1/4]

_A type-alias mapping to the co/contra-variant counterpart._ 
```C++
using X::Dual =  X;
```




<hr>



### typedef Dual [1/4]

_A type-alias mapping to the co/contra-variant counterpart._ 
```C++
using X::Dual =  X;
```




<hr>
## Public Static Attributes Documentation




### variable IS\_CONTRAVARIANT 

_A boolean indicating if dimension describes a contravariant coordinate._ 
```C++
static bool constexpr X::IS_CONTRAVARIANT;
```




<hr>



### variable IS\_COVARIANT 

_A boolean indicating if dimension describes a covariant coordinate._ 
```C++
static bool constexpr X::IS_COVARIANT;
```




<hr>



### variable PERIODIC 

_Define periodicity of the dimension. Here, not periodic._ 
```C++
static bool constexpr X::PERIODIC;
```



A boolean indicating if the dimension is periodic. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/geometry/geometry.hpp`

