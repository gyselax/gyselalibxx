

# Class Moments::Impl

**template &lt;class Grid1D, class MemorySpace&gt;**



[**ClassList**](annotated.md) **>** [**Moments**](classMoments.md) **>** [**Impl**](classMoments_1_1Impl.md)



[_**Impl**_](classMoments_1_1Impl.md) _object storing attributes in_`MemorySpace` _._

* `#include <moments.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**Moments**](classMoments.md) | [**discrete\_dimension\_type**](#typedef-discrete_dimension_type)  <br>_alias of the discrete dimension_  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Impl**](#function-impl-23) ([**Impl**](classMoments_1_1Impl.md)&lt; Grid1D, OMemorySpace &gt; const & impl) <br>_Conversion constructor between different memory spaces._  |
|   | [**Impl**](#function-impl-33) () <br>_Main constructor._  |




























## Public Types Documentation




### typedef discrete\_dimension\_type 

_alias of the discrete dimension_ 
```C++
using Moments::Impl< Grid1D, MemorySpace >::discrete_dimension_type =  Moments;
```




<hr>
## Public Functions Documentation




### function Impl [2/3]

_Conversion constructor between different memory spaces._ 
```C++
template<class OMemorySpace>
inline explicit Moments::Impl::Impl (
    Impl < Grid1D, OMemorySpace > const & impl
) 
```





**Parameters:**


* `impl` object from `OMemorySpace` that will be used to initialise this object on `MemorySpace` 




        

<hr>



### function Impl [3/3]

_Main constructor._ 
```C++
inline Moments::Impl::Impl () 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/speciesinfo/moments.hpp`

