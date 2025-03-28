

# File mask\_tanh.hpp



[**FileList**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**mask\_tanh.hpp**](mask__tanh_8hpp.md)

[Go to the source code of this file](mask__tanh_8hpp_source.md)



* `#include "ddc_aliases.hpp"`
* `#include "geometry.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| enum  | [**MaskType**](#enum-masktype)  <br>_An enum class that allows choosing between two types of masks._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  host\_t&lt; DFieldMemX &gt; | [**mask\_tanh**](#function-mask_tanh) (IdxRangeX const & gridx, double extent, double stiffness, MaskType const type, bool normalised) <br>_Constructs a mask function defined with hyperbolic tangents._  |




























## Public Types Documentation




### enum MaskType 

_An enum class that allows choosing between two types of masks._ 
```C++
enum MaskType {
    Normal,
    Inverted
};
```




<hr>
## Public Functions Documentation




### function mask\_tanh 

_Constructs a mask function defined with hyperbolic tangents._ 
```C++
host_t< DFieldMemX > mask_tanh (
    IdxRangeX const & gridx,
    double extent,
    double stiffness,
    MaskType const type,
    bool normalised
) 
```



Consider the index range [xmin, xmax], and {xleft, xright} the transition coordinates defined using the extent parameter.


If type = 'normal' the mask equals one inside the [xleft, xright] interval and zero outside.


If type = 'inverted' the mask equals zero inside the [xleft, xright] interval and one outside.


If normalised = true, the mask is normalised so that its integral equals one




**Parameters:**


* `gridx` The mesh in the x direction. 
* `extent` A parameter that sets the extent of the mask. 
* `stiffness` A parameter that sets the stiffness of the mask. 
* `type` A MaskType parameter that defines the type of the mask. 
* `normalised` A boolean that equals true if the integral of the mask must be equal to one. 



**Returns:**

A Dfield containing the mask. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/mask_tanh.hpp`

