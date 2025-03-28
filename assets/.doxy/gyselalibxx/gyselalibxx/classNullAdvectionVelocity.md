

# Class NullAdvectionVelocity

**template &lt;class IdxRangeFdistribu, class IdxRangeSpatial&gt;**



[**ClassList**](annotated.md) **>** [**NullAdvectionVelocity**](classNullAdvectionVelocity.md)



_This is a class which imitates a velocity advection. It inherits from IAdvectionV and can be used as an advection operator but does not actually modify the distribution function. This can be useful for debugging purposes._ 

* `#include <nulladvectionvx.hpp>`



Inherits the following classes: IAdvectionV< IdxRangeFdistribu, IdxRangeSpatial >


































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**NullAdvectionVelocity**](#function-nulladvectionvelocity) () = default<br> |
|  DField&lt; IdxRangeFdistribu &gt; | [**operator()**](#function-operator) (DField&lt; IdxRangeFdistribu &gt; allfdistribu, DConstField&lt; IdxRangeSpatial &gt; electric\_field, double dt) override const<br>_Do nothing instead of advecting fdistribu along GridV for a duration dt._  |
|   | [**~NullAdvectionVelocity**](#function-nulladvectionvelocity) () override<br> |




























## Public Functions Documentation




### function NullAdvectionVelocity 

```C++
NullAdvectionVelocity::NullAdvectionVelocity () = default
```




<hr>



### function operator() 

_Do nothing instead of advecting fdistribu along GridV for a duration dt._ 
```C++
inline DField< IdxRangeFdistribu > NullAdvectionVelocity::operator() (
    DField< IdxRangeFdistribu > allfdistribu,
    DConstField< IdxRangeSpatial > electric_field,
    double dt
) override const
```





**Parameters:**


* `allfdistribu` Reference to an array containing the value of the distribution function. 
* `electric_field` The electric field which derives from electrostatic potential and is the advection speed. 
* `dt` Time step



**Returns:**

A reference to the allfdistribu array containing the value of the function at the coordinates. 





        

<hr>



### function ~NullAdvectionVelocity 

```C++
NullAdvectionVelocity::~NullAdvectionVelocity () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/nulladvectionvx.hpp`

