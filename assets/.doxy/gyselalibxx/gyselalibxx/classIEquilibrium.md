

# Class IEquilibrium



[**ClassList**](annotated.md) **>** [**IEquilibrium**](classIEquilibrium.md)



_An abstract class for initialising a distribution function in (species,vpar,mu)._ [More...](#detailed-description)

* `#include <iequilibrium.hpp>`





Inherited by the following classes: [BumpontailEquilibrium](classBumpontailEquilibrium.md),  [MaxwellianEquilibrium](classMaxwellianEquilibrium.md),  [MaxwellianEquilibrium](classMaxwellianEquilibrium.md),  [MaxwellianEquilibrium](classMaxwellianEquilibrium.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVparMu | [**operator()**](#function-operator) (DFieldSpVparMu allfequilibrium) const = 0<br>_Operator for initialising an equilibrium distribution function._  |
| virtual DFieldSpVx | [**operator()**](#function-operator_1) (DFieldSpVx allfequilibrium) const = 0<br>_Operator for initialising a distribution function that does not depend on space._  |
| virtual DFieldSpVxVy | [**operator()**](#function-operator_2) (DFieldSpVxVy allfequilibrium) const = 0<br>_Operator for initialising a distribution function that does not depend on space._  |
| virtual  | [**~IEquilibrium**](#function-iequilibrium-13) () = default<br> |
| virtual  | [**~IEquilibrium**](#function-iequilibrium-13) () = default<br> |
| virtual  | [**~IEquilibrium**](#function-iequilibrium-13) () = default<br> |




























## Detailed Description


An abstract class for initialising the equilibrium state of the distribution function. The equilibrium state does not depend on spatial dimensions.


An abstract class for initialising a distribution function that does not depend on space. 


    
## Public Functions Documentation




### function operator() 

_Operator for initialising an equilibrium distribution function._ 
```C++
virtual DFieldSpVparMu IEquilibrium::operator() (
    DFieldSpVparMu allfequilibrium
) const = 0
```





**Parameters:**


* `allfequilibrium` On input: the uninitialized distribution function. On output: the initialised distribution function. 



**Returns:**

The initialised equilibrium distribution function. 





        

<hr>



### function operator() 

_Operator for initialising a distribution function that does not depend on space._ 
```C++
virtual DFieldSpVx IEquilibrium::operator() (
    DFieldSpVx allfequilibrium
) const = 0
```





**Parameters:**


* `allfequilibrium` On input: the uninitialized distribution function. On output: the initialised distribution function. 



**Returns:**

The initialised distribution function. 





        

<hr>



### function operator() 

_Operator for initialising a distribution function that does not depend on space._ 
```C++
virtual DFieldSpVxVy IEquilibrium::operator() (
    DFieldSpVxVy allfequilibrium
) const = 0
```





**Parameters:**


* `allfequilibrium` On input: the uninitialized distribution function. On output: the initialised distribution function. 



**Returns:**

The initialised distribution function. 





        

<hr>



### function ~IEquilibrium [1/3]

```C++
virtual IEquilibrium::~IEquilibrium () = default
```




<hr>



### function ~IEquilibrium [1/3]

```C++
virtual IEquilibrium::~IEquilibrium () = default
```




<hr>



### function ~IEquilibrium [1/3]

```C++
virtual IEquilibrium::~IEquilibrium () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryVparMu/initialisation/iequilibrium.hpp`

