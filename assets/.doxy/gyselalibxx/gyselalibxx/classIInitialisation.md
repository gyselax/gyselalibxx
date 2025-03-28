

# Class IInitialisation



[**ClassList**](annotated.md) **>** [**IInitialisation**](classIInitialisation.md)



_An abstract class that allows for initialising a distribution function._ 

* `#include <iinitialisation.hpp>`





Inherited by the following classes: [NoPerturbInitialisation](classNoPerturbInitialisation.md),  [RestartInitialisation](classRestartInitialisation.md),  [SingleModePerturbInitialisation](classSingleModePerturbInitialisation.md),  [SingleModePerturbInitialisation](classSingleModePerturbInitialisation.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVparMu | [**operator()**](#function-operator) (DFieldSpVparMu allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual DFieldSpXVx | [**operator()**](#function-operator_1) (DFieldSpXVx allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual DFieldSpXYVxVy | [**operator()**](#function-operator_2) (DFieldSpXYVxVy allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual  | [**~IInitialisation**](#function-iinitialisation-13) () = default<br> |
| virtual  | [**~IInitialisation**](#function-iinitialisation-13) () = default<br> |
| virtual  | [**~IInitialisation**](#function-iinitialisation-13) () = default<br> |




























## Public Functions Documentation




### function operator() 

_Operator for initialising a distribution function._ 
```C++
virtual DFieldSpVparMu IInitialisation::operator() (
    DFieldSpVparMu allfdistribu
) const = 0
```





**Parameters:**


* `allfdistribu` On input: the uninitialized distribution function. On output: the initialised distribution function. 



**Returns:**

The initialised distribution function. 





        

<hr>



### function operator() 

_Operator for initialising a distribution function._ 
```C++
virtual DFieldSpXVx IInitialisation::operator() (
    DFieldSpXVx allfdistribu
) const = 0
```





**Parameters:**


* `allfdistribu` On input: the uninitialized distribution function. On output: the initialised distribution function. 



**Returns:**

The initialised distribution function. 





        

<hr>



### function operator() 

_Operator for initialising a distribution function._ 
```C++
virtual DFieldSpXYVxVy IInitialisation::operator() (
    DFieldSpXYVxVy allfdistribu
) const = 0
```





**Parameters:**


* `allfdistribu` On input: the uninitialized distribution function. On output: the initialised distribution function. 



**Returns:**

The initialised distribution function. 





        

<hr>



### function ~IInitialisation [1/3]

```C++
virtual IInitialisation::~IInitialisation () = default
```




<hr>



### function ~IInitialisation [1/3]

```C++
virtual IInitialisation::~IInitialisation () = default
```




<hr>



### function ~IInitialisation [1/3]

```C++
virtual IInitialisation::~IInitialisation () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryVparMu/initialisation/iinitialisation.hpp`

