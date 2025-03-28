

# Class NoPerturbInitialisation



[**ClassList**](annotated.md) **>** [**NoPerturbInitialisation**](classNoPerturbInitialisation.md)



_Initialisation operator with no perturbation, i.e the distribution function equal to the Maxwellian._ 

* `#include <noperturbinitialisation.hpp>`



Inherits the following classes: [IInitialisation](classIInitialisation.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**NoPerturbInitialisation**](#function-noperturbinitialisation) (DConstFieldSpVparMu fequilibrium) <br>_Creates an instance of the_ [_**NoPerturbInitialisation**_](classNoPerturbInitialisation.md) _class._ |
| virtual DFieldSpVparMu | [**operator()**](#function-operator) (DFieldSpVparMu allfdistribu) override const<br>_Initialises the distribution function as as a Maxwellian._  |
|   | [**~NoPerturbInitialisation**](#function-noperturbinitialisation) () override<br> |


## Public Functions inherited from IInitialisation

See [IInitialisation](classIInitialisation.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVparMu | [**operator()**](classIInitialisation.md#function-operator) (DFieldSpVparMu allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual DFieldSpXVx | [**operator()**](classIInitialisation.md#function-operator_1) (DFieldSpXVx allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual DFieldSpXYVxVy | [**operator()**](classIInitialisation.md#function-operator_2) (DFieldSpXYVxVy allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |






















































## Public Functions Documentation




### function NoPerturbInitialisation 

_Creates an instance of the_ [_**NoPerturbInitialisation**_](classNoPerturbInitialisation.md) _class._
```C++
NoPerturbInitialisation::NoPerturbInitialisation (
    DConstFieldSpVparMu fequilibrium
) 
```





**Parameters:**


* `fequilibrium` A Maxwellian. 




        

<hr>



### function operator() 

_Initialises the distribution function as as a Maxwellian._ 
```C++
virtual DFieldSpVparMu NoPerturbInitialisation::operator() (
    DFieldSpVparMu allfdistribu
) override const
```





**Parameters:**


* `allfdistribu` The initialised distribution function. 



**Returns:**

The initialised distribution function. 





        
Implements [*IInitialisation::operator()*](classIInitialisation.md#function-operator)


<hr>



### function ~NoPerturbInitialisation 

```C++
NoPerturbInitialisation::~NoPerturbInitialisation () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryVparMu/initialisation/noperturbinitialisation.hpp`

