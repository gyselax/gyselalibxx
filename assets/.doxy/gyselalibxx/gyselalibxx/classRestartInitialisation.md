

# Class RestartInitialisation



[**ClassList**](annotated.md) **>** [**RestartInitialisation**](classRestartInitialisation.md)



_A class that initialises the distribution function from a previous simulation._ [More...](#detailed-description)

* `#include <restartinitialisation.hpp>`



Inherits the following classes: [IInitialisation](classIInitialisation.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**RestartInitialisation**](#function-restartinitialisation) (int iter\_start, double & time\_start) <br>_Create an initialisation object._  |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu) override const<br>_Triggers a PDI event to fill the distribution function with values from a hdf5 file._  |
|   | [**~RestartInitialisation**](#function-restartinitialisation) () override<br> |


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






















































## Detailed Description


A class that triggers a PDI event to read the values of a distribution function saved in a hdf5 file. These values are copied to the field that represents the distribution function. 


    
## Public Functions Documentation




### function RestartInitialisation 

_Create an initialisation object._ 
```C++
RestartInitialisation::RestartInitialisation (
    int iter_start,
    double & time_start
) 
```





**Parameters:**


* `iter_start` An integer representing the number of iteration already performed to produce the distribution function used to initialise the current simulation. 
* `time_start` The physical time corresponding to iter\_start. 




        

<hr>



### function operator() 

_Triggers a PDI event to fill the distribution function with values from a hdf5 file._ 
```C++
virtual DFieldSpXVx RestartInitialisation::operator() (
    DFieldSpXVx allfdistribu
) override const
```





**Parameters:**


* `allfdistribu` The distribution function initialised with the values read from an external file. 



**Returns:**

The initialised distribution function. 





        
Implements [*IInitialisation::operator()*](classIInitialisation.md#function-operator_1)


<hr>



### function ~RestartInitialisation 

```C++
RestartInitialisation::~RestartInitialisation () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/initialisation/restartinitialisation.hpp`

