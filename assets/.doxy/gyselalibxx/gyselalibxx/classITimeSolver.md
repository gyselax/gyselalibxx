

# Class ITimeSolver



[**ClassList**](annotated.md) **>** [**ITimeSolver**](classITimeSolver.md)



_An abstract class for solving a Boltzmann-Poisson system of equations._ [More...](#detailed-description)

* `#include <itimesolver.hpp>`





Inherited by the following classes: [PredCorr](classPredCorr.md),  [PredCorr](classPredCorr.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, double time\_start, double dt, int steps=1) const = 0<br>_Operator for solving the Boltzmann-Poisson system._  |
| virtual DFieldSpVxVyXY | [**operator()**](#function-operator_1) (DFieldSpVxVyXY allfdistribu, double dt, int steps=1) const = 0<br>_Solves the Vlasov-Poisson system._  |
| virtual  | [**~ITimeSolver**](#function-itimesolver-12) () = default<br> |
| virtual  | [**~ITimeSolver**](#function-itimesolver-12) () = default<br> |




























## Detailed Description


An abstract class for solving a Vlasov-Poisson system of equations. 


    
## Public Functions Documentation




### function operator() 

_Operator for solving the Boltzmann-Poisson system._ 
```C++
virtual DFieldSpXVx ITimeSolver::operator() (
    DFieldSpXVx allfdistribu,
    double time_start,
    double dt,
    int steps=1
) const = 0
```





**Parameters:**


* `allfdistribu` On input : the initial value of the distribution function. On output : the value of the distribution function after solving the Boltzmann-Poisson system a given number of iterations. 
* `time_start` The physical time at the start of the simulation. 
* `dt` The timestep. 
* `steps` The number of iterations to be performed by the solver. 



**Returns:**

The distribution function after solving the system. 





        

<hr>



### function operator() 

_Solves the Vlasov-Poisson system._ 
```C++
virtual DFieldSpVxVyXY ITimeSolver::operator() (
    DFieldSpVxVyXY allfdistribu,
    double dt,
    int steps=1
) const = 0
```





**Parameters:**


* `allfdistribu` On input : the initial value of the distribution function. On output : the value of the distribution function after solving the Vlasov-Poisson system a given number of iterations. 
* `dt` The timestep. 
* `steps` The number of iterations to be performed by the predictor-corrector. 



**Returns:**

The distribution function after solving the system. 





        

<hr>



### function ~ITimeSolver [1/2]

```C++
virtual ITimeSolver::~ITimeSolver () = default
```




<hr>



### function ~ITimeSolver [1/2]

```C++
virtual ITimeSolver::~ITimeSolver () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/time_integration/itimesolver.hpp`

