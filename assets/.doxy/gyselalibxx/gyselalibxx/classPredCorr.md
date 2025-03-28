

# Class PredCorr



[**ClassList**](annotated.md) **>** [**PredCorr**](classPredCorr.md)



_A class that solves a Boltzmann-Poisson system of equations using a predictor-corrector scheme._ [More...](#detailed-description)

* `#include <predcorr.hpp>`



Inherits the following classes: [ITimeSolver](classITimeSolver.md),  [ITimeSolver](classITimeSolver.md)










































































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PredCorr**](#function-predcorr-12) ([**IBoltzmannSolver**](classIBoltzmannSolver.md) const & boltzmann\_solver, [**IQNSolver**](classIQNSolver.md) const & poisson\_solver) <br>_Creates an instance of the predictor-corrector class._  |
|   | [**PredCorr**](#function-predcorr-22) ([**IVlasovSolver**](classIVlasovSolver.md) const & vlasov\_solver, [**IQNSolver**](classIQNSolver.md) const & poisson\_solver) <br>_Creates an instance of the predictor-corrector class._  |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, double time\_start, double dt, int steps=1) override const<br>_Solves the Boltzmann-Poisson system._  |
| virtual DFieldSpVxVyXY | [**operator()**](#function-operator_1) (DFieldSpVxVyXY allfdistribu, double dt, int steps=1) override const<br>_Solves the Vlasov-Poisson system._  |
|   | [**~PredCorr**](#function-predcorr-12) () override<br> |
|   | [**~PredCorr**](#function-predcorr-12) () override<br> |


## Public Functions inherited from ITimeSolver

See [ITimeSolver](classITimeSolver.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classITimeSolver.md#function-operator) (DFieldSpXVx allfdistribu, double time\_start, double dt, int steps=1) const = 0<br>_Operator for solving the Boltzmann-Poisson system._  |
| virtual DFieldSpVxVyXY | [**operator()**](classITimeSolver.md#function-operator_1) (DFieldSpVxVyXY allfdistribu, double dt, int steps=1) const = 0<br>_Solves the Vlasov-Poisson system._  |
| virtual  | [**~ITimeSolver**](classITimeSolver.md#function-itimesolver-12) () = default<br> |
| virtual  | [**~ITimeSolver**](classITimeSolver.md#function-itimesolver-12) () = default<br> |


## Public Functions inherited from ITimeSolver

See [ITimeSolver](classITimeSolver.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classITimeSolver.md#function-operator) (DFieldSpXVx allfdistribu, double time\_start, double dt, int steps=1) const = 0<br>_Operator for solving the Boltzmann-Poisson system._  |
| virtual DFieldSpVxVyXY | [**operator()**](classITimeSolver.md#function-operator_1) (DFieldSpVxVyXY allfdistribu, double dt, int steps=1) const = 0<br>_Solves the Vlasov-Poisson system._  |
| virtual  | [**~ITimeSolver**](classITimeSolver.md#function-itimesolver-12) () = default<br> |
| virtual  | [**~ITimeSolver**](classITimeSolver.md#function-itimesolver-12) () = default<br> |
















































































## Detailed Description


A class that solves a Vlasov-Poisson system of equations using a predictor-corrector scheme.


A class that solves a Boltzmann-Poisson system with a predictor corrector scheme. This scheme consists in estimating the electric potential after a time interval of a half-timestep. This potential is then used to compute the value of the distribution function at time t+dt, where dt is the timestep.


A class that solves a Vlasov-Poisson system with a predictor corrector scheme. This scheme consists in estimating the electric potential after a time interval of a half-timestep. This potential is then used to compute the value of the distribution function at time t+dt, where dt is the timestep. 


    
## Public Functions Documentation




### function PredCorr [1/2]

_Creates an instance of the predictor-corrector class._ 
```C++
PredCorr::PredCorr (
    IBoltzmannSolver const & boltzmann_solver,
    IQNSolver const & poisson_solver
) 
```





**Parameters:**


* `boltzmann_solver` A solver for a Boltzmann equation. 
* `poisson_solver` A solver for a Quasi-Neutrality equation. 




        

<hr>



### function PredCorr [2/2]

_Creates an instance of the predictor-corrector class._ 
```C++
PredCorr::PredCorr (
    IVlasovSolver const & vlasov_solver,
    IQNSolver const & poisson_solver
) 
```





**Parameters:**


* `vlasov_solver` A solver for a Boltzmann equation. 
* `poisson_solver` A solver for a Poisson equation. 




        

<hr>



### function operator() 

_Solves the Boltzmann-Poisson system._ 
```C++
virtual DFieldSpXVx PredCorr::operator() (
    DFieldSpXVx allfdistribu,
    double time_start,
    double dt,
    int steps=1
) override const
```





**Parameters:**


* `allfdistribu` On input : the initial value of the distribution function. On output : the value of the distribution function after solving the Boltzmann-Poisson system a given number of iterations. 
* `time_start` The physical time at the start of the simulation. 
* `dt` The timestep. 
* `steps` The number of iterations to be performed by the predictor-corrector. 



**Returns:**

The distribution function after solving the system. 





        
Implements [*ITimeSolver::operator()*](classITimeSolver.md#function-operator)


<hr>



### function operator() 

_Solves the Vlasov-Poisson system._ 
```C++
virtual DFieldSpVxVyXY PredCorr::operator() (
    DFieldSpVxVyXY allfdistribu,
    double dt,
    int steps=1
) override const
```





**Parameters:**


* `allfdistribu` On input : the initial value of the distribution function. On output : the value of the distribution function after solving the Vlasov-Poisson system a given number of iterations. 
* `dt` The timestep. 
* `steps` The number of iterations to be performed by the predictor-corrector. 



**Returns:**

The distribution function after solving the system. 





        
Implements [*ITimeSolver::operator()*](classITimeSolver.md#function-operator_1)


<hr>



### function ~PredCorr [1/2]

```C++
PredCorr::~PredCorr () override
```




<hr>



### function ~PredCorr [1/2]

```C++
PredCorr::~PredCorr () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/time_integration/predcorr.hpp`

