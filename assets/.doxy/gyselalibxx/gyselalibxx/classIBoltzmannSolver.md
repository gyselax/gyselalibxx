

# Class IBoltzmannSolver



[**ClassList**](annotated.md) **>** [**IBoltzmannSolver**](classIBoltzmannSolver.md)



_An abstract class for solving a Boltzmann equation._ 

* `#include <iboltzmannsolver.hpp>`





Inherited by the following classes: [SplitRightHandSideSolver](classSplitRightHandSideSolver.md),  [SplitVlasovSolver](classSplitVlasovSolver.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, DConstFieldX efield, double dt) const = 0<br>_Operator for solving the Boltzmann equation on one timestep._  |
| virtual  | [**~IBoltzmannSolver**](#function-iboltzmannsolver) () = default<br> |




























## Public Functions Documentation




### function operator() 

_Operator for solving the Boltzmann equation on one timestep._ 
```C++
virtual DFieldSpXVx IBoltzmannSolver::operator() (
    DFieldSpXVx allfdistribu,
    DConstFieldX efield,
    double dt
) const = 0
```





**Parameters:**


* `allfdistribu` On input : the initial value of the distribution function. On output : the value of the distribution function after solving the Boltzmann equation on one timestep. 
* `efield` The electric field computed at every spatial position. 
* `dt` The timestep. 



**Returns:**

The distribution function after solving the Boltzmann equation. 





        

<hr>



### function ~IBoltzmannSolver 

```C++
virtual IBoltzmannSolver::~IBoltzmannSolver () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/boltzmann/iboltzmannsolver.hpp`

