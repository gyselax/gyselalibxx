

# Class SplitRightHandSideSolver



[**ClassList**](annotated.md) **>** [**SplitRightHandSideSolver**](classSplitRightHandSideSolver.md)



_A class that solves a Boltzmann equation using Strang's splitting._ [More...](#detailed-description)

* `#include <splitrighthandsidesolver.hpp>`



Inherits the following classes: [IBoltzmannSolver](classIBoltzmannSolver.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplitRightHandSideSolver**](#function-splitrighthandsidesolver) ([**IBoltzmannSolver**](classIBoltzmannSolver.md) const & vlasov\_solver, std::vector&lt; std::reference\_wrapper&lt; [**IRightHandSide**](classIRightHandSide.md) const &gt; &gt; rhs) <br>_Creates an instance of the split boltzmann solver class._  |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, DConstFieldX electric\_field, double dt) override const<br>_Solves a Boltzmann equation on a timestep dt._  |
|   | [**~SplitRightHandSideSolver**](#function-splitrighthandsidesolver) () override<br> |


## Public Functions inherited from IBoltzmannSolver

See [IBoltzmannSolver](classIBoltzmannSolver.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classIBoltzmannSolver.md#function-operator) (DFieldSpXVx allfdistribu, DConstFieldX efield, double dt) const = 0<br>_Operator for solving the Boltzmann equation on one timestep._  |
| virtual  | [**~IBoltzmannSolver**](classIBoltzmannSolver.md#function-iboltzmannsolver) () = default<br> |






















































## Detailed Description


The solver splits the Boltzmann equation and separates the advective part from the source part. The sources refers to any operator that appears on the right-hand-side of Boltzmann's equation (typically, every operator except the advections). The splitting involves solving all the source terms on a dt/2 timestep, then solving the advections on a dt timestep using a Vlasov solver, then solving the sources again on dt/2 in reverse order. 


    
## Public Functions Documentation




### function SplitRightHandSideSolver 

_Creates an instance of the split boltzmann solver class._ 
```C++
SplitRightHandSideSolver::SplitRightHandSideSolver (
    IBoltzmannSolver const & vlasov_solver,
    std::vector< std::reference_wrapper< IRightHandSide const > > rhs
) 
```





**Parameters:**


* `vlasov_solver` A solver for the associated Vlasov equation (the boltzmann equation with no sources). 
* `rhs` A vector containing all of the source terms of the considered Boltzmann equation. 




        

<hr>



### function operator() 

_Solves a Boltzmann equation on a timestep dt._ 
```C++
virtual DFieldSpXVx SplitRightHandSideSolver::operator() (
    DFieldSpXVx allfdistribu,
    DConstFieldX electric_field,
    double dt
) override const
```





**Parameters:**


* `allfdistribu` On input: the initial value of the distribution function. On output: the value of the distribution function after solving the Boltzmann equation. 
* `electric_field` The electric field computed at all spatial positions. 
* `dt` The timestep. 



**Returns:**

The distribution function after solving the Boltzmann equation. 





        
Implements [*IBoltzmannSolver::operator()*](classIBoltzmannSolver.md#function-operator)


<hr>



### function ~SplitRightHandSideSolver 

```C++
SplitRightHandSideSolver::~SplitRightHandSideSolver () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/boltzmann/splitrighthandsidesolver.hpp`

