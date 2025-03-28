

# Class IRightHandSide



[**ClassList**](annotated.md) **>** [**IRightHandSide**](classIRightHandSide.md)



_An abstract class representing a source in Boltzmann equation._ 

* `#include <irighthandside.hpp>`





Inherited by the following classes: [CollisionsInter](classCollisionsInter.md),  [CollisionsIntra](classCollisionsIntra.md),  [KineticSource](classKineticSource.md),  [KrookSourceAdaptive](classKrookSourceAdaptive.md),  [KrookSourceConstant](classKrookSourceConstant.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, double dt) const = 0<br>_Operator for applying the source term on the distribution function._  |
| virtual  | [**~IRightHandSide**](#function-irighthandside) () = default<br> |




























## Public Functions Documentation




### function operator() 

_Operator for applying the source term on the distribution function._ 
```C++
virtual DFieldSpXVx IRightHandSide::operator() (
    DFieldSpXVx allfdistribu,
    double dt
) const = 0
```



The source  acts on the distribution function following the evolution equation df/dt = S 

**Parameters:**


* `allfdistribu` On input: the initial value of the distribution function. On output: the value of the distribution function after solving the source evolution equation on one timestep. 
* `dt` The timestep. 



**Returns:**

The distribution function after solving the source evolution equation. 





        

<hr>



### function ~IRightHandSide 

```C++
virtual IRightHandSide::~IRightHandSide () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/irighthandside.hpp`

