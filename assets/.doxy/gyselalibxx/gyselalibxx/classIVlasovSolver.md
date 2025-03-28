

# Class IVlasovSolver



[**ClassList**](annotated.md) **>** [**IVlasovSolver**](classIVlasovSolver.md)



_An abstract class for solving a Vlasov equation._ 

* `#include <ivlasovsolver.hpp>`





Inherited by the following classes: [MpiSplitVlasovSolver](classMpiSplitVlasovSolver.md),  [SplitVlasovSolver](classSplitVlasovSolver.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVxVyXY | [**operator()**](#function-operator) (DFieldSpVxVyXY allfdistribu, DConstFieldXY efield\_x, DConstFieldXY efield\_y, double dt) const = 0<br>_Solves a Vlasov equation on a timestep dt._  |
| virtual  | [**~IVlasovSolver**](#function-ivlasovsolver) () = default<br> |




























## Public Functions Documentation




### function operator() 

_Solves a Vlasov equation on a timestep dt._ 
```C++
virtual DFieldSpVxVyXY IVlasovSolver::operator() (
    DFieldSpVxVyXY allfdistribu,
    DConstFieldXY efield_x,
    DConstFieldXY efield_y,
    double dt
) const = 0
```





**Parameters:**


* `allfdistribu` On input : the initial value of the distribution function. On output : the value of the distribution function after solving the Vlasov equation. 
* `efield_x` The electric field in the x direction computed at all spatial positions. 
* `efield_y` The electric field in the y direction computed at all spatial positions. 
* `dt` The timestep.



**Returns:**

The distribution function after solving the Vlasov equation. 





        

<hr>



### function ~IVlasovSolver 

```C++
virtual IVlasovSolver::~IVlasovSolver () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXYVxVy/vlasov/ivlasovsolver.hpp`

