

# Class SplitVlasovSolver



[**ClassList**](annotated.md) **>** [**SplitVlasovSolver**](classSplitVlasovSolver.md)



_A class that solves a Vlasov equation using Strang's splitting._ [More...](#detailed-description)

* `#include <splitvlasovsolver.hpp>`



Inherits the following classes: [IBoltzmannSolver](classIBoltzmannSolver.md),  [IVlasovSolver](classIVlasovSolver.md)










































































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplitVlasovSolver**](#function-splitvlasovsolver-12) ([**IAdvectionSpatial**](classIAdvectionSpatial.md)&lt; [**GeometryXVx**](classGeometryXVx.md), [**GridX**](structGridX.md) &gt; const & advec\_x, [**IAdvectionVelocity**](classIAdvectionVelocity.md)&lt; [**GeometryXVx**](classGeometryXVx.md), [**GridVx**](structGridVx.md) &gt; const & advec\_vx) <br>_Creates an instance of the split vlasov solver class._  |
|   | [**SplitVlasovSolver**](#function-splitvlasovsolver-22) ([**IAdvectionSpatial**](classIAdvectionSpatial.md)&lt; [**GeometryVxVyXY**](classGeometryVxVyXY.md), [**GridX**](structGridX.md) &gt; const & advec\_x, [**IAdvectionSpatial**](classIAdvectionSpatial.md)&lt; [**GeometryVxVyXY**](classGeometryVxVyXY.md), [**GridY**](structGridY.md) &gt; const & advec\_y, [**IAdvectionVelocity**](classIAdvectionVelocity.md)&lt; [**GeometryVxVyXY**](classGeometryVxVyXY.md), [**GridVx**](structGridVx.md) &gt; const & advec\_vx, [**IAdvectionVelocity**](classIAdvectionVelocity.md)&lt; [**GeometryVxVyXY**](classGeometryVxVyXY.md), [**GridVy**](structGridVy.md) &gt; const & advec\_vy) <br>_Creates an instance of the split vlasov solver class._  |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, DConstFieldX electric\_field, double dt) override const<br>_Solves a Vlasov equation on a timestep dt._  |
| virtual DFieldSpVxVyXY | [**operator()**](#function-operator_1) (DFieldSpVxVyXY allfdistribu, DConstFieldXY electric\_field\_x, DConstFieldXY electric\_field\_y, double dt) override const<br>_Solves a Vlasov equation on a timestep dt._  |
|   | [**~SplitVlasovSolver**](#function-splitvlasovsolver-12) () override<br> |
|   | [**~SplitVlasovSolver**](#function-splitvlasovsolver-12) () override<br> |


## Public Functions inherited from IBoltzmannSolver

See [IBoltzmannSolver](classIBoltzmannSolver.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classIBoltzmannSolver.md#function-operator) (DFieldSpXVx allfdistribu, DConstFieldX efield, double dt) const = 0<br>_Operator for solving the Boltzmann equation on one timestep._  |
| virtual  | [**~IBoltzmannSolver**](classIBoltzmannSolver.md#function-iboltzmannsolver) () = default<br> |


## Public Functions inherited from IVlasovSolver

See [IVlasovSolver](classIVlasovSolver.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVxVyXY | [**operator()**](classIVlasovSolver.md#function-operator) (DFieldSpVxVyXY allfdistribu, DConstFieldXY efield\_x, DConstFieldXY efield\_y, double dt) const = 0<br>_Solves a Vlasov equation on a timestep dt._  |
| virtual  | [**~IVlasovSolver**](classIVlasovSolver.md#function-ivlasovsolver) () = default<br> |
















































































## Detailed Description


The Vlasov equation is split between two advection equations along the [**X**](structX.md) and [**Vx**](structVx.md) directions. The splitting involves solving the x-direction advection first on a time interval of length dt/2, then the vx-direction advection on a time dt, and then x-direction again on dt/2.


The Vlasov equation is split between four advection equations along the [**X**](structX.md), [**Y**](structY.md), [**Vx**](structVx.md) and [**Vy**](structVy.md) directions. The splitting involves solving the advections in the [**X**](structX.md), [**Y**](structY.md), and [**Vx**](structVx.md) directions first on a time interval of length dt/2, then the Vy-direction advection on a time dt, and finally the [**X**](structX.md), [**Y**](structY.md), and [**Vx**](structVx.md) directions again in reverse order on dt/2. 


    
## Public Functions Documentation




### function SplitVlasovSolver [1/2]

_Creates an instance of the split vlasov solver class._ 
```C++
SplitVlasovSolver::SplitVlasovSolver (
    IAdvectionSpatial < GeometryXVx , GridX > const & advec_x,
    IAdvectionVelocity < GeometryXVx , GridVx > const & advec_vx
) 
```





**Parameters:**


* `advec_x` An advection operator along the x direction. 
* `advec_vx` An advection operator along the vx direction. 




        

<hr>



### function SplitVlasovSolver [2/2]

_Creates an instance of the split vlasov solver class._ 
```C++
SplitVlasovSolver::SplitVlasovSolver (
    IAdvectionSpatial < GeometryVxVyXY , GridX > const & advec_x,
    IAdvectionSpatial < GeometryVxVyXY , GridY > const & advec_y,
    IAdvectionVelocity < GeometryVxVyXY , GridVx > const & advec_vx,
    IAdvectionVelocity < GeometryVxVyXY , GridVy > const & advec_vy
) 
```





**Parameters:**


* `advec_x` An advection operator along the x direction. 
* `advec_y` An advection operator along the y direction. 
* `advec_vx` An advection operator along the vx direction. 
* `advec_vy` An advection operator along the vy direction. 




        

<hr>



### function operator() 

_Solves a Vlasov equation on a timestep dt._ 
```C++
virtual DFieldSpXVx SplitVlasovSolver::operator() (
    DFieldSpXVx allfdistribu,
    DConstFieldX electric_field,
    double dt
) override const
```





**Parameters:**


* `allfdistribu` On input : the initial value of the distribution function. On output : the value of the distribution function after solving the Vlasov equation. 
* `electric_field` The electric field computed at all spatial positions. 
* `dt` The timestep. 



**Returns:**

The distribution function after solving the Vlasov equation. 





        
Implements [*IBoltzmannSolver::operator()*](classIBoltzmannSolver.md#function-operator)


<hr>



### function operator() 

_Solves a Vlasov equation on a timestep dt._ 
```C++
virtual DFieldSpVxVyXY SplitVlasovSolver::operator() (
    DFieldSpVxVyXY allfdistribu,
    DConstFieldXY electric_field_x,
    DConstFieldXY electric_field_y,
    double dt
) override const
```





**Parameters:**


* `allfdistribu` On input : the initial value of the distribution function. On output : the value of the distribution function after solving the Vlasov equation. 
* `electric_field_x` The electric field in the x direction computed at all spatial positions. 
* `electric_field_y` The electric field in the y direction computed at all spatial positions. 
* `dt` The timestep.



**Returns:**

The distribution function after solving the Vlasov equation. 





        
Implements [*IVlasovSolver::operator()*](classIVlasovSolver.md#function-operator)


<hr>



### function ~SplitVlasovSolver [1/2]

```C++
SplitVlasovSolver::~SplitVlasovSolver () override
```




<hr>



### function ~SplitVlasovSolver [1/2]

```C++
SplitVlasovSolver::~SplitVlasovSolver () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/boltzmann/splitvlasovsolver.hpp`

