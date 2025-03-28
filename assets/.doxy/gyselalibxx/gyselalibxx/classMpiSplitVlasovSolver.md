

# Class MpiSplitVlasovSolver



[**ClassList**](annotated.md) **>** [**MpiSplitVlasovSolver**](classMpiSplitVlasovSolver.md)



_A class that solves a Vlasov equation using Strang's splitting on an MPI distributed mesh._ [More...](#detailed-description)

* `#include <mpisplitvlasovsolver.hpp>`



Inherits the following classes: [IVlasovSolver](classIVlasovSolver.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MpiSplitVlasovSolver**](#function-mpisplitvlasovsolver) ([**IAdvectionSpatial**](classIAdvectionSpatial.md)&lt; [**GeometryVxVyXY**](classGeometryVxVyXY.md), [**GridX**](structGridX.md) &gt; const & advec\_x, [**IAdvectionSpatial**](classIAdvectionSpatial.md)&lt; [**GeometryVxVyXY**](classGeometryVxVyXY.md), [**GridY**](structGridY.md) &gt; const & advec\_y, [**IAdvectionVelocity**](classIAdvectionVelocity.md)&lt; [**GeometryXYVxVy**](classGeometryXYVxVy.md), [**GridVx**](structGridVx.md) &gt; const & advec\_vx, [**IAdvectionVelocity**](classIAdvectionVelocity.md)&lt; [**GeometryXYVxVy**](classGeometryXYVxVy.md), [**GridVy**](structGridVy.md) &gt; const & advec\_vy, [**MPITransposeAllToAll**](classMPITransposeAllToAll.md)&lt; [**X2DSplit**](classMPILayout.md), [**V2DSplit**](classMPILayout.md) &gt; const & transpose) <br>_Creates an instance of the split vlasov solver class._  |
| virtual DFieldSpVxVyXY | [**operator()**](#function-operator) (DFieldSpVxVyXY allfdistribu, DConstFieldXY electric\_field\_x, DConstFieldXY electric\_field\_y, double dt) override const<br>_Solves a Vlasov equation on a timestep dt._  |
|   | [**~MpiSplitVlasovSolver**](#function-mpisplitvlasovsolver) () override<br> |


## Public Functions inherited from IVlasovSolver

See [IVlasovSolver](classIVlasovSolver.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVxVyXY | [**operator()**](classIVlasovSolver.md#function-operator) (DFieldSpVxVyXY allfdistribu, DConstFieldXY efield\_x, DConstFieldXY efield\_y, double dt) const = 0<br>_Solves a Vlasov equation on a timestep dt._  |
| virtual  | [**~IVlasovSolver**](classIVlasovSolver.md#function-ivlasovsolver) () = default<br> |






















































## Detailed Description


The Vlasov equation is split between four advection equations along the [**X**](structX.md), [**Y**](structY.md), [**Vx**](structVx.md) and [**Vy**](structVy.md) directions. The splitting involves solving the advections in the [**X**](structX.md), [**Y**](structY.md), and [**Vx**](structVx.md) directions first on a time interval of length dt/2, then the Vy-direction advection on a time dt, and finally the [**X**](structX.md), [**Y**](structY.md), and [**Vx**](structVx.md) directions again in reverse order on dt/2. 


    
## Public Functions Documentation




### function MpiSplitVlasovSolver 

_Creates an instance of the split vlasov solver class._ 
```C++
MpiSplitVlasovSolver::MpiSplitVlasovSolver (
    IAdvectionSpatial < GeometryVxVyXY , GridX > const & advec_x,
    IAdvectionSpatial < GeometryVxVyXY , GridY > const & advec_y,
    IAdvectionVelocity < GeometryXYVxVy , GridVx > const & advec_vx,
    IAdvectionVelocity < GeometryXYVxVy , GridVy > const & advec_vy,
    MPITransposeAllToAll < X2DSplit , V2DSplit > const & transpose
) 
```





**Parameters:**


* `advec_x` An advection operator along the x direction. 
* `advec_y` An advection operator along the y direction. 
* `advec_vx` An advection operator along the vx direction. 
* `advec_vy` An advection operator along the vy direction. 
* `transpose` A MPI transpose operator to move between layouts. 




        

<hr>



### function operator() 

_Solves a Vlasov equation on a timestep dt._ 
```C++
virtual DFieldSpVxVyXY MpiSplitVlasovSolver::operator() (
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



### function ~MpiSplitVlasovSolver 

```C++
MpiSplitVlasovSolver::~MpiSplitVlasovSolver () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXYVxVy/vlasov/mpisplitvlasovsolver.hpp`

