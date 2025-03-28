

# Class PredCorrRK2XY

**template &lt;class PoissonSolver, class AdvectionX, class AdvectionY&gt;**



[**ClassList**](annotated.md) **>** [**PredCorrRK2XY**](classPredCorrRK2XY.md)



_Predictor-corrector based on_ [_**RK2**_](classRK2.md) _for the guiding-centre model._[More...](#detailed-description)

* `#include <predcorr_RK2.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**PredCorrRK2XY**](#function-predcorrrk2xy) (PoissonSolver const & poisson\_solver, AdvectionX const & advection\_x, AdvectionY const & advection\_y) <br>_Instantiate the predictor-corrector._  |
|  void | [**operator()**](#function-operator) (DFieldXY allfdistribu, double const dt, int const nbiter) <br>_Apply the predictor-corrector method on several time steps._  |
|   | [**~PredCorrRK2XY**](#function-predcorrrk2xy) () = default<br> |




























## Detailed Description


It solves in time the following guiding-centre equations system:



* ,
* ,
* .




This method is mainly a Runge-Kutta 2 method:


for ,


First, it advects on a half time step:
* 1./2. From , it computes  with a [**FFTPoissonSolver**](classFFTPoissonSolver.md);
* 3. From  and , it computes  with a [**BslAdvection1D**](classBslAdvection1D.md) on ;




Secondly, it advects on a full time step:
* 4./5. From , it computes  with a [**FFTPoissonSolver**](classFFTPoissonSolver.md);
* 6. From  and , it computes  with a BslAdvectionRP on .






**Template parameters:**


* `PoissonSolver` Type of the Poisson solver applied in the method. 
* `AdvectionX` Type of the 1D advection operator applied to advect along [**X**](structX.md). 
* `AdvectionY` Type of the 1D advection operator applied to advect along [**Y**](structY.md). 




    
## Public Functions Documentation




### function PredCorrRK2XY 

_Instantiate the predictor-corrector._ 
```C++
inline PredCorrRK2XY::PredCorrRK2XY (
    PoissonSolver const & poisson_solver,
    AdvectionX const & advection_x,
    AdvectionY const & advection_y
) 
```





**Parameters:**


* `poisson_solver` Poisson solver also computing the electric field. 
 
* `advection_x` 1D advection operator along  direction. 
* `advection_y` 1D advection operator along  direction. 




        

<hr>



### function operator() 

_Apply the predictor-corrector method on several time steps._ 
```C++
inline void PredCorrRK2XY::operator() (
    DFieldXY allfdistribu,
    double const dt,
    int const nbiter
) 
```



Along the simulation, the data are saved in an output folder.




**Parameters:**


* `allfdistribu` Initial function . 
* `dt` Time step. 
* `nbiter` Number of time steps. 




        

<hr>



### function ~PredCorrRK2XY 

```C++
PredCorrRK2XY::~PredCorrRK2XY () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXY/time_integration/predcorr_RK2.hpp`

