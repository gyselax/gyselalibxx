

# Class BslPredCorrRTheta

**template &lt;class Mapping, class FootFinder&gt;**



[**ClassList**](annotated.md) **>** [**BslPredCorrRTheta**](classBslPredCorrRTheta.md)



_Predictor-corrector for the Vlasov-Poisson equations._ [More...](#detailed-description)

* `#include <bsl_predcorr.hpp>`



Inherits the following classes: [ITimeSolverRTheta](classITimeSolverRTheta.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslPredCorrRTheta**](#function-bslpredcorrrtheta) (Mapping const & mapping, [**BslAdvectionRTheta**](classBslAdvectionRTheta.md)&lt; FootFinder, Mapping &gt; const & advection\_solver, SplineRThetaBuilder\_host const & builder, SplineRThetaEvaluatorNullBound\_host const & rhs\_evaluator, [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), SplineRThetaEvaluatorNullBound &gt; const & poisson\_solver) <br>_Instantiate a_ [_**BslPredCorrRTheta**_](classBslPredCorrRTheta.md) _._ |
| virtual host\_t&lt; DFieldRTheta &gt; | [**operator()**](#function-operator) (host\_t&lt; DFieldRTheta &gt; allfdistribu, double const dt, int const steps) override const<br>_Solves on_  _the equations system._ |


## Public Functions inherited from ITimeSolverRTheta

See [ITimeSolverRTheta](classITimeSolverRTheta.md)

| Type | Name |
| ---: | :--- |
| virtual host\_t&lt; DFieldRTheta &gt; | [**operator()**](classITimeSolverRTheta.md#function-operator) (host\_t&lt; DFieldRTheta &gt; allfdistribu, double const dt, int const steps=1) const = 0<br>_Solves on_  _the equations system._ |
| virtual  | [**~ITimeSolverRTheta**](classITimeSolverRTheta.md#function-itimesolverrtheta) () = default<br> |
















































## Protected Functions inherited from ITimeSolverRTheta

See [ITimeSolverRTheta](classITimeSolverRTheta.md)

| Type | Name |
| ---: | :--- |
|  void | [**display\_time\_difference**](classITimeSolverRTheta.md#function-display_time_difference) (std::string const & title, std::chrono::time\_point&lt; std::chrono::system\_clock &gt; const & start\_time, std::chrono::time\_point&lt; std::chrono::system\_clock &gt; const & end\_time) const<br>_Displays the time difference between two given times and a title._  |






## Detailed Description


It solves in time the following Vlasov-Poisson equations system:



* ,
* ,
* ,




we write .


This method is mainly a Runge-Kutta 2 method:


for ,


First, it advects on a half time step:
* 1. From , it computes  with a [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md);
* 2. From , it computes  with a [**AdvectionFieldFinder**](classAdvectionFieldFinder.md);
* 3. From  and , it computes  with a [**BslAdvectionRTheta**](classBslAdvectionRTheta.md) on ;




Secondly, it advects on a full time step:
* 4. From , it computes  with a [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md);
* 5. From , it computes  with a [**AdvectionFieldFinder**](classAdvectionFieldFinder.md);
* 6. From  and , it computes  with a [**BslAdvectionRTheta**](classBslAdvectionRTheta.md) on .






**Template parameters:**


* `Mapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates. 
* `FootFinder` A IFootFinder class. 




    
## Public Functions Documentation




### function BslPredCorrRTheta 

_Instantiate a_ [_**BslPredCorrRTheta**_](classBslPredCorrRTheta.md) _._
```C++
inline BslPredCorrRTheta::BslPredCorrRTheta (
    Mapping const & mapping,
    BslAdvectionRTheta < FootFinder, Mapping > const & advection_solver,
    SplineRThetaBuilder_host const & builder,
    SplineRThetaEvaluatorNullBound_host const & rhs_evaluator,
    PolarSplineFEMPoissonLikeSolver < GridR , GridTheta , PolarBSplinesRTheta , SplineRThetaEvaluatorNullBound > const & poisson_solver
) 
```





**Parameters:**


* `mapping` The mapping function from the logical index range to the physical index range. 
* `advection_solver` The advection operator. 
* `builder` The spline builder for the computation of the RHS and the advection field. 
* `rhs_evaluator` The evaluator of B-splines for the RHS. 
* `poisson_solver` The PDE solver which computes the electrical potential. 




        

<hr>



### function operator() 

_Solves on_  _the equations system._
```C++
inline virtual host_t< DFieldRTheta > BslPredCorrRTheta::operator() (
    host_t< DFieldRTheta > allfdistribu,
    double const dt,
    int const steps
) override const
```





**Parameters:**


* `allfdistribu` On input: the initial condition. On output: the solution at . 
* `dt` The time step. 
* `steps` The number  of time interactions.



**Returns:**

A Field toward allfdistribu. 





        
Implements [*ITimeSolverRTheta::operator()*](classITimeSolverRTheta.md#function-operator)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/time_solver/bsl_predcorr.hpp`

