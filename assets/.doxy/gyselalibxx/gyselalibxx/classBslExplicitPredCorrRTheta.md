

# Class BslExplicitPredCorrRTheta

**template &lt;class LogicalToPhysicalMapping, class LogicalToPseudoPhysicalMapping&gt;**



[**ClassList**](annotated.md) **>** [**BslExplicitPredCorrRTheta**](classBslExplicitPredCorrRTheta.md)



_A second order explicit predictor-corrector for the Vlasov-Poisson equations._ [More...](#detailed-description)

* `#include <bsl_predcorr_second_order_explicit.hpp>`



Inherits the following classes: [ITimeSolverRTheta](classITimeSolverRTheta.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslExplicitPredCorrRTheta**](#function-bslexplicitpredcorrrtheta) (LogicalToPhysicalMapping const & logical\_to\_physical, LogicalToPseudoPhysicalMapping const & logical\_to\_pseudo\_physical, [**BslAdvectionRTheta**](classBslAdvectionPolar.md) const & advection\_solver, IdxRangeRTheta const & grid, SplineRThetaBuilder\_host const & builder, [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md)&lt; [**GridR**](structGridR.md), [**GridTheta**](structGridTheta.md), [**PolarBSplinesRTheta**](structPolarBSplinesRTheta.md), SplineRThetaEvaluatorNullBound &gt; const & poisson\_solver, SplineRThetaEvaluatorConstBound\_host const & advection\_evaluator) <br>_Instantiate a_ [_**BslExplicitPredCorrRTheta**_](classBslExplicitPredCorrRTheta.md) _._ |
| virtual host\_t&lt; DFieldRTheta &gt; | [**operator()**](#function-operator) (host\_t&lt; DFieldRTheta &gt; allfdistribu, double const dt, int const steps) const<br>_Solves on_  _the equations system._ |


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


The second order explicit predictor-corrector is also detailed in Edoardo Zoni's article ([https://doi.org/10.1016/j.jcp.2019.108889](https://doi.org/10.1016/j.jcp.2019.108889)).


for ,


First, it predicts:
* 1. From , it computes  with a [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md);
* 2. From , it computes  with a [**AdvectionFieldFinder**](classAdvectionFieldFinder.md);
* 3. From  and , it computes  with a [**BslAdvectionPolar**](classBslAdvectionPolar.md) on ;




We write  the characteristic feet such that .


Secondly, it corrects:
* 4. From , it computes  with a [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md);
* 5. From , it computes  with a [**AdvectionFieldFinder**](classAdvectionFieldFinder.md);
* 6. From  and , it computes  with a [**BslAdvectionPolar**](classBslAdvectionPolar.md) on .




(With  the characteristic feet such that .)




**Template parameters:**


* `LogicalToPhysicalMapping` A class describing a mapping from curvilinear coordinates to Cartesian coordinates. 
* `LogicalToPseudoPhysicalMapping` A class describing a mapping from curvilinear coordinates to pseudo-Cartesian coordinates. 




    
## Public Functions Documentation




### function BslExplicitPredCorrRTheta 

_Instantiate a_ [_**BslExplicitPredCorrRTheta**_](classBslExplicitPredCorrRTheta.md) _._
```C++
inline BslExplicitPredCorrRTheta::BslExplicitPredCorrRTheta (
    LogicalToPhysicalMapping const & logical_to_physical,
    LogicalToPseudoPhysicalMapping const & logical_to_pseudo_physical,
    BslAdvectionRTheta const & advection_solver,
    IdxRangeRTheta const & grid,
    SplineRThetaBuilder_host const & builder,
    PolarSplineFEMPoissonLikeSolver < GridR , GridTheta , PolarBSplinesRTheta , SplineRThetaEvaluatorNullBound > const & poisson_solver,
    SplineRThetaEvaluatorConstBound_host const & advection_evaluator
) 
```





**Parameters:**


* `logical_to_physical` The mapping from the logical domain to the physical domain. 
* `logical_to_pseudo_physical` The mapping from the logical domain to the pseudo-physical domain. 
* `advection_solver` The advection operator with an [**Euler**](classEuler.md) method. 
* `grid` The index range on which the functions are defined. 
* `builder` A spline builder to get the spline representation of the advection field and the RHS. 
* `poisson_solver` The PDE solver which computes the electrical potential. 
* `advection_evaluator` An evaluator of B-splines for the spline advection field. 




        

<hr>



### function operator() 

_Solves on_  _the equations system._
```C++
inline virtual host_t< DFieldRTheta > BslExplicitPredCorrRTheta::operator() (
    host_t< DFieldRTheta > allfdistribu,
    double const dt,
    int const steps
) const
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/time_solver/bsl_predcorr_second_order_explicit.hpp`

