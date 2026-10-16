

# Class PolarSplineFEMPoissonLikeSolver::CoeffEvaluator

**template &lt;class Evaluator, class Coeff&gt;**



[**ClassList**](annotated.md) **>** [**PolarSplineFEMPoissonLikeSolver**](classPolarSplineFEMPoissonLikeSolver.md) **>** [**CoeffEvaluator**](classPolarSplineFEMPoissonLikeSolver_1_1CoeffEvaluator.md)



_A wrapper that binds an evaluator with its batched coefficient field to present a single callable_ `double operator()(IdxBatch, CoordRTheta)` _._[More...](#detailed-description)

* `#include <polar_spline_fem_poisson_like_solver.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CoeffEvaluator**](#function-coeffevaluator) (Evaluator const & evaluator, Coeff coeff) <br>_Constructor of an evaluator wrapper._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**operator()**](#function-operator) (IdxBatch idx\_batch, CoordRTheta const & coord) const<br>_Evaluate the interpolation at the specified batch and coordinate._  |




























## Detailed Description


This allows the \(\alpha\) and \(\beta\) coefficients to be passed to `PolarSplineFEMPoissonLikeAssembler`, which expects a generic, batch-index-aware callable.




**Template parameters:**


* `Evaluator` The type of the 2D evaluator. 
* `Coeff` The type of the batched spline coefficient field. 




    
## Public Functions Documentation




### function CoeffEvaluator 

_Constructor of an evaluator wrapper._ 
```C++
inline PolarSplineFEMPoissonLikeSolver::CoeffEvaluator::CoeffEvaluator (
    Evaluator const & evaluator,
    Coeff coeff
) 
```




<hr>



### function operator() 

_Evaluate the interpolation at the specified batch and coordinate._ 
```C++
inline KOKKOS_INLINE_FUNCTION double PolarSplineFEMPoissonLikeSolver::CoeffEvaluator::operator() (
    IdxBatch idx_batch,
    CoordRTheta const & coord
) const
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/polar_spline_fem_poisson_like_solver.hpp`

