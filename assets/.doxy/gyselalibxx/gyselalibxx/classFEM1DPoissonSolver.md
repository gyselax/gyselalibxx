

# Class FEM1DPoissonSolver

**template &lt;class SplineBuilder, class SplineEvaluator&gt;**



[**ClassList**](annotated.md) **>** [**FEM1DPoissonSolver**](classFEM1DPoissonSolver.md)



[More...](#detailed-description)

* `#include <fem_1d_poisson_solver.hpp>`



Inherits the following classes: [IPoissonSolver](classIPoissonSolver.md)












## Classes

| Type | Name |
| ---: | :--- |
| struct | [**GridPDEDimQ**](structFEM1DPoissonSolver_1_1GridPDEDimQ.md) <br>_The grid of quadrature points along the PDEDim direction._  |
| struct | [**HiddenFEMBSplines**](structFEM1DPoissonSolver_1_1HiddenFEMBSplines.md) <br> |










































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**FEM1DPoissonSolver**](#function-fem1dpoissonsolver) (SplineBuilder const & spline\_builder, SplineEvaluator const & spline\_evaluator) <br> |
|  field\_type | [**operator()**](#function-operator) (field\_type phi, field\_type rho) override const<br>_An operator which calculates the solution_  _to Poisson's equation:_ _._ |
|  field\_type | [**operator()**](#function-operator_1) (field\_type phi, vector\_field\_type E, field\_type rho) override const<br>_An operator which calculates the solution_  _to Poisson's equation and its derivative:_ __ _._ |
|  void | [**solve\_matrix\_system**](#function-solve_matrix_system) (BatchedFEMBSplinesCoeff phi\_spline\_coef, field\_type rho) const<br> |
























































## Detailed Description


A class to solve the following equation:  using a Finite Element Method.




**Template parameters:**


* `SplineEvaluator` An evaluator which can be used to evaluate splines. 




    
## Public Functions Documentation




### function FEM1DPoissonSolver 

```C++
inline FEM1DPoissonSolver::FEM1DPoissonSolver (
    SplineBuilder const & spline_builder,
    SplineEvaluator const & spline_evaluator
) 
```



Construct the FemQNSolver operator.




**Parameters:**


* `spline_builder` A spline builder which calculates the coefficients of a spline representation. 
* `spline_evaluator` A spline evaluator which provides the value of a spline representation from its coefficients. 




        

<hr>



### function operator() 

_An operator which calculates the solution_  _to Poisson's equation:_ _._
```C++
inline field_type FEM1DPoissonSolver::operator() (
    field_type phi,
    field_type rho
) override const
```





**Parameters:**


* `phi` The solution to Poisson's equation. 
* `rho` The right-hand side of Poisson's equation.



**Returns:**

A reference to the solution to Poisson's equation. 





        

<hr>



### function operator() 

_An operator which calculates the solution_  _to Poisson's equation and its derivative:_ __ _._
```C++
inline field_type FEM1DPoissonSolver::operator() (
    field_type phi,
    vector_field_type E,
    field_type rho
) override const
```





**Parameters:**


* `phi` The solution to Poisson's equation. 
* `E` The negative derivative of the solution to Poisson's equation. 
* `rho` The right-hand side of Poisson's equation.



**Returns:**

A reference to the solution to Poisson's equation. 





        

<hr>



### function solve\_matrix\_system 

```C++
inline void FEM1DPoissonSolver::solve_matrix_system (
    BatchedFEMBSplinesCoeff phi_spline_coef,
    field_type rho
) const
```



[SHOULD BE PRIVATE (Kokkos limitation)]




**Parameters:**


* `phi_spline_coef` The spline coefficients which will describe the result on the FEM basis. 
* `rho` The function on the right hand side of the equation at the interpolation points. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/fem_1d_poisson_solver.hpp`

