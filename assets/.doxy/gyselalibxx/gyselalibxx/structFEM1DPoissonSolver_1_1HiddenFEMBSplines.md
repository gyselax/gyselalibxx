

# Struct FEM1DPoissonSolver::HiddenFEMBSplines



[**ClassList**](annotated.md) **>** [**FEM1DPoissonSolver**](classFEM1DPoissonSolver.md) **>** [**HiddenFEMBSplines**](structFEM1DPoissonSolver_1_1HiddenFEMBSplines.md)



[More...](#detailed-description)

* `#include <fem_1d_poisson_solver.hpp>`



Inherits the following classes: ddc::NonUniformBSplines< PDEDim, InputBSplines::degree()>






























































## Detailed Description


The internal type of the FEM basis. This is used if the provided BSplines are uniform as non-uniform BSplines are required to correctly enforce Dirichlet boundary conditions.


This type should be private but is public due to Kokkos restrictions. 


    

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/fem_1d_poisson_solver.hpp`

