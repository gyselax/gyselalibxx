

# Class GMGPolarTools::HomogeneousDirichletBoundaryConditions



[**ClassList**](annotated.md) **>** [**GMGPolarTools**](namespaceGMGPolarTools.md) **>** [**HomogeneousDirichletBoundaryConditions**](classGMGPolarTools_1_1HomogeneousDirichletBoundaryConditions.md)



_Homogeneous Dirichlet boundary conditions satisfying the GMGPolar BoundaryConditions concept._ 

* `#include <gmg_polar_poisson_like_solver.hpp>`







































## Public Static Functions

| Type | Name |
| ---: | :--- |
|  KOKKOS\_INLINE\_FUNCTION double | [**u\_D**](#function-u_d) (const double & r, const double & theta) <br>_The value of the solution on the boundary._  |
|  KOKKOS\_INLINE\_FUNCTION double | [**u\_D\_Interior**](#function-u_d_interior) (const double & r, const double & theta) <br>_The value of the solution on the inner boundary (at r=rmin). Required for the concept, not needed here._  |


























## Public Static Functions Documentation




### function u\_D 

_The value of the solution on the boundary._ 
```C++
static inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::HomogeneousDirichletBoundaryConditions::u_D (
    const double & r,
    const double & theta
) 
```




<hr>



### function u\_D\_Interior 

_The value of the solution on the inner boundary (at r=rmin). Required for the concept, not needed here._ 
```C++
static inline KOKKOS_INLINE_FUNCTION double GMGPolarTools::HomogeneousDirichletBoundaryConditions::u_D_Interior (
    const double & r,
    const double & theta
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/pde_solvers/gmg_polar_poisson_like_solver.hpp`

