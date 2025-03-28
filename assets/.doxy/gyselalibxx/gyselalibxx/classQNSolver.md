

# Class QNSolver



[**ClassList**](annotated.md) **>** [**QNSolver**](classQNSolver.md)



_An operator which solves the Quasi-Neutrality equation using a fast Fourier transform._ [More...](#detailed-description)

* `#include <qnsolver.hpp>`



Inherits the following classes: [IQNSolver](classIQNSolver.md),  [IQNSolver](classIQNSolver.md)










































































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**QNSolver**](#function-qnsolver-12) ([**PoissonSolver**](classIPoissonSolver.md) const & solve\_poisson, [**IChargeDensityCalculator**](classIChargeDensityCalculator.md) const & compute\_rho) <br> |
|   | [**QNSolver**](#function-qnsolver-12) ([**PoissonSolver**](classIPoissonSolver.md) const & solve\_poisson, [**IChargeDensityCalculator**](classIChargeDensityCalculator.md) const & compute\_rho) <br> |
| virtual void | [**operator()**](#function-operator) (DFieldX electrostatic\_potential, DFieldX electric\_field, DConstFieldSpXVx allfdistribu) override const<br> |
| virtual void | [**operator()**](#function-operator_1) (DFieldXY electrostatic\_potential, DFieldXY electric\_field\_x, DFieldXY electric\_field\_y, DConstFieldSpVxVyXY allfdistribu) override const<br> |
|   | [**~QNSolver**](#function-qnsolver-12) () override<br> |
|   | [**~QNSolver**](#function-qnsolver-12) () override<br> |


## Public Functions inherited from IQNSolver

See [IQNSolver](classIQNSolver.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIQNSolver.md#function-operator) (host\_t&lt; DFieldRTheta &gt; electrostatic\_potential, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; &gt; electric\_field, host\_t&lt; DConstFieldRTheta &gt; allfdistribu) const = 0<br>_Compute the electrical potential and the electric field from the Quasi-Neutrality equation._  |
| virtual void | [**operator()**](classIQNSolver.md#function-operator_1) (DFieldX electrostatic\_potential, DFieldX electric\_field, DConstFieldSpXVx allfdistribu) const = 0<br> |
| virtual void | [**operator()**](classIQNSolver.md#function-operator_2) (DFieldXY electrostatic\_potential, DFieldXY electric\_field\_x, DFieldXY electric\_field\_y, DConstFieldSpVxVyXY allfdistribu) const = 0<br> |
| virtual  | [**~IQNSolver**](classIQNSolver.md#function-iqnsolver-13) () = default<br> |
| virtual  | [**~IQNSolver**](classIQNSolver.md#function-iqnsolver-13) () = default<br> |
| virtual  | [**~IQNSolver**](classIQNSolver.md#function-iqnsolver-13) () = default<br> |


## Public Functions inherited from IQNSolver

See [IQNSolver](classIQNSolver.md)

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](classIQNSolver.md#function-operator) (host\_t&lt; DFieldRTheta &gt; electrostatic\_potential, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; &gt; electric\_field, host\_t&lt; DConstFieldRTheta &gt; allfdistribu) const = 0<br>_Compute the electrical potential and the electric field from the Quasi-Neutrality equation._  |
| virtual void | [**operator()**](classIQNSolver.md#function-operator_1) (DFieldX electrostatic\_potential, DFieldX electric\_field, DConstFieldSpXVx allfdistribu) const = 0<br> |
| virtual void | [**operator()**](classIQNSolver.md#function-operator_2) (DFieldXY electrostatic\_potential, DFieldXY electric\_field\_x, DFieldXY electric\_field\_y, DConstFieldSpVxVyXY allfdistribu) const = 0<br> |
| virtual  | [**~IQNSolver**](classIQNSolver.md#function-iqnsolver-13) () = default<br> |
| virtual  | [**~IQNSolver**](classIQNSolver.md#function-iqnsolver-13) () = default<br> |
| virtual  | [**~IQNSolver**](classIQNSolver.md#function-iqnsolver-13) () = default<br> |
















































































## Detailed Description


An operator which solves the Quasi-Neutrality equation:  using a fast Fourier transform on a periodic index range. This operator only works for equidistant points.


The electric field,  is calculated using a spline interpolation implemented in ElectricField. 


    
## Public Functions Documentation




### function QNSolver [1/2]

```C++
QNSolver::QNSolver (
    PoissonSolver const & solve_poisson,
    IChargeDensityCalculator const & compute_rho
) 
```



Construct the FftQNSolver operator.




**Parameters:**


* `solve_poisson` The operator which solves the Poisson solver. 
* `compute_rho` The operator which calculates the charge density, the right hand side of the equation. 




        

<hr>



### function QNSolver [1/2]

```C++
QNSolver::QNSolver (
    PoissonSolver const & solve_poisson,
    IChargeDensityCalculator const & compute_rho
) 
```



Construct the [**QNSolver**](classQNSolver.md) operator.




**Parameters:**


* `solve_poisson` The operator which solves the Poisson solver. 
* `compute_rho` The operator which calculates the charge density, the right hand side of the equation. 




        

<hr>



### function operator() 

```C++
virtual void QNSolver::operator() (
    DFieldX electrostatic_potential,
    DFieldX electric_field,
    DConstFieldSpXVx allfdistribu
) override const
```



The operator which solves the equation using the method described by the class.




**Parameters:**


* `electrostatic_potential` The electrostatic potential, the result of the poisson solver. 
* `electric_field` The electric field, the derivative of the electrostatic potential. 
* `allfdistribu` The distribution function. 




        
Implements [*IQNSolver::operator()*](classIQNSolver.md#function-operator_1)


<hr>



### function operator() 

```C++
virtual void QNSolver::operator() (
    DFieldXY electrostatic_potential,
    DFieldXY electric_field_x,
    DFieldXY electric_field_y,
    DConstFieldSpVxVyXY allfdistribu
) override const
```



The operator which solves the equation using the method described by the class.




**Parameters:**


* `electrostatic_potential` The electrostatic potential, the result of the poisson solver. 
* `electric_field_x` The x-component of the electric field, the gradient of the electrostatic potential. 
* `electric_field_y` The y-component of the electric field, the gradient of the electrostatic potential. 
* `allfdistribu` The distribution function. 




        
Implements [*IQNSolver::operator()*](classIQNSolver.md#function-operator_2)


<hr>



### function ~QNSolver [1/2]

```C++
QNSolver::~QNSolver () override
```




<hr>



### function ~QNSolver [1/2]

```C++
QNSolver::~QNSolver () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/poisson/qnsolver.hpp`

