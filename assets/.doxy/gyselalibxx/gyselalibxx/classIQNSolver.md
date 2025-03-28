

# Class IQNSolver



[**ClassList**](annotated.md) **>** [**IQNSolver**](classIQNSolver.md)



_Base class for a Quasi-Neutrality solver._ [More...](#detailed-description)

* `#include <iqnsolver.hpp>`





Inherited by the following classes: [NullQNSolver](classNullQNSolver.md),  [NullQNSolver](classNullQNSolver.md),  [QNSolver](classQNSolver.md),  [QNSolver](classQNSolver.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual void | [**operator()**](#function-operator) (host\_t&lt; DFieldRTheta &gt; electrostatic\_potential, host\_t&lt; [**DVectorFieldRTheta**](classVectorField.md)&lt; [**X**](structX.md), [**Y**](structY.md) &gt; &gt; electric\_field, host\_t&lt; DConstFieldRTheta &gt; allfdistribu) const = 0<br>_Compute the electrical potential and the electric field from the Quasi-Neutrality equation._  |
| virtual void | [**operator()**](#function-operator_1) (DFieldX electrostatic\_potential, DFieldX electric\_field, DConstFieldSpXVx allfdistribu) const = 0<br> |
| virtual void | [**operator()**](#function-operator_2) (DFieldXY electrostatic\_potential, DFieldXY electric\_field\_x, DFieldXY electric\_field\_y, DConstFieldSpVxVyXY allfdistribu) const = 0<br> |
| virtual  | [**~IQNSolver**](#function-iqnsolver-13) () = default<br> |
| virtual  | [**~IQNSolver**](#function-iqnsolver-13) () = default<br> |
| virtual  | [**~IQNSolver**](#function-iqnsolver-13) () = default<br> |




























## Detailed Description


An operator which solves the Quasi-Neutrality equation using a fast Fourier transform.


An operator which solves the Quasi-Neutrality equation.


An operator which solves the Quasi-Neutrality equation: 


An operator which solves the Quasi-Neutrality equation:  using a fast Fourier transform on a periodic index range. This operator only works for equidistant points. 


    
## Public Functions Documentation




### function operator() 

_Compute the electrical potential and the electric field from the Quasi-Neutrality equation._ 
```C++
virtual void IQNSolver::operator() (
    host_t< DFieldRTheta > electrostatic_potential,
    host_t< DVectorFieldRTheta < X , Y > > electric_field,
    host_t< DConstFieldRTheta > allfdistribu
) const = 0
```





**Parameters:**


* `electrostatic_potential` The solution of the Quasi-Neutrality equation. 
* `electric_field` The electric field . 
* `allfdistribu` The rhs of the Quasi-Neutrality equation. 




        

<hr>



### function operator() 

```C++
virtual void IQNSolver::operator() (
    DFieldX electrostatic_potential,
    DFieldX electric_field,
    DConstFieldSpXVx allfdistribu
) const = 0
```



The operator which solves the equation using the method described by the class.




**Parameters:**


* `electrostatic_potential` The electrostatic potential, the result of the poisson solver. 
* `electric_field` The electric field, the derivative of the electrostatic potential. 
* `allfdistribu` The distribution function. 




        

<hr>



### function operator() 

```C++
virtual void IQNSolver::operator() (
    DFieldXY electrostatic_potential,
    DFieldXY electric_field_x,
    DFieldXY electric_field_y,
    DConstFieldSpVxVyXY allfdistribu
) const = 0
```



The operator which solves the equation using the method described by the class.




**Parameters:**


* `electrostatic_potential` The electrostatic potential, the result of the poisson solver. 
* `electric_field_x` The x-component of the electric field, the gradient of the electrostatic potential. 
* `electric_field_y` The y-component of the electric field, the gradient of the electrostatic potential. 
* `allfdistribu` The distribution function. 




        

<hr>



### function ~IQNSolver [1/3]

```C++
virtual IQNSolver::~IQNSolver () = default
```




<hr>



### function ~IQNSolver [1/3]

```C++
virtual IQNSolver::~IQNSolver () = default
```




<hr>



### function ~IQNSolver [1/3]

```C++
virtual IQNSolver::~IQNSolver () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/poisson/iqnsolver.hpp`

