

# Class IChargeDensityCalculator



[**ClassList**](annotated.md) **>** [**IChargeDensityCalculator**](classIChargeDensityCalculator.md)



_A class which computes charges density._ [More...](#detailed-description)

* `#include <ichargedensitycalculator.hpp>`





Inherited by the following classes: [ChargeDensityCalculator](classChargeDensityCalculator.md),  [ChargeDensityCalculator](classChargeDensityCalculator.md),  [MpiChargeDensityCalculator](classMpiChargeDensityCalculator.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual DFieldX | [**operator()**](#function-operator) (DFieldX rho, DConstFieldSpXVx allfdistribu) const = 0<br> |
| virtual void | [**operator()**](#function-operator_1) (DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const = 0<br> |




























## Detailed Description


A class which calculates the charge density.


A class which calculates the charge density. This is then used as the right hand side of the Quasi-Neutrality equation.


A class which computes charges density by solving the equation:  where  is the charge of the species  and  is the distribution function. 


    
## Public Functions Documentation




### function operator() 

```C++
virtual DFieldX IChargeDensityCalculator::operator() (
    DFieldX rho,
    DConstFieldSpXVx allfdistribu
) const = 0
```



Calculate the charge density rho from the distribution function.




**Parameters:**


* `rho` The charge density. 
* `allfdistribu` The distribution function.



**Returns:**

rho The charge density. 





        

<hr>



### function operator() 

```C++
virtual void IChargeDensityCalculator::operator() (
    DFieldXY rho,
    DConstFieldSpVxVyXY allfdistribu
) const = 0
```



Calculate the charge density rho from the distribution function.


Calculate the charge density by calculating the spline representation of slices of the distribution function at each spatial point along the velocity direction. This representation is then integrated and multiplied by the charge to find the charge density.




**Parameters:**


* `rho` The charge density. 
* `allfdistribu` The distribution function. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/poisson/ichargedensitycalculator.hpp`

