

# Class ChargeDensityCalculator



[**ClassList**](annotated.md) **>** [**ChargeDensityCalculator**](classChargeDensityCalculator.md)



_A class which computes charges density with Kokkos._ [More...](#detailed-description)

* `#include <chargedensitycalculator.hpp>`



Inherits the following classes: [IChargeDensityCalculator](classIChargeDensityCalculator.md),  [IChargeDensityCalculator](classIChargeDensityCalculator.md)










































































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ChargeDensityCalculator**](#function-chargedensitycalculator-12) (DConstFieldVx coeffs) <br>_Create a_ [_**ChargeDensityCalculator**_](classChargeDensityCalculator.md) _object._ |
|   | [**ChargeDensityCalculator**](#function-chargedensitycalculator-22) (DConstFieldVxVy coeffs) <br>_Create a_ [_**ChargeDensityCalculator**_](classChargeDensityCalculator.md) _object._ |
| virtual DFieldX | [**operator()**](#function-operator) (DFieldX rho, DConstFieldSpXVx allfdistribu) const<br>_Computes the charge density rho from the distribution function._  |
| virtual void | [**operator()**](#function-operator_1) (DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const<br>_Computes the charge density rho from the distribution function._  |


## Public Functions inherited from IChargeDensityCalculator

See [IChargeDensityCalculator](classIChargeDensityCalculator.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldX | [**operator()**](classIChargeDensityCalculator.md#function-operator) (DFieldX rho, DConstFieldSpXVx allfdistribu) const = 0<br> |
| virtual void | [**operator()**](classIChargeDensityCalculator.md#function-operator_1) (DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const = 0<br> |


## Public Functions inherited from IChargeDensityCalculator

See [IChargeDensityCalculator](classIChargeDensityCalculator.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldX | [**operator()**](classIChargeDensityCalculator.md#function-operator) (DFieldX rho, DConstFieldSpXVx allfdistribu) const = 0<br> |
| virtual void | [**operator()**](classIChargeDensityCalculator.md#function-operator_1) (DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const = 0<br> |
















































































## Detailed Description


A class which computes charges density by solving the equation:  where  is the charge of the species  and  is the distribution function.


A class which computes charges density by solving the equation:  where  is the charge of the species  and  is the distribution function. 


    
## Public Functions Documentation




### function ChargeDensityCalculator [1/2]

_Create a_ [_**ChargeDensityCalculator**_](classChargeDensityCalculator.md) _object._
```C++
explicit ChargeDensityCalculator::ChargeDensityCalculator (
    DConstFieldVx coeffs
) 
```





**Parameters:**


* `coeffs` The coefficients of the quadrature. 




        

<hr>



### function ChargeDensityCalculator [2/2]

_Create a_ [_**ChargeDensityCalculator**_](classChargeDensityCalculator.md) _object._
```C++
explicit ChargeDensityCalculator::ChargeDensityCalculator (
    DConstFieldVxVy coeffs
) 
```





**Parameters:**


* `coeffs` The coefficients of the quadrature. 




        

<hr>



### function operator() 

_Computes the charge density rho from the distribution function._ 
```C++
virtual DFieldX ChargeDensityCalculator::operator() (
    DFieldX rho,
    DConstFieldSpXVx allfdistribu
) const
```





**Parameters:**


* `rho` 
* `allfdistribu` 



**Returns:**

rho The charge density. 





        
Implements [*IChargeDensityCalculator::operator()*](classIChargeDensityCalculator.md#function-operator)


<hr>



### function operator() 

_Computes the charge density rho from the distribution function._ 
```C++
virtual void ChargeDensityCalculator::operator() (
    DFieldXY rho,
    DConstFieldSpVxVyXY allfdistribu
) const
```





**Parameters:**


* `rho` 
* `allfdistribu` 




        
Implements [*IChargeDensityCalculator::operator()*](classIChargeDensityCalculator.md#function-operator_1)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/poisson/chargedensitycalculator.hpp`

