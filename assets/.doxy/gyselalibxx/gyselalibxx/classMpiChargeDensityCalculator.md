

# Class MpiChargeDensityCalculator



[**ClassList**](annotated.md) **>** [**MpiChargeDensityCalculator**](classMpiChargeDensityCalculator.md)



_A class which computes charges density with Kokkos._ [More...](#detailed-description)

* `#include <mpichargedensitycalculator.hpp>`



Inherits the following classes: [IChargeDensityCalculator](classIChargeDensityCalculator.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MpiChargeDensityCalculator**](#function-mpichargedensitycalculator) (MPI\_Comm comm, [**IChargeDensityCalculator**](classIChargeDensityCalculator.md) const & local\_charge\_density\_calculator) <br>_Create a_ [_**MpiChargeDensityCalculator**_](classMpiChargeDensityCalculator.md) _object._ |
| virtual void | [**operator()**](#function-operator) (DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const<br>_Computes the charge density rho from the distribution function._  |


## Public Functions inherited from IChargeDensityCalculator

See [IChargeDensityCalculator](classIChargeDensityCalculator.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldX | [**operator()**](classIChargeDensityCalculator.md#function-operator) (DFieldX rho, DConstFieldSpXVx allfdistribu) const = 0<br> |
| virtual void | [**operator()**](classIChargeDensityCalculator.md#function-operator_1) (DFieldXY rho, DConstFieldSpVxVyXY allfdistribu) const = 0<br> |






















































## Detailed Description


A class which computes charges density by solving the equation:  where  is the charge of the species  and  is the distribution function. 


    
## Public Functions Documentation




### function MpiChargeDensityCalculator 

_Create a_ [_**MpiChargeDensityCalculator**_](classMpiChargeDensityCalculator.md) _object._
```C++
explicit MpiChargeDensityCalculator::MpiChargeDensityCalculator (
    MPI_Comm comm,
    IChargeDensityCalculator const & local_charge_density_calculator
) 
```





**Parameters:**


* `comm` The MPI communicator across which the calculation is carried out. 
* `local_charge_density_calculator` An operator which calculates the density locally on a given MPI node. The results from this operator will then be combined using MPI. 




        

<hr>



### function operator() 

_Computes the charge density rho from the distribution function._ 
```C++
virtual void MpiChargeDensityCalculator::operator() (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXYVxVy/poisson/mpichargedensitycalculator.hpp`

