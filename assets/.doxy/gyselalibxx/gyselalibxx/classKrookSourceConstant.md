

# Class KrookSourceConstant



[**ClassList**](annotated.md) **>** [**KrookSourceConstant**](classKrookSourceConstant.md)



_A class that describes a source of particles._ [More...](#detailed-description)

* `#include <krook_source_constant.hpp>`



Inherits the following classes: [IRightHandSide](classIRightHandSide.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**KrookSourceConstant**](#function-krooksourceconstant-12) (IdxRangeX const & gridx, IdxRangeVx const & gridv, RhsType const type, double extent, double stiffness, double amplitude, double density, double temperature) <br>_Creates an instance of the_ [_**KrookSourceConstant**_](classKrookSourceConstant.md) _class._ |
|   | [**KrookSourceConstant**](#function-krooksourceconstant-22) ([**KrookSourceConstant**](classKrookSourceConstant.md) &&) = default<br>_Creates an instance of the_ [_**KrookSourceConstant**_](classKrookSourceConstant.md) _class._ |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, double dt) override const<br>_Update the distribution function following the_ [_**KrookSourceConstant**_](classKrookSourceConstant.md) _operator._ |
|   | [**~KrookSourceConstant**](#function-krooksourceconstant) () override<br> |


## Public Functions inherited from IRightHandSide

See [IRightHandSide](classIRightHandSide.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classIRightHandSide.md#function-operator) (DFieldSpXVx allfdistribu, double dt) const = 0<br>_Operator for applying the source term on the distribution function._  |
| virtual  | [**~IRightHandSide**](classIRightHandSide.md#function-irighthandside) () = default<br> |






















































## Detailed Description


The [**KrookSourceConstant**](classKrookSourceConstant.md) class solves the following evolution equation: df/dt = -amplitude \* mask \* (f - ftarget)


mask defines the spatial region where the operator is active.


ftarget is a maxwellian characterised by density and temperature, and a zero fluid velocity.


amplitude is a constant for both species.


The solution of the evolution equation is therefore : f(t+dt) = ftarget + (f(t)-ftarget)\*exp(-amplitude\*mask\*dt) 


    
## Public Functions Documentation




### function KrookSourceConstant [1/2]

_Creates an instance of the_ [_**KrookSourceConstant**_](classKrookSourceConstant.md) _class._
```C++
KrookSourceConstant::KrookSourceConstant (
    IdxRangeX const & gridx,
    IdxRangeVx const & gridv,
    RhsType const type,
    double extent,
    double stiffness,
    double amplitude,
    double density,
    double temperature
) 
```





**Parameters:**


* `gridx` The mesh in the x direction. 
* `gridv` The mesh in the vx direction. 
* `type` A RhsType parameter that defines the region where the operator is active. If type = Source, the mask equals one in the central zone of the plasma of width extent; If type = Sink, the mask equals zero in the central zone of the plasma of width extent; 
* `extent` A parameter that sets the extent of the source. 
* `stiffness` A parameter that sets the stiffness of the source extent. 
* `amplitude` A parameter that sets the the amplitude of the source. 
* `density` A parameter that sets the density of the Maxwellian ftarget. 
* `temperature` A parameter that sets the temperature of the Maxwellian ftarget. 




        

<hr>



### function KrookSourceConstant [2/2]

_Creates an instance of the_ [_**KrookSourceConstant**_](classKrookSourceConstant.md) _class._
```C++
KrookSourceConstant::KrookSourceConstant (
    KrookSourceConstant &&
) = default
```




<hr>



### function operator() 

_Update the distribution function following the_ [_**KrookSourceConstant**_](classKrookSourceConstant.md) _operator._
```C++
virtual DFieldSpXVx KrookSourceConstant::operator() (
    DFieldSpXVx allfdistribu,
    double dt
) override const
```



Update the distribution function for both electrons and ions to show how it is modified following the effect of the [**KrookSourceConstant**](classKrookSourceConstant.md) operator.




**Parameters:**


* `allfdistribu` The distribution function. 
* `dt` The time step.



**Returns:**

A field referencing the distribution function passed as argument. 





        
Implements [*IRightHandSide::operator()*](classIRightHandSide.md#function-operator)


<hr>



### function ~KrookSourceConstant 

```C++
KrookSourceConstant::~KrookSourceConstant () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/krook_source_constant.hpp`

