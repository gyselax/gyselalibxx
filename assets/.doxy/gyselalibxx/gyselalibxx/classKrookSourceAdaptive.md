

# Class KrookSourceAdaptive



[**ClassList**](annotated.md) **>** [**KrookSourceAdaptive**](classKrookSourceAdaptive.md)



_A class that describes a source of particles._ [More...](#detailed-description)

* `#include <krook_source_adaptive.hpp>`



Inherits the following classes: [IRightHandSide](classIRightHandSide.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**KrookSourceAdaptive**](#function-krooksourceadaptive-12) (IdxRangeX const & gridx, IdxRangeVx const & gridvx, RhsType const type, double extent, double stiffness, double amplitude, double density, double temperature) <br>_Creates an instance of the_ [_**KrookSourceAdaptive**_](classKrookSourceAdaptive.md) _class._ |
|   | [**KrookSourceAdaptive**](#function-krooksourceadaptive-22) ([**KrookSourceAdaptive**](classKrookSourceAdaptive.md) &&) = default<br>_Creates an instance of the_ [_**KrookSourceAdaptive**_](classKrookSourceAdaptive.md) _class._ |
|  void | [**get\_amplitudes**](#function-get_amplitudes) (DFieldSpX amplitudes, DConstFieldSpXVx allfdistribu) const<br>_Computes the amplitude coefficient of the_ [_**KrookSourceAdaptive**_](classKrookSourceAdaptive.md) _operator._ |
|  void | [**get\_derivative**](#function-get_derivative) (DFieldSpXVx df, DConstFieldSpXVx f, DConstFieldSpXVx f0) const<br>_Computes the expression of the time derivative of the distribution function._  |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, double dt) override const<br>_Update the distribution function following the_ [_**KrookSourceAdaptive**_](classKrookSourceAdaptive.md) _operator._ |
|   | [**~KrookSourceAdaptive**](#function-krooksourceadaptive) () override<br> |


## Public Functions inherited from IRightHandSide

See [IRightHandSide](classIRightHandSide.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classIRightHandSide.md#function-operator) (DFieldSpXVx allfdistribu, double dt) const = 0<br>_Operator for applying the source term on the distribution function._  |
| virtual  | [**~IRightHandSide**](classIRightHandSide.md#function-irighthandside) () = default<br> |






















































## Detailed Description


The [**KrookSourceAdaptive**](classKrookSourceAdaptive.md) class solves the following evolution equation: df/dt = -amplitude \* mask \* (f - ftarget)


mask defines the spatial region where the operator is active.


ftarget is a maxwellian characterised by density and temperature, and a zero fluid velocity.


amplitude depends on space, time and the considered species so that:
* amplitude(ions) = m\_amplitude = constant 

* amplitude(electrons, x, t) = m\_amplitude (density\_ions(x,t) - m\_density) / (density\_electrons(x,t) - m\_density) 
 so that the operator conserves locally the charge.




The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/main/doc/geometryXVx/krook_source.pdf). 


    
## Public Functions Documentation




### function KrookSourceAdaptive [1/2]

_Creates an instance of the_ [_**KrookSourceAdaptive**_](classKrookSourceAdaptive.md) _class._
```C++
KrookSourceAdaptive::KrookSourceAdaptive (
    IdxRangeX const & gridx,
    IdxRangeVx const & gridvx,
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
* `gridvx` The mesh in the vx direction. 
* `type` A RhsType parameter that defines the region where the operator is active. If type = Source, the mask equals one in the central zone of the plasma of width extent; If type = Sink, the mask equals zero in the central zone of the plasma of width extent; 
* `extent` A parameter that sets the extent of the source. 
* `stiffness` A parameter that sets the stiffness of the source extent. 
* `amplitude` A parameter that sets the amplitude of the source. 
* `density` A parameter that sets the density of the Maxwellian ftarget. 
* `temperature` A parameter that sets the temperature of the Maxwellian ftarget. 




        

<hr>



### function KrookSourceAdaptive [2/2]

_Creates an instance of the_ [_**KrookSourceAdaptive**_](classKrookSourceAdaptive.md) _class._
```C++
KrookSourceAdaptive::KrookSourceAdaptive (
    KrookSourceAdaptive &&
) = default
```




<hr>



### function get\_amplitudes 

_Computes the amplitude coefficient of the_ [_**KrookSourceAdaptive**_](classKrookSourceAdaptive.md) _operator._
```C++
void KrookSourceAdaptive::get_amplitudes (
    DFieldSpX amplitudes,
    DConstFieldSpXVx allfdistribu
) const
```



This coefficient depends on the considered species and ensures that the operator conserves the charge locally.




**Parameters:**


* `amplitudes` A field that contains on output the amplitude coefficients for each species. 
* `allfdistribu` The distribution function. 




        

<hr>



### function get\_derivative 

_Computes the expression of the time derivative of the distribution function._ 
```C++
void KrookSourceAdaptive::get_derivative (
    DFieldSpXVx df,
    DConstFieldSpXVx f,
    DConstFieldSpXVx f0
) const
```



The expression is df = -amplitude \* mask \* (f - ftarget). This function is used for the time integrator ([**RK2**](classRK2.md) for instance).




**Parameters:**


* `df` The time derivative. 
* `f` The distribution function. 
* `f0` An optional parameter. 




        

<hr>



### function operator() 

_Update the distribution function following the_ [_**KrookSourceAdaptive**_](classKrookSourceAdaptive.md) _operator._
```C++
virtual DFieldSpXVx KrookSourceAdaptive::operator() (
    DFieldSpXVx allfdistribu,
    double dt
) override const
```



Update the distribution function for both electrons and ions to show how it is modified following the effect of the [**KrookSourceAdaptive**](classKrookSourceAdaptive.md) operator.




**Parameters:**


* `allfdistribu` The distribution function. 
* `dt` The time step.



**Returns:**

A field referencing the distribution function passed as argument. 





        
Implements [*IRightHandSide::operator()*](classIRightHandSide.md#function-operator)


<hr>



### function ~KrookSourceAdaptive 

```C++
KrookSourceAdaptive::~KrookSourceAdaptive () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/krook_source_adaptive.hpp`

