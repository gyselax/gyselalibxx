

# Class KineticSource



[**ClassList**](annotated.md) **>** [**KineticSource**](classKineticSource.md)



_A class that describes a source of particles._ [More...](#detailed-description)

* `#include <kinetic_source.hpp>`



Inherits the following classes: [IRightHandSide](classIRightHandSide.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**KineticSource**](#function-kineticsource) (IdxRangeX const & gridx, IdxRangeVx const & gridv, double extent, double stiffness, double amplitude, double density, double energy, double temperature) <br>_Creates an instance of the_ [_**KineticSource**_](classKineticSource.md) _class._ |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, double dt) override const<br>_Update the distribution function following the_ [_**KineticSource**_](classKineticSource.md) _operator._ |
|   | [**~KineticSource**](#function-kineticsource) () override<br> |


## Public Functions inherited from IRightHandSide

See [IRightHandSide](classIRightHandSide.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classIRightHandSide.md#function-operator) (DFieldSpXVx allfdistribu, double dt) const = 0<br>_Operator for applying the source term on the distribution function._  |
| virtual  | [**~IRightHandSide**](classIRightHandSide.md#function-irighthandside) () = default<br> |






















































## Detailed Description


The [**KineticSource**](classKineticSource.md) class solves the following evolution equation: df/dt = S\_kin Where S\_kin = spatial\_extent(x) \* velocity\_shape(v) Since S\_kin does not depend on time, we have f(t+dt) = f(t) + S\_kin\*dt as a solution of this evolution equation.


spatial\_extent defines the location where the source is active. spatial\_extent is normalised, so that its integral along the spatial direction equals one. It has a hyperbolic tangent shape. It is equal to one in a central zone of the plasma of width defined by the extent parameter.


velocity\_shape defines the velocity profile of the source in the parallel velocity direction. It is the sum of a source that injects only density, and a source that injects only energy. If the density and energy parameters are equal to one (usual case), the resulting velocity\_shape is maxwellian.


The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/devel/doc/geometryXVx/kinetic_source.pdf). 


    
## Public Functions Documentation




### function KineticSource 

_Creates an instance of the_ [_**KineticSource**_](classKineticSource.md) _class._
```C++
KineticSource::KineticSource (
    IdxRangeX const & gridx,
    IdxRangeVx const & gridv,
    double extent,
    double stiffness,
    double amplitude,
    double density,
    double energy,
    double temperature
) 
```





**Parameters:**


* `gridx` The mesh in the x direction. 
* `gridv` The mesh in the vx direction. 
* `extent` A parameter that sets the spatial extent of the source. 
* `stiffness` A parameter that sets stiffness of the source extent. 
* `amplitude` A parameter that sets the amplitude of the source. 
* `density` A parameter that sets the density of the source. 
* `energy` A parameter that sets the energy of the source. 
* `temperature` A parameter that sets the temperature of the source. 




        

<hr>



### function operator() 

_Update the distribution function following the_ [_**KineticSource**_](classKineticSource.md) _operator._
```C++
virtual DFieldSpXVx KineticSource::operator() (
    DFieldSpXVx allfdistribu,
    double dt
) override const
```



Update the distribution function for both electrons and ions to show how it is modified following the effect of the [**KineticSource**](classKineticSource.md) operator.




**Parameters:**


* `allfdistribu` The distribution function. 
* `dt` The time step over which the collisions occur.



**Returns:**

A field referencing the distribution function passed as argument. 





        
Implements [*IRightHandSide::operator()*](classIRightHandSide.md#function-operator)


<hr>



### function ~KineticSource 

```C++
KineticSource::~KineticSource () override
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/kinetic_source.hpp`

