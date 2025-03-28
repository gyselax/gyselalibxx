

# Class CollisionsInter



[**ClassList**](annotated.md) **>** [**CollisionsInter**](classCollisionsInter.md)



_Class describing the inter-species collision operator._ [More...](#detailed-description)

* `#include <collisions_inter.hpp>`



Inherits the following classes: [IRightHandSide](classIRightHandSide.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CollisionsInter**](#function-collisionsinter) (IdxRangeSpXVx const & mesh, double nustar0) <br>_The constructor for the operator._  |
|  void | [**get\_derivative**](#function-get_derivative) (DFieldSpXVx df, DConstFieldSpXVx allfdistribu) const<br>_Computes the expression of the time derivative of the distribution function._  |
|  double | [**get\_nustar0**](#function-get_nustar0) () const<br>_Get the collision coefficient._  |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, double dt) override const<br>_Update the distribution function for inter-species collision._  |


## Public Functions inherited from IRightHandSide

See [IRightHandSide](classIRightHandSide.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classIRightHandSide.md#function-operator) (DFieldSpXVx allfdistribu, double dt) const = 0<br>_Operator for applying the source term on the distribution function._  |
| virtual  | [**~IRightHandSide**](classIRightHandSide.md#function-irighthandside) () = default<br> |






















































## Detailed Description


The inter-species collision operator accounts for momentum and energy transfer between the maxwellian parts of the distribution function of different species. It is solved using a explicit time integrator ([**RK2**](classRK2.md) for instance).


The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/main/doc/geometryXVx/collisions_intra_inter.pdf). 


    
## Public Functions Documentation




### function CollisionsInter 

_The constructor for the operator._ 
```C++
CollisionsInter::CollisionsInter (
    IdxRangeSpXVx const & mesh,
    double nustar0
) 
```





**Parameters:**


* `mesh` The index range on which the operator will act. 
* `nustar0` The normalised collisionality. 




        

<hr>



### function get\_derivative 

_Computes the expression of the time derivative of the distribution function._ 
```C++
void CollisionsInter::get_derivative (
    DFieldSpXVx df,
    DConstFieldSpXVx allfdistribu
) const
```



The expression is df = C, where C is the inter species collision operator. This function is made for the time integrator that is used to solve the collision operator ([**RK2**](classRK2.md) for instance).




**Parameters:**


* `df` The time derivative. 
* `allfdistribu` The distribution function. 




        

<hr>



### function get\_nustar0 

_Get the collision coefficient._ 
```C++
double CollisionsInter::get_nustar0 () const
```





**Returns:**

The collisionality. 





        

<hr>



### function operator() 

_Update the distribution function for inter-species collision._ 
```C++
virtual DFieldSpXVx CollisionsInter::operator() (
    DFieldSpXVx allfdistribu,
    double dt
) override const
```



Update the distribution function for both electrons and ions to show how it is modified following collisions between the various species. This operator only handles collisions between particles of different species.




**Parameters:**


* `allfdistribu` The distribution function. 
* `dt` The time step over which the collisions occur.



**Returns:**

A field referencing the distribution function passed as argument. 





        
Implements [*IRightHandSide::operator()*](classIRightHandSide.md#function-operator)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/collisions_inter.hpp`

