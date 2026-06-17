

# Class BslAdvectionPolar

**template &lt;class FootFinder, class LogicalToPhysicalMapping, class Builder2D, class Evaluator2D&gt;**



[**ClassList**](annotated.md) **>** [**BslAdvectionPolar**](classBslAdvectionPolar.md)



_Define an advection operator on 2D_ \((r, \theta)\) _domain._[More...](#detailed-description)

* `#include <bsl_advection_polar.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BslAdvectionPolar**](#function-bsladvectionpolar) (Builder2D const & builder\_2d, Evaluator2D const & evaluator\_2d, FootFinder & foot\_finder, LogicalToPhysicalMapping const & logical\_to\_physical\_mapping) <br>_Instantiate an advection operator._  |
|  DFieldFDistribu | [**operator()**](#function-operator) (DFieldFDistribu allfdistribu, [**DVectorConstFieldAdvection**](classVectorField.md) advection\_field, double dt) const<br>_Advect a function over a time step dt with the given advection field along the physical directions._  |
|   | [**~BslAdvectionPolar**](#function-bsladvectionpolar) () = default<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  void | [**unify\_value\_at\_centre\_pt**](#function-unify_value_at_centre_pt) (Field&lt; T, IdxRangeBatched, MemorySpace &gt; values) <br>_Replace the value at_ \((r=0, \theta)\) _point by the value at_\((r=0,0)\) _for all_\(\theta\) _._ |


























## Detailed Description


The advection operator uses a backward semi-Lagrangian method. The method is based on the property that the solution is constant along the characteristics.


For the following equation: 
\[\partial_t f(t,x) + A(t, x) \cdot \nabla_x f(t,x) = 0,\]



we write the characteristics: 
\[\partial_t X(t; s, x) = A(t, X(t; s, x)), \qquad \text{ with } X(s; s, x) = x.\]



Then the property gives us: 
\[f(t, x) = f(0, X(t; 0, x)), \quad \forall t.\]



So the first step of the advection operator is to compute the feet of the characteristics \(X(t; t+\Delta t, x_i)\) for each mesh point \(x_i\).


For the second step, we interpolate the function at the computed feet of the characteristics, and obtain the function at the next time step: \(f(t + \Delta t, x) = f(t, X(t; t+\Delta, x))\).


Different time integration methods are implemented to solve the equation of the characteristics. They are defined in the IPolarFootFinder class.


The feet can be advected on different domains (physical domain or pseudo-physical domain) which are determined in the [**PolarFootFinder**](classPolarFootFinder.md) operator.


The interpolation of the function is always done in the logical domain, where the B-splines are defined.




**See also:** IPolarFootFinder 



    
## Public Functions Documentation




### function BslAdvectionPolar 

_Instantiate an advection operator._ 
```C++
inline BslAdvectionPolar::BslAdvectionPolar (
    Builder2D const & builder_2d,
    Evaluator2D const & evaluator_2d,
    FootFinder & foot_finder,
    LogicalToPhysicalMapping const & logical_to_physical_mapping
) 
```





**Parameters:**


* `builder_2d` The 2D builder used to compute interpolations coefficients from the function values at interpolation points. 
* `evaluator_2d` The 2D evaluator used to evaluate the interpolating function at the feet of the characteristics. 
* `foot_finder` An IFootFinder which computes the feet of the characteristics. 
* `logical_to_physical_mapping` The mapping function from the logical domain to the physical domain. 




        

<hr>



### function operator() 

_Advect a function over a time step dt with the given advection field along the physical directions._ 
```C++
inline DFieldFDistribu BslAdvectionPolar::operator() (
    DFieldFDistribu allfdistribu,
    DVectorConstFieldAdvection advection_field,
    double dt
) const
```





**Parameters:**


* `allfdistribu` A Field containing the values of the function we want to advect. 
* `advection_field` A field of vectors defined on the Cartesian basis containing the values of the advection field at each point on the logical grid. 
* `dt` A time step used.



**Returns:**

A Field to allfdistribu advected on the time step given. 





        

<hr>



### function ~BslAdvectionPolar 

```C++
BslAdvectionPolar::~BslAdvectionPolar () = default
```




<hr>
## Public Static Functions Documentation




### function unify\_value\_at\_centre\_pt 

_Replace the value at_ \((r=0, \theta)\) _point by the value at_\((r=0,0)\) _for all_\(\theta\) _._
```C++
template<class T>
static inline void BslAdvectionPolar::unify_value_at_centre_pt (
    Field< T, IdxRangeBatched, MemorySpace > values
) 
```



For polar geometry, to ensure continuity at the centre point, we have to be sure that all the points for \(r = 0\) have the same value. As the computation of the values of a table can induces machine errors, this function is useful to reset the values at the central point at the same value.




**Parameters:**


* `values` The table of values we want to unify at the central point. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/bsl_advection_polar.hpp`

