

# Class CollisionsIntra



[**ClassList**](annotated.md) **>** [**CollisionsIntra**](classCollisionsIntra.md)



_Class describing the intra-species collision operator._ [More...](#detailed-description)

* `#include <collisions_intra.hpp>`



Inherits the following classes: [IRightHandSide](classIRightHandSide.md)












## Classes

| Type | Name |
| ---: | :--- |
| struct | [**GhostedVx**](structCollisionsIntra_1_1GhostedVx.md) <br> |
| struct | [**GhostedVxStaggered**](structCollisionsIntra_1_1GhostedVxStaggered.md) <br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GhostedVx**](structCollisionsIntra_1_1GhostedVx.md) &gt; | [**IdxRangeSpXVx\_ghosted**](#typedef-idxrangespxvx_ghosted)  <br> |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GhostedVxStaggered**](structCollisionsIntra_1_1GhostedVxStaggered.md) &gt; | [**IdxRangeSpXVx\_ghosted\_staggered**](#typedef-idxrangespxvx_ghosted_staggered)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GhostedVx**](structCollisionsIntra_1_1GhostedVx.md) &gt; | [**IdxSpXVx\_ghosted**](#typedef-idxspxvx_ghosted)  <br> |
| typedef Idx&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GhostedVxStaggered**](structCollisionsIntra_1_1GhostedVxStaggered.md) &gt; | [**IdxSpXVx\_ghosted\_staggered**](#typedef-idxspxvx_ghosted_staggered)  <br> |
| typedef Idx&lt; [**GhostedVx**](structCollisionsIntra_1_1GhostedVx.md) &gt; | [**IdxVx\_ghosted**](#typedef-idxvx_ghosted)  <br> |
| typedef Idx&lt; [**GhostedVxStaggered**](structCollisionsIntra_1_1GhostedVxStaggered.md) &gt; | [**IdxVx\_ghosted\_staggered**](#typedef-idxvx_ghosted_staggered)  <br> |








































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CollisionsIntra**](#function-collisionsintra) (IdxRangeSpXVx const & mesh, double nustar0) <br>_The constructor for the operator._  |
|  void | [**compute\_matrix\_coeff**](#function-compute_matrix_coeff) (DFieldSpXVx AA, DFieldSpXVx BB, DFieldSpXVx CC, Field&lt; double, [**IdxRangeSpXVx\_ghosted**](classCollisionsIntra.md#typedef-idxrangespxvx_ghosted) &gt; Dcoll, Field&lt; double, [**IdxRangeSpXVx\_ghosted\_staggered**](classCollisionsIntra.md#typedef-idxrangespxvx_ghosted_staggered) &gt; Dcoll\_staggered, Field&lt; double, [**IdxRangeSpXVx\_ghosted**](classCollisionsIntra.md#typedef-idxrangespxvx_ghosted) &gt; Nucoll, double deltat) const<br>_Compute the coefficients of the tridiagonal matrix of the collision operator linear system._  |
|  void | [**compute\_rhs\_vector**](#function-compute_rhs_vector) (DFieldSpXVx RR, DConstFieldSpXVx AA, DConstFieldSpXVx BB, DConstFieldSpXVx CC, DConstFieldSpXVx allfdistribu, double fthresh) const<br>_Compute the right-hand-side of the collision operator linear system._  |
|  void | [**fill\_matrix\_with\_coeff**](#function-fill_matrix_with_coeff) ([**Matrix\_Banded**](classMatrix__Banded.md) & matrix, host\_t&lt; DConstFieldVx &gt; AA, host\_t&lt; DConstFieldVx &gt; BB, host\_t&lt; DConstFieldVx &gt; CC) const<br>_Fill the tridiagonal matrix of the collision operator linear system with the required coefficients._  |
|  IdxRange&lt; [**GhostedVx**](structCollisionsIntra_1_1GhostedVx.md) &gt; const & | [**get\_gridvx\_ghosted**](#function-get_gridvx_ghosted) () const<br>_Get the ghosted vx mesh used for computing finite differences centred derivatives._  |
|  IdxRange&lt; [**GhostedVxStaggered**](structCollisionsIntra_1_1GhostedVxStaggered.md) &gt; const & | [**get\_gridvx\_ghosted\_staggered**](#function-get_gridvx_ghosted_staggered) () const<br>_Get the ghosted and staggered vx mesh used for computing finite differences centred derivatives._  |
|  IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), [**GhostedVx**](structCollisionsIntra_1_1GhostedVx.md) &gt; const & | [**get\_mesh\_ghosted**](#function-get_mesh_ghosted) () const<br>_Get a mesh containing the species, spatial and the ghosted vx mesh._  |
|  double | [**get\_nustar0**](#function-get_nustar0) () const<br>_Get the collision coefficient._  |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu, double dt) override const<br>_Update the distribution function for intra-species collision._  |


## Public Functions inherited from IRightHandSide

See [IRightHandSide](classIRightHandSide.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpXVx | [**operator()**](classIRightHandSide.md#function-operator) (DFieldSpXVx allfdistribu, double dt) const = 0<br>_Operator for applying the source term on the distribution function._  |
| virtual  | [**~IRightHandSide**](classIRightHandSide.md#function-irighthandside) () = default<br> |






















































## Detailed Description


The intra-species collision operator can be written as  Where Dcoll and Nucoll are the collision operator diffusion and advection coefficients that depend on the absolute value of the velocity and possibly on space. f\_a is the distribution function. The Nucoll and Dcoll coefficients are adjusted at each time and spatial position so that the collision operator conserves the density, momentum and energy of the considered species it acts on.


The evolution equation for the collisions  is solved using a Crank-Nicolson finite difference scheme adapted for non-uniform meshes. The derivatives are computed using centred second-order finite differences. To use the same derivatives formula at the edges of the vx mesh, we introduce a ghosted vx mesh with additional points at the edges. To further improve the accuracy of the computations of the derivatives, we also introduce another vx mesh that is staggered with respect to the initial mesh. The points of the staggered mesh lie at the middle of the cells of the initial vx mesh (where a cell is defined by two adjacent points). The Crank-Nicolson scheme leads to the formulation of a linear system that needs to be resolved at each spatial position of the simulation box. Note that this linear system depends on the considered spatial position.


The complete description of the operator can be found in [rhs docs](https://github.com/gyselax/gyselalibxx/blob/devel/doc/geometryXVx/collisions_intra_inter.pdf). 


    
## Public Types Documentation




### typedef IdxRangeSpXVx\_ghosted 

```C++
using CollisionsIntra::IdxRangeSpXVx_ghosted =  IdxRange<Species, GridX, GhostedVx>;
```



A type representing a mesh for species, space and ghosted vx mesh. 


        

<hr>



### typedef IdxRangeSpXVx\_ghosted\_staggered 

```C++
using CollisionsIntra::IdxRangeSpXVx_ghosted_staggered =  IdxRange<Species, GridX, GhostedVxStaggered>;
```



A type representing a mesh for species, space and ghosted staggered vx mesh. 


        

<hr>



### typedef IdxSpXVx\_ghosted 

```C++
using CollisionsIntra::IdxSpXVx_ghosted =  Idx<Species, GridX, GhostedVx>;
```



A type representing a species, space and ghosted vx index. 


        

<hr>



### typedef IdxSpXVx\_ghosted\_staggered 

```C++
using CollisionsIntra::IdxSpXVx_ghosted_staggered =  Idx<Species, GridX, GhostedVxStaggered>;
```



A type representing a species, space and ghosted staggered vx index. 


        

<hr>



### typedef IdxVx\_ghosted 

```C++
using CollisionsIntra::IdxVx_ghosted =  Idx<GhostedVx>;
```



A type representing a ghosted vx index. 


        

<hr>



### typedef IdxVx\_ghosted\_staggered 

```C++
using CollisionsIntra::IdxVx_ghosted_staggered =  Idx<GhostedVxStaggered>;
```



A type representing a ghosted staggered vx index. 


        

<hr>
## Public Functions Documentation




### function CollisionsIntra 

_The constructor for the operator._ 
```C++
CollisionsIntra::CollisionsIntra (
    IdxRangeSpXVx const & mesh,
    double nustar0
) 
```





**Parameters:**


* `mesh` The index range on which the operator will act. 
* `nustar0` The normalised collisionality. 




        

<hr>



### function compute\_matrix\_coeff 

_Compute the coefficients of the tridiagonal matrix of the collision operator linear system._ 
```C++
void CollisionsIntra::compute_matrix_coeff (
    DFieldSpXVx AA,
    DFieldSpXVx BB,
    DFieldSpXVx CC,
    Field< double, IdxRangeSpXVx_ghosted > Dcoll,
    Field< double, IdxRangeSpXVx_ghosted_staggered > Dcoll_staggered,
    Field< double, IdxRangeSpXVx_ghosted > Nucoll,
    double deltat
) const
```





**Parameters:**


* `AA` A vector representing the lower diagonal of the matrix of the linear system. 
* `BB` A vector representing the diagonal of the matrix of the linear system. 
* `CC` A vector representing the upper diagonal of the matrix of the linear system. 
* `Dcoll` The velocity-dependent diffusion coefficient of the collision operator. 
* `Dcoll_staggered` The velocity-dependent diffusion coefficient of the collision operator computed on the staggered vx mesh. 
* `Nucoll` The velocity-dependent advection coefficient of the collision operator. 
* `deltat` The time step. 




        

<hr>



### function compute\_rhs\_vector 

_Compute the right-hand-side of the collision operator linear system._ 
```C++
void CollisionsIntra::compute_rhs_vector (
    DFieldSpXVx RR,
    DConstFieldSpXVx AA,
    DConstFieldSpXVx BB,
    DConstFieldSpXVx CC,
    DConstFieldSpXVx allfdistribu,
    double fthresh
) const
```





**Parameters:**


* `RR` A vector representing the right-hand-side of the system. 
* `AA` A vector representing the lower diagonal of the matrix of the linear system. 
* `BB` A vector representing the diagonal of the matrix of the linear system. 
* `CC` A vector representing the upper diagonal of the matrix of the linear system. 
* `allfdistribu` The distribution function. 
* `fthresh` A constant value used for imposing Dirichlet boundary condition to solve the linear system. 
 




        

<hr>



### function fill\_matrix\_with\_coeff 

_Fill the tridiagonal matrix of the collision operator linear system with the required coefficients._ 
```C++
void CollisionsIntra::fill_matrix_with_coeff (
    Matrix_Banded & matrix,
    host_t< DConstFieldVx > AA,
    host_t< DConstFieldVx > BB,
    host_t< DConstFieldVx > CC
) const
```





**Parameters:**


* `matrix` A banded (tridiagonal) matrix. 
* `AA` A vector containing the lower diagonal coefficients of the matrix. 
* `BB` A vector containing the diagonal coefficients of the matrix. 
* `CC` A vector containing the upper diagonal coefficients of the matrix. 




        

<hr>



### function get\_gridvx\_ghosted 

_Get the ghosted vx mesh used for computing finite differences centred derivatives._ 
```C++
IdxRange< GhostedVx > const & CollisionsIntra::get_gridvx_ghosted () const
```





**Returns:**

The ghosted vx mesh. 





        

<hr>



### function get\_gridvx\_ghosted\_staggered 

_Get the ghosted and staggered vx mesh used for computing finite differences centred derivatives._ 
```C++
IdxRange< GhostedVxStaggered > const & CollisionsIntra::get_gridvx_ghosted_staggered () const
```





**Returns:**

The ghosted and staggered vx mesh. 





        

<hr>



### function get\_mesh\_ghosted 

_Get a mesh containing the species, spatial and the ghosted vx mesh._ 
```C++
IdxRange< Species , GridX , GhostedVx > const & CollisionsIntra::get_mesh_ghosted () const
```





**Returns:**

The species, spatial, and ghosted vx mesh. 





        

<hr>



### function get\_nustar0 

_Get the collision coefficient._ 
```C++
double CollisionsIntra::get_nustar0 () const
```





**Returns:**

The collisionality coefficient. 





        

<hr>



### function operator() 

_Update the distribution function for intra-species collision._ 
```C++
virtual DFieldSpXVx CollisionsIntra::operator() (
    DFieldSpXVx allfdistribu,
    double dt
) override const
```



Update the distribution function for both electrons and ions to show how it is modified following collisions within various species. This operator only handles collisions between particles of the same species.




**Parameters:**


* `allfdistribu` The distribution function. 
* `dt` The time step over which the collisions occur.



**Returns:**

A field referencing the distribution function passed as argument. 





        
Implements [*IRightHandSide::operator()*](classIRightHandSide.md#function-operator)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/collisions_intra.hpp`

