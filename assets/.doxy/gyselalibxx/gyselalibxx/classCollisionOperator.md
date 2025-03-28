

# Class CollisionOperator

**template &lt;class [**CollisionConfiguration**](classCollisionConfiguration.md)&gt;**



[**ClassList**](annotated.md) **>** [**CollisionOperator**](classCollisionOperator.md)



_A class which computes the collision operator in (Sp,vpar,mu)._ 

* `#include <collision_operator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**CollisionConfiguration**](classCollisionConfiguration.md) | [**CollisionConfigurationType**](#typedef-collisionconfigurationtype)  <br>_The type of collision configuration. Namely, what kind of geometry the operator will ... operate on._  |
| typedef typename [**CollisionConfigurationType::IdxRangeDistributionFunctionType**](classCollisionConfiguration.md#typedef-idxrangedistributionfunctiontype) | [**IdxRangeDistributionFunctionType**](#typedef-idxrangedistributionfunctiontype)  <br>_The distribution function that we should expect from the operator user/caller._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CollisionOperator**](#function-collisionoperator) ([**CollisionConfigurationType**](classCollisionOperator.md#typedef-collisionconfigurationtype) const & collision\_configuration) <br>_Construct a new Collision Operator object._  |
|  void | [**operator()**](#function-operator) (DField&lt; [**IdxRangeDistributionFunctionType**](classCollisionOperator.md#typedef-idxrangedistributionfunctiontype) &gt; all\_f\_distribution, double deltat\_coll) const<br>_Apply the collision operator to the distribution functions of all species on all species._  |
|   | [**~CollisionOperator**](#function-collisionoperator) () <br>_Destroy the Collision operator object._  |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  ::koliop\_Operator | [**m\_operator\_handle**](#variable-m_operator_handle)  <br>_Opaque C type representing the operator._  |




















## Public Types Documentation




### typedef CollisionConfigurationType 

_The type of collision configuration. Namely, what kind of geometry the operator will ... operate on._ 
```C++
using CollisionOperator< CollisionConfiguration >::CollisionConfigurationType =  CollisionConfiguration;
```




<hr>



### typedef IdxRangeDistributionFunctionType 

_The distribution function that we should expect from the operator user/caller._ 
```C++
using CollisionOperator< CollisionConfiguration >::IdxRangeDistributionFunctionType =  typename CollisionConfigurationType::IdxRangeDistributionFunctionType;
```




<hr>
## Public Functions Documentation




### function CollisionOperator 

_Construct a new Collision Operator object._ 
```C++
inline explicit CollisionOperator::CollisionOperator (
    CollisionConfigurationType const & collision_configuration
) 
```





**Parameters:**


* `collision_configuration` This parametrise the operator. It contains the data required for its initialisation. 




        

<hr>



### function operator() 

_Apply the collision operator to the distribution functions of all species on all species._ 
```C++
inline void CollisionOperator::operator() (
    DField< IdxRangeDistributionFunctionType > all_f_distribution,
    double deltat_coll
) const
```





**Parameters:**


* `all_f_distribution` All the distribution function, depending on the CollisionConfigurationType, the existence of the theta, r, phi dimension may vary. At most, we have (species, phi, r, theta, vpar, mu) in layout right. 
* `deltat_coll` Collision time step. 




        

<hr>



### function ~CollisionOperator 

_Destroy the Collision operator object._ 
```C++
inline CollisionOperator::~CollisionOperator () 
```




<hr>
## Protected Attributes Documentation




### variable m\_operator\_handle 

_Opaque C type representing the operator._ 
```C++
::koliop_Operator CollisionOperator< CollisionConfiguration >::m_operator_handle;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/collisions/collision_operator.hpp`

