

# Class CollisionOperator

**template &lt;class CollisionConfiguration&gt;**



[**ClassList**](annotated.md) **>** [**CollisionOperator**](classCollisionOperator.md)



_A class which computes the collision operator in (Sp,vpar,mu)._ 

* `#include <collision_operator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename CollisionConfiguration::CollisionConfigurationDataType | [**CollisionConfigurationDataType**](#typedef-collisionconfigurationdatatype)  <br>[_**CollisionConfigurationDataType**_](classCollisionOperator.md#typedef-collisionconfigurationdatatype) _._ |
| typedef typename CollisionConfiguration::IdxRangeDistributionFunctionType | [**IdxRangeDistributionFunctionType**](#typedef-idxrangedistributionfunctiontype)  <br>[_**IdxRangeDistributionFunctionType**_](classCollisionOperator.md#typedef-idxrangedistributionfunctiontype) _._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CollisionOperator**](#function-collisionoperator) (CollisionConfiguration const & collision\_configuration) <br>_Construct a new Collision Operator object._  |
|  void | [**operator()**](#function-operator) (DField&lt; [**IdxRangeDistributionFunctionType**](classCollisionOperator.md#typedef-idxrangedistributionfunctiontype) &gt; all\_f\_distribution, double deltat\_coll) <br>_Apply the collision operator to the distribution functions of all species on all species._  |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  CollisionConfiguration const  & | [**m\_collision\_configuration**](#variable-m_collision_configuration)  <br>_The collision configuration._  |
|  koliop::collisions::Operator | [**m\_operator**](#variable-m_operator)  <br>_Opaque C type representing the operator._  |




















## Public Types Documentation




### typedef CollisionConfigurationDataType 

[_**CollisionConfigurationDataType**_](classCollisionOperator.md#typedef-collisionconfigurationdatatype) _._
```C++
using CollisionOperator< CollisionConfiguration >::CollisionConfigurationDataType = 
            typename CollisionConfiguration::CollisionConfigurationDataType;
```




<hr>



### typedef IdxRangeDistributionFunctionType 

[_**IdxRangeDistributionFunctionType**_](classCollisionOperator.md#typedef-idxrangedistributionfunctiontype) _._
```C++
using CollisionOperator< CollisionConfiguration >::IdxRangeDistributionFunctionType = 
            typename CollisionConfiguration::IdxRangeDistributionFunctionType;
```




<hr>
## Public Functions Documentation




### function CollisionOperator 

_Construct a new Collision Operator object._ 
```C++
inline explicit CollisionOperator::CollisionOperator (
    CollisionConfiguration const & collision_configuration
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
) 
```





**Parameters:**


* `all_f_distribution` All the distribution function, depending on the CollisionConfiguration, the existence of the special dimension may vary. At most, we have (species, tor3, tor1, tor2, vpar, mu) in layout right. 
* `deltat_coll` Collision time step. 




        

<hr>
## Protected Attributes Documentation




### variable m\_collision\_configuration 

_The collision configuration._ 
```C++
CollisionConfiguration const& CollisionOperator< CollisionConfiguration >::m_collision_configuration;
```




<hr>



### variable m\_operator 

_Opaque C type representing the operator._ 
```C++
koliop::collisions::Operator CollisionOperator< CollisionConfiguration >::m_operator;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/collisions/collision_operator.hpp`

