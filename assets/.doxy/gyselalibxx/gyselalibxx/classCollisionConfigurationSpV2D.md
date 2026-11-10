

# Class CollisionConfigurationSpV2D



[**ClassList**](annotated.md) **>** [**CollisionConfigurationSpV2D**](classCollisionConfigurationSpV2D.md)



_Class to collect information to initialise the collision operator for a SpVparMu geometry._ 

* `#include <collision_configuration.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef detail::CollisionConfigurationData&lt; [**IdxRangeDistributionFunctionType**](classCollisionConfigurationSpV2D.md#typedef-idxrangedistributionfunctiontype), detail::InternalSpoofGridTor3, detail::InternalSpoofGridTor2, detail::InternalSpoofGridTor1, [**GridVpar**](structGridVpar.md), [**GridMu**](structGridMu.md) &gt; | [**CollisionConfigurationDataType**](#typedef-collisionconfigurationdatatype)  <br>_Container for the operator input data._  |
| typedef IdxRange&lt; [**Species**](structSpecies.md), [**GridVpar**](structGridVpar.md), [**GridMu**](structGridMu.md) &gt; | [**IdxRangeDistributionFunctionType**](#typedef-idxrangedistributionfunctiontype)  <br>_Distribution function index range, define the distribution function layout (layout right)._  |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CollisionConfigurationSpV2D**](#function-collisionconfigurationspv2d) (PC\_tree\_t const & yaml\_input\_file, [**IdxRangeDistributionFunctionType**](classCollisionConfigurationSpV2D.md#typedef-idxrangedistributionfunctiontype) index\_range\_fdistribution, DConstField&lt; IdxRange&lt; [**GridMu**](structGridMu.md) &gt; &gt; coeff\_intdmu, DConstField&lt; IdxRange&lt; [**GridVpar**](structGridVpar.md) &gt; &gt; coeff\_intdvpar, double B\_norm, DConstField&lt; IdxRange&lt; [**Species**](structSpecies.md), [**GridVpar**](structGridVpar.md) &gt; &gt; Bstar\_s) <br>_The constructor for the_ [_**CollisionConfigurationSpV2D**_](classCollisionConfigurationSpV2D.md) _class._ |
|  const [**CollisionConfigurationDataType**](classCollisionConfigurationSpV2D.md#typedef-collisionconfigurationdatatype) & | [**configuration**](#function-configuration) () const<br>_Can be used to obtain the quantities the operator needs for its initialisation, modification is not authorised._  |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  [**CollisionConfigurationDataType**](classCollisionConfigurationSpV2D.md#typedef-collisionconfigurationdatatype) | [**m\_operator\_quantities**](#variable-m_operator_quantities)  <br>_Container for the arrays the quantities the operator needs._  |




















## Public Types Documentation




### typedef CollisionConfigurationDataType 

_Container for the operator input data._ 
```C++
using CollisionConfigurationSpV2D::CollisionConfigurationDataType =  detail::CollisionConfigurationData<
            IdxRangeDistributionFunctionType,
            detail::InternalSpoofGridTor3,
            detail::InternalSpoofGridTor2,
            detail::InternalSpoofGridTor1,
            GridVpar,
            GridMu>;
```




<hr>



### typedef IdxRangeDistributionFunctionType 

_Distribution function index range, define the distribution function layout (layout right)._ 
```C++
using CollisionConfigurationSpV2D::IdxRangeDistributionFunctionType =  IdxRange<Species, GridVpar, GridMu>;
```




<hr>
## Public Functions Documentation




### function CollisionConfigurationSpV2D 

_The constructor for the_ [_**CollisionConfigurationSpV2D**_](classCollisionConfigurationSpV2D.md) _class._
```C++
inline CollisionConfigurationSpV2D::CollisionConfigurationSpV2D (
    PC_tree_t const & yaml_input_file,
    IdxRangeDistributionFunctionType index_range_fdistribution,
    DConstField< IdxRange< GridMu > > coeff_intdmu,
    DConstField< IdxRange< GridVpar > > coeff_intdvpar,
    double B_norm,
    DConstField< IdxRange< Species , GridVpar > > Bstar_s
) 
```





**Parameters:**


* `yaml_input_file` yaml namelist file object. 
* `index_range_fdistribution` index range of the whole distribution function. 
* `coeff_intdmu` quadrature coefficients. 
* `coeff_intdvpar` quadrature coefficients. 
* `B_norm` The norm of the magnetic field. 
* `Bstar_s` Bstar defined on the coordinates. 




        

<hr>



### function configuration 

_Can be used to obtain the quantities the operator needs for its initialisation, modification is not authorised._ 
```C++
inline const CollisionConfigurationDataType & CollisionConfigurationSpV2D::configuration () const
```





**Returns:**

const [**CollisionConfigurationDataType**](classCollisionConfigurationSpV2D.md#typedef-collisionconfigurationdatatype)& 





        

<hr>
## Protected Attributes Documentation




### variable m\_operator\_quantities 

_Container for the arrays the quantities the operator needs._ 
```C++
CollisionConfigurationDataType CollisionConfigurationSpV2D::m_operator_quantities;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryVparMu/collisions/collision_configuration.hpp`

