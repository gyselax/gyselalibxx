

# Class CollisionConfiguration

**template &lt;class GridSp, class [**GridVpar**](structGridVpar.md), class [**GridMu**](structGridMu.md)&gt;**



[**ClassList**](annotated.md) **>** [**CollisionConfiguration**](classCollisionConfiguration.md)



_Class to collect information to initialise the collision operator for a SpVparMu geometry._ [More...](#detailed-description)

* `#include <collision_configuration.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef detail::CollisionConfigurationData&lt; [**IdxRangeDistributionFunctionType**](classCollisionConfiguration.md#typedef-idxrangedistributionfunctiontype), [**GridSpType**](classCollisionConfiguration.md#typedef-gridsptype), [**GridPhiType**](classCollisionConfiguration.md#typedef-gridphitype), [**GridRType**](classCollisionConfiguration.md#typedef-gridrtype), [**GridThetaType**](classCollisionConfiguration.md#typedef-gridthetatype), [**GridVparType**](classCollisionConfiguration.md#typedef-gridvpartype), [**GridMuType**](classCollisionConfiguration.md#typedef-gridmutype) &gt; | [**CollisionConfigurationDataType**](#typedef-collisionconfigurationdatatype)  <br>_Container for the operator input data._  |
| typedef [**GridMu**](structGridMu.md) | [**GridMuType**](#typedef-gridmutype)  <br>[_**Mu**_](structMu.md) _grid._ |
| typedef detail::InternalSpoofGridPhi | [**GridPhiType**](#typedef-gridphitype)  <br>_Phi grid._  |
| typedef detail::InternalSpoofGridR | [**GridRType**](#typedef-gridrtype)  <br>[_**R**_](structR.md) _grid._ |
| typedef GridSp | [**GridSpType**](#typedef-gridsptype)  <br>_Sp grid._  |
| typedef detail::InternalSpoofGridTheta | [**GridThetaType**](#typedef-gridthetatype)  <br>[_**Theta**_](structTheta.md) _grid._ |
| typedef [**GridVpar**](structGridVpar.md) | [**GridVparType**](#typedef-gridvpartype)  <br>[_**Vpar**_](structVpar.md) _grid._ |
| typedef IdxRange&lt; [**GridSpType**](classCollisionConfiguration.md#typedef-gridsptype), [**GridVparType**](classCollisionConfiguration.md#typedef-gridvpartype), [**GridMuType**](classCollisionConfiguration.md#typedef-gridmutype) &gt; | [**IdxRangeDistributionFunctionType**](#typedef-idxrangedistributionfunctiontype)  <br>_Distribution function index range._  |
| typedef IdxRange&lt; [**GridMuType**](classCollisionConfiguration.md#typedef-gridmutype) &gt; | [**IdxRangeMuType**](#typedef-idxrangemutype)  <br>[_**Mu**_](structMu.md) _index range._ |
| typedef IdxRange&lt; [**GridSpType**](classCollisionConfiguration.md#typedef-gridsptype), [**GridVparType**](classCollisionConfiguration.md#typedef-gridvpartype) &gt; | [**IdxRangeSpVparType**](#typedef-idxrangespvpartype)  <br>_Sp,_ [_**Vpar**_](structVpar.md) _index range._ |
| typedef IdxRange&lt; [**GridVparType**](classCollisionConfiguration.md#typedef-gridvpartype) &gt; | [**IdxRangeVparType**](#typedef-idxrangevpartype)  <br>[_**Vpar**_](structVpar.md) _index range._ |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CollisionConfiguration**](#function-collisionconfiguration) (PC\_tree\_t const & yaml\_input\_file, [**IdxRangeDistributionFunctionType**](classCollisionConfiguration.md#typedef-idxrangedistributionfunctiontype) index\_range\_fdistribution, DConstField&lt; [**IdxRangeMuType**](classCollisionConfiguration.md#typedef-idxrangemutype) &gt; coeff\_intdmu, DConstField&lt; [**IdxRangeVparType**](classCollisionConfiguration.md#typedef-idxrangevpartype) &gt; coeff\_intdvpar, double B\_norm, DConstField&lt; [**IdxRangeSpVparType**](classCollisionConfiguration.md#typedef-idxrangespvpartype) &gt; Bstar\_s) <br>_The constructor for the_ [_**CollisionConfiguration**_](classCollisionConfiguration.md) _class._ |
|  const [**CollisionConfigurationDataType**](classCollisionConfiguration.md#typedef-collisionconfigurationdatatype) & | [**configuration**](#function-configuration) () const<br>_Can be used to obtain the quantities the operator needs for its initialisation, modification is not authorised._  |








## Protected Attributes

| Type | Name |
| ---: | :--- |
|  double | [**m\_nustar0**](#variable-m_nustar0)  <br>_Value of nustar0 = to nustar0\_rpeak read in the YAML input file._  |
|  [**CollisionConfigurationDataType**](classCollisionConfiguration.md#typedef-collisionconfigurationdatatype) | [**m\_operator\_quantities**](#variable-m_operator_quantities)  <br>_Container for the arrays the quantities the operator needs._  |
















## Protected Functions

| Type | Name |
| ---: | :--- |
|  void | [**do\_configuration\_data\_initialisation**](#function-do_configuration_data_initialisation) (double input\_B\_norm, DConstField&lt; [**IdxRangeSpVparType**](classCollisionConfiguration.md#typedef-idxrangespvpartype) &gt; input\_Bstar\_s) <br>_Does the initialisation specific to this geometry._  |




## Detailed Description


NOTE: Thanks to C++17 template constructor argument type deduction, we do not need to give the 3 template arguments. Giving a correct index\_range\_fdistribution is sufficient. 


    
## Public Types Documentation




### typedef CollisionConfigurationDataType 

_Container for the operator input data._ 
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::CollisionConfigurationDataType =  detail::CollisionConfigurationData< IdxRangeDistributionFunctionType, GridSpType, GridPhiType, GridRType, GridThetaType, GridVparType, GridMuType>;
```




<hr>



### typedef GridMuType 

[_**Mu**_](structMu.md) _grid._
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::GridMuType =  GridMu;
```




<hr>



### typedef GridPhiType 

_Phi grid._ 
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::GridPhiType =  detail::InternalSpoofGridPhi;
```




<hr>



### typedef GridRType 

[_**R**_](structR.md) _grid._
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::GridRType =  detail::InternalSpoofGridR;
```




<hr>



### typedef GridSpType 

_Sp grid._ 
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::GridSpType =  GridSp;
```




<hr>



### typedef GridThetaType 

[_**Theta**_](structTheta.md) _grid._
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::GridThetaType =  detail::InternalSpoofGridTheta;
```




<hr>



### typedef GridVparType 

[_**Vpar**_](structVpar.md) _grid._
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::GridVparType =  GridVpar;
```




<hr>



### typedef IdxRangeDistributionFunctionType 

_Distribution function index range._ 
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::IdxRangeDistributionFunctionType =  IdxRange<GridSpType, GridVparType, GridMuType>;
```




<hr>



### typedef IdxRangeMuType 

[_**Mu**_](structMu.md) _index range._
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::IdxRangeMuType =  IdxRange<GridMuType>;
```




<hr>



### typedef IdxRangeSpVparType 

_Sp,_ [_**Vpar**_](structVpar.md) _index range._
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::IdxRangeSpVparType =  IdxRange<GridSpType, GridVparType>;
```




<hr>



### typedef IdxRangeVparType 

[_**Vpar**_](structVpar.md) _index range._
```C++
using CollisionConfiguration< GridSp, GridVpar, GridMu >::IdxRangeVparType =  IdxRange<GridVparType>;
```




<hr>
## Public Functions Documentation




### function CollisionConfiguration 

_The constructor for the_ [_**CollisionConfiguration**_](classCollisionConfiguration.md) _class._
```C++
inline CollisionConfiguration::CollisionConfiguration (
    PC_tree_t const & yaml_input_file,
    IdxRangeDistributionFunctionType index_range_fdistribution,
    DConstField< IdxRangeMuType > coeff_intdmu,
    DConstField< IdxRangeVparType > coeff_intdvpar,
    double B_norm,
    DConstField< IdxRangeSpVparType > Bstar_s
) 
```





**Parameters:**


* `yaml_input_file` yaml namelist file object. 
* `index_range_fdistribution` index range of the whole distribution function. 
* `coeff_intdmu` quadrature coefficients. 
* `coeff_intdvpar` quadrature coefficients. 
* `B_norm` The norm of the magnetic field defined over the polar slice. 
* `Bstar_s` Bstar defined on the coordinates (species,r,theta,vpar). 




        

<hr>



### function configuration 

_Can be used to obtain the quantities the operator needs for its initialisation, modification is not authorised._ 
```C++
inline const CollisionConfigurationDataType & CollisionConfiguration::configuration () const
```





**Returns:**

const CollisionConfigurationDataType& 





        

<hr>
## Protected Attributes Documentation




### variable m\_nustar0 

_Value of nustar0 = to nustar0\_rpeak read in the YAML input file._ 
```C++
double CollisionConfiguration< GridSp, GridVpar, GridMu >::m_nustar0;
```




<hr>



### variable m\_operator\_quantities 

_Container for the arrays the quantities the operator needs._ 
```C++
CollisionConfigurationDataType CollisionConfiguration< GridSp, GridVpar, GridMu >::m_operator_quantities;
```




<hr>
## Protected Functions Documentation




### function do\_configuration\_data\_initialisation 

_Does the initialisation specific to this geometry._ 
```C++
inline void CollisionConfiguration::do_configuration_data_initialisation (
    double input_B_norm,
    DConstField< IdxRangeSpVparType > input_Bstar_s
) 
```





**Parameters:**


* `input_B_norm` 
* `input_Bstar_s` 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryVparMu/collisions/collision_configuration.hpp`

