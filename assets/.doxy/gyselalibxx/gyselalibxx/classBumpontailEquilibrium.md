

# Class BumpontailEquilibrium



[**ClassList**](annotated.md) **>** [**BumpontailEquilibrium**](classBumpontailEquilibrium.md)



_A class that initialises the distribution function as a sum of two Maxwellian functions._ [More...](#detailed-description)

* `#include <bumpontailequilibrium.hpp>`



Inherits the following classes: [IEquilibrium](classIEquilibrium.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**BumpontailEquilibrium**](#function-bumpontailequilibrium) (host\_t&lt; DFieldMemSp &gt; epsilon\_bot, host\_t&lt; DFieldMemSp &gt; temperature\_bot, host\_t&lt; DFieldMemSp &gt; mean\_velocity\_bot) <br>_Creates an instance of the_ [_**BumpontailEquilibrium**_](classBumpontailEquilibrium.md) _class._ |
|  void | [**compute\_twomaxwellian**](#function-compute_twomaxwellian) (DFieldVx fMaxwellian, double epsilon\_bot, double temperature\_bot, double mean\_velocity\_bot) const<br>_Compute a distribution function defined as a sum of two Maxwellians. This distribution function can be written as $f(x,v) = f1(v) + f2(v) $ with $f1(v) = (1-epsilon)/(sqrt(2\*PI))\*exp(-v\*\*2/2)$ $f2(v) = epsilon/sqrt(2\*PI\*T0)[exp(-(v-v0)\*\*2/2\*T0)$._  |
|  host\_t&lt; DConstFieldSp &gt; | [**epsilon\_bot**](#function-epsilon_bot) () const<br>_A method for accessing the m\_epsilon\_bot member variable of the class._  |
|  host\_t&lt; DConstFieldSp &gt; | [**mean\_velocity\_bot**](#function-mean_velocity_bot) () const<br>_A method for accessing the m\_mean\_velocity\_bot member variable of the class._  |
| virtual DFieldSpVx | [**operator()**](#function-operator) (DFieldSpVx allfequilibrium) override const<br>_Initialises the distribution function as the sum of a bulk and a bump-on-tail Maxwellians._  |
|  host\_t&lt; DConstFieldSp &gt; | [**temperature\_bot**](#function-temperature_bot) () const<br>_A method for accessing the m\_temperature\_bot member variable of the class._  |
|   | [**~BumpontailEquilibrium**](#function-bumpontailequilibrium) () override<br> |


## Public Functions inherited from IEquilibrium

See [IEquilibrium](classIEquilibrium.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVparMu | [**operator()**](classIEquilibrium.md#function-operator) (DFieldSpVparMu allfequilibrium) const = 0<br>_Operator for initialising an equilibrium distribution function._  |
| virtual DFieldSpVx | [**operator()**](classIEquilibrium.md#function-operator_1) (DFieldSpVx allfequilibrium) const = 0<br>_Operator for initialising a distribution function that does not depend on space._  |
| virtual DFieldSpVxVy | [**operator()**](classIEquilibrium.md#function-operator_2) (DFieldSpVxVy allfequilibrium) const = 0<br>_Operator for initialising a distribution function that does not depend on space._  |
| virtual  | [**~IEquilibrium**](classIEquilibrium.md#function-iequilibrium-13) () = default<br> |
| virtual  | [**~IEquilibrium**](classIEquilibrium.md#function-iequilibrium-13) () = default<br> |
| virtual  | [**~IEquilibrium**](classIEquilibrium.md#function-iequilibrium-13) () = default<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  [**BumpontailEquilibrium**](classBumpontailEquilibrium.md) | [**init\_from\_input**](#function-init_from_input) (IdxRangeSp idx\_range\_kinsp, PC\_tree\_t const & yaml\_input\_file) <br>_Read the density, temperature and mean velocity required to initialise the bump-on-tail Maxwellian in a YAML input file._  |




















































## Detailed Description


This class initialises the distribution function as a sum of two Maxwellian, enabling the study of the so-called bump-on-tail instability. One of the Maxwellians represents the bulk of the distribution function that has no mean velocity, and the other Maxwellian corresponds to high velocity particles. The second Maxwellian is referred to as the "bump-on-tail" Maxwellian. 


    
## Public Functions Documentation




### function BumpontailEquilibrium 

_Creates an instance of the_ [_**BumpontailEquilibrium**_](classBumpontailEquilibrium.md) _class._
```C++
BumpontailEquilibrium::BumpontailEquilibrium (
    host_t< DFieldMemSp > epsilon_bot,
    host_t< DFieldMemSp > temperature_bot,
    host_t< DFieldMemSp > mean_velocity_bot
) 
```





**Parameters:**


* `epsilon_bot` A parameter that represents the density of the bump-on-tail Maxwellian for each species. 
* `temperature_bot` A parameter that represents the temperature of the bump-on-tail Maxwellian for each species. 
* `mean_velocity_bot` A parameter that represents the mean velocity of the bump-on-tail Maxwellian for each species. 




        

<hr>



### function compute\_twomaxwellian 

_Compute a distribution function defined as a sum of two Maxwellians. This distribution function can be written as $f(x,v) = f1(v) + f2(v) $ with $f1(v) = (1-epsilon)/(sqrt(2\*PI))\*exp(-v\*\*2/2)$ $f2(v) = epsilon/sqrt(2\*PI\*T0)[exp(-(v-v0)\*\*2/2\*T0)$._ 
```C++
void BumpontailEquilibrium::compute_twomaxwellian (
    DFieldVx fMaxwellian,
    double epsilon_bot,
    double temperature_bot,
    double mean_velocity_bot
) const
```





**Parameters:**


* `fMaxwellian` The initial distribution function. 
* `epsilon_bot` A parameter that represents the density of the bump-on-tail Maxwellian. 
* `temperature_bot` A parameter that represents the temperature of the bump-on-tail Maxwellian. 
* `mean_velocity_bot` A parameter that represents the mean velocity of the bump-on-tail Maxwellian. 




        

<hr>



### function epsilon\_bot 

_A method for accessing the m\_epsilon\_bot member variable of the class._ 
```C++
inline host_t< DConstFieldSp > BumpontailEquilibrium::epsilon_bot () const
```





**Returns:**

a field containing the m\_epsilon\_bot variable. 





        

<hr>



### function mean\_velocity\_bot 

_A method for accessing the m\_mean\_velocity\_bot member variable of the class._ 
```C++
inline host_t< DConstFieldSp > BumpontailEquilibrium::mean_velocity_bot () const
```





**Returns:**

a field containing the m\_velocity\_bot variable. 





        

<hr>



### function operator() 

_Initialises the distribution function as the sum of a bulk and a bump-on-tail Maxwellians._ 
```C++
virtual DFieldSpVx BumpontailEquilibrium::operator() (
    DFieldSpVx allfequilibrium
) override const
```





**Parameters:**


* `allfequilibrium` The initialised distribution function. 



**Returns:**

The initialised distribution function. 





        
Implements [*IEquilibrium::operator()*](classIEquilibrium.md#function-operator_1)


<hr>



### function temperature\_bot 

_A method for accessing the m\_temperature\_bot member variable of the class._ 
```C++
inline host_t< DConstFieldSp > BumpontailEquilibrium::temperature_bot () const
```





**Returns:**

a field containing the m\_temperature\_bot variable. 





        

<hr>



### function ~BumpontailEquilibrium 

```C++
BumpontailEquilibrium::~BumpontailEquilibrium () override
```




<hr>
## Public Static Functions Documentation




### function init\_from\_input 

_Read the density, temperature and mean velocity required to initialise the bump-on-tail Maxwellian in a YAML input file._ 
```C++
static BumpontailEquilibrium BumpontailEquilibrium::init_from_input (
    IdxRangeSp idx_range_kinsp,
    PC_tree_t const & yaml_input_file
) 
```





**Parameters:**


* `idx_range_kinsp` Index range for the kinetic species 
* `yaml_input_file` YAML input file 



**Returns:**

an instance of Maxwellian distribution function. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/initialisation/bumpontailequilibrium.hpp`

