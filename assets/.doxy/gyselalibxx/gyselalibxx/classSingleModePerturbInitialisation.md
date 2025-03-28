

# Class SingleModePerturbInitialisation



[**ClassList**](annotated.md) **>** [**SingleModePerturbInitialisation**](classSingleModePerturbInitialisation.md)



_A class that initialises the distribution function as a perturbed Maxwellian._ [More...](#detailed-description)

* `#include <singlemodeperturbinitialisation.hpp>`



Inherits the following classes: [IInitialisation](classIInitialisation.md),  [IInitialisation](classIInitialisation.md)










































































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SingleModePerturbInitialisation**](#function-singlemodeperturbinitialisation-12) (DConstFieldSpVx fequilibrium, host\_t&lt; IFieldMemSp &gt; init\_perturb\_mode, host\_t&lt; DFieldMemSp &gt; init\_perturb\_amplitude) <br>_Creates an instance of the_ [_**SingleModePerturbInitialisation**_](classSingleModePerturbInitialisation.md) _class._ |
|   | [**SingleModePerturbInitialisation**](#function-singlemodeperturbinitialisation-22) (DConstFieldSpVxVy fequilibrium, host\_t&lt; IFieldMemSp &gt; init\_perturb\_mode, host\_t&lt; DFieldMemSp &gt; init\_perturb\_amplitude) <br>_Creates an instance of the_ [_**SingleModePerturbInitialisation**_](classSingleModePerturbInitialisation.md) _class._ |
| virtual DFieldSpXVx | [**operator()**](#function-operator) (DFieldSpXVx allfdistribu) override const<br>_Initialises the distribution function as as a perturbed Maxwellian._  |
| virtual DFieldSpXYVxVy | [**operator()**](#function-operator_1) (DFieldSpXYVxVy allfdistribu) override const<br>_Initialises the distribution function as as a perturbed Maxwellian._  |
|  void | [**perturbation\_initialisation**](#function-perturbation_initialisation-12) (DFieldX perturbation, int const perturb\_mode, double const perturb\_amplitude) const<br>_Initialisation of the perturbation._  |
|  void | [**perturbation\_initialisation**](#function-perturbation_initialisation-22) (DFieldXY perturbation, int const perturb\_mode, double const perturb\_amplitude) const<br>_Initialisation of the perturbation._  |
|   | [**~SingleModePerturbInitialisation**](#function-singlemodeperturbinitialisation-12) () override<br> |
|   | [**~SingleModePerturbInitialisation**](#function-singlemodeperturbinitialisation-12) () override<br> |


## Public Functions inherited from IInitialisation

See [IInitialisation](classIInitialisation.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVparMu | [**operator()**](classIInitialisation.md#function-operator) (DFieldSpVparMu allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual DFieldSpXVx | [**operator()**](classIInitialisation.md#function-operator_1) (DFieldSpXVx allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual DFieldSpXYVxVy | [**operator()**](classIInitialisation.md#function-operator_2) (DFieldSpXYVxVy allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |


## Public Functions inherited from IInitialisation

See [IInitialisation](classIInitialisation.md)

| Type | Name |
| ---: | :--- |
| virtual DFieldSpVparMu | [**operator()**](classIInitialisation.md#function-operator) (DFieldSpVparMu allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual DFieldSpXVx | [**operator()**](classIInitialisation.md#function-operator_1) (DFieldSpXVx allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual DFieldSpXYVxVy | [**operator()**](classIInitialisation.md#function-operator_2) (DFieldSpXYVxVy allfdistribu) const = 0<br>_Operator for initialising a distribution function._  |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |
| virtual  | [**~IInitialisation**](classIInitialisation.md#function-iinitialisation-13) () = default<br> |


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  [**SingleModePerturbInitialisation**](classSingleModePerturbInitialisation.md) | [**init\_from\_input**](#function-init_from_input-12) (DConstFieldSpVx allfequilibrium, IdxRangeSp idx\_range\_kinsp, PC\_tree\_t const & yaml\_input\_file) <br>_Read init\_perturb\_mode and init\_perturb amplitude in a YAML input file to initialise the perturbation._  |
|  [**SingleModePerturbInitialisation**](classSingleModePerturbInitialisation.md) | [**init\_from\_input**](#function-init_from_input-22) (DConstFieldSpVxVy allfequilibrium, IdxRangeSp idx\_range\_kinsp, PC\_tree\_t const & yaml\_input\_file) <br>_Read init\_perturb\_mode and init\_perturb amplitude in a YAML input file to initialise the perturbation._  |














































































## Detailed Description


Initialisation operator with a sinusoidal perturbation of a Maxwellian. This initialises all species.


A class that initialises the distribution function as a perturbed Maxwellian defined as $f = f\_{maxw}(v) \* (1 + perturb(x))$, where $f\_{maxw}(v)$ is a Maxwellian, and $perturb(x)$ is a sinusoidal perturbation. 


    
## Public Functions Documentation




### function SingleModePerturbInitialisation [1/2]

_Creates an instance of the_ [_**SingleModePerturbInitialisation**_](classSingleModePerturbInitialisation.md) _class._
```C++
SingleModePerturbInitialisation::SingleModePerturbInitialisation (
    DConstFieldSpVx fequilibrium,
    host_t< IFieldMemSp > init_perturb_mode,
    host_t< DFieldMemSp > init_perturb_amplitude
) 
```





**Parameters:**


* `fequilibrium` A Maxwellian. 
* `init_perturb_mode` The perturbation mode. 
* `init_perturb_amplitude` The perturbation amplitude. 




        

<hr>



### function SingleModePerturbInitialisation [2/2]

_Creates an instance of the_ [_**SingleModePerturbInitialisation**_](classSingleModePerturbInitialisation.md) _class._
```C++
SingleModePerturbInitialisation::SingleModePerturbInitialisation (
    DConstFieldSpVxVy fequilibrium,
    host_t< IFieldMemSp > init_perturb_mode,
    host_t< DFieldMemSp > init_perturb_amplitude
) 
```





**Parameters:**


* `fequilibrium` A Maxwellian. 
* `init_perturb_mode` The perturbation mode. 
* `init_perturb_amplitude` The perturbation amplitude. 




        

<hr>



### function operator() 

_Initialises the distribution function as as a perturbed Maxwellian._ 
```C++
virtual DFieldSpXVx SingleModePerturbInitialisation::operator() (
    DFieldSpXVx allfdistribu
) override const
```





**Parameters:**


* `allfdistribu` The initialised distribution function. 



**Returns:**

The initialised distribution function. 





        
Implements [*IInitialisation::operator()*](classIInitialisation.md#function-operator_1)


<hr>



### function operator() 

_Initialises the distribution function as as a perturbed Maxwellian._ 
```C++
virtual DFieldSpXYVxVy SingleModePerturbInitialisation::operator() (
    DFieldSpXYVxVy allfdistribu
) override const
```





**Parameters:**


* `allfdistribu` The initialised distribution function. 



**Returns:**

The initialised distribution function. 





        
Implements [*IInitialisation::operator()*](classIInitialisation.md#function-operator_2)


<hr>



### function perturbation\_initialisation [1/2]

_Initialisation of the perturbation._ 
```C++
void SingleModePerturbInitialisation::perturbation_initialisation (
    DFieldX perturbation,
    int const perturb_mode,
    double const perturb_amplitude
) const
```





**Parameters:**


* `perturbation` On input: an uninitialized array On output: an array containing a values that has a sinusoidal variation with given amplitude and mode. 
* `perturb_mode` The mode of the perturbation. 
* `perturb_amplitude` The amplitude of the perturbation. 




        

<hr>



### function perturbation\_initialisation [2/2]

_Initialisation of the perturbation._ 
```C++
void SingleModePerturbInitialisation::perturbation_initialisation (
    DFieldXY perturbation,
    int const perturb_mode,
    double const perturb_amplitude
) const
```





**Parameters:**


* `perturbation` On input: an uninitialized array On output: an array containing a values that has a sinusoidal variation with given amplitude and mode. 
* `perturb_mode` The mode of the perturbation. 
* `perturb_amplitude` The amplitude of the perturbation. 




        

<hr>



### function ~SingleModePerturbInitialisation [1/2]

```C++
SingleModePerturbInitialisation::~SingleModePerturbInitialisation () override
```




<hr>



### function ~SingleModePerturbInitialisation [1/2]

```C++
SingleModePerturbInitialisation::~SingleModePerturbInitialisation () override
```




<hr>
## Public Static Functions Documentation




### function init\_from\_input [1/2]

_Read init\_perturb\_mode and init\_perturb amplitude in a YAML input file to initialise the perturbation._ 
```C++
static SingleModePerturbInitialisation SingleModePerturbInitialisation::init_from_input (
    DConstFieldSpVx allfequilibrium,
    IdxRangeSp idx_range_kinsp,
    PC_tree_t const & yaml_input_file
) 
```





**Parameters:**


* `allfequilibrium` equilibrium distribution function. 
* `idx_range_kinsp` Index range for the kinetic species. 
* `yaml_input_file` YAML input file. 



**Returns:**

an instance of [**SingleModePerturbInitialisation**](classSingleModePerturbInitialisation.md) class. 





        

<hr>



### function init\_from\_input [2/2]

_Read init\_perturb\_mode and init\_perturb amplitude in a YAML input file to initialise the perturbation._ 
```C++
static SingleModePerturbInitialisation SingleModePerturbInitialisation::init_from_input (
    DConstFieldSpVxVy allfequilibrium,
    IdxRangeSp idx_range_kinsp,
    PC_tree_t const & yaml_input_file
) 
```





**Parameters:**


* `allfequilibrium` equilibrium distribution function. 
* `idx_range_kinsp` Index range for the kinetic species. 
* `yaml_input_file` YAML input file. 



**Returns:**

an instance of [**SingleModePerturbInitialisation**](classSingleModePerturbInitialisation.md) class. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/initialisation/singlemodeperturbinitialisation.hpp`

