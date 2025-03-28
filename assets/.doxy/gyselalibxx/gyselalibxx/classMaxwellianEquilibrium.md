

# Class MaxwellianEquilibrium



[**ClassList**](annotated.md) **>** [**MaxwellianEquilibrium**](classMaxwellianEquilibrium.md)



_Equilibrium operator as Maxwellian. This initialises all species._ [More...](#detailed-description)

* `#include <maxwellianequilibrium.hpp>`



Inherits the following classes: [IEquilibrium](classIEquilibrium.md),  [IEquilibrium](classIEquilibrium.md),  [IEquilibrium](classIEquilibrium.md)






























































































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MaxwellianEquilibrium**](#function-maxwellianequilibrium-13) (host\_t&lt; DFieldMemSp &gt; mass, host\_t&lt; DFieldMemSp &gt; density\_eq, host\_t&lt; DFieldMemSp &gt; temperature\_eq, host\_t&lt; DFieldMemSp &gt; mean\_velocity\_eq, double magnetic\_field) <br>_The constructor for the_ [_**MaxwellianEquilibrium**_](classMaxwellianEquilibrium.md) _class._ |
|   | [**MaxwellianEquilibrium**](#function-maxwellianequilibrium-23) (host\_t&lt; DFieldMemSp &gt; density\_eq, host\_t&lt; DFieldMemSp &gt; temperature\_eq, host\_t&lt; DFieldMemSp &gt; mean\_velocity\_eq) <br>_The constructor for the_ [_**MaxwellianEquilibrium**_](classMaxwellianEquilibrium.md) _class._ |
|   | [**MaxwellianEquilibrium**](#function-maxwellianequilibrium-23) (host\_t&lt; DFieldMemSp &gt; density\_eq, host\_t&lt; DFieldMemSp &gt; temperature\_eq, host\_t&lt; DFieldMemSp &gt; mean\_velocity\_eq) <br>_The constructor for the_ [_**MaxwellianEquilibrium**_](classMaxwellianEquilibrium.md) _class._ |
|  host\_t&lt; DConstFieldSp &gt; | [**density\_eq**](#function-density_eq-13) () const<br>_A method for accessing the m\_density\_eq member variable of the class._  |
|  host\_t&lt; DConstFieldSp &gt; | [**density\_eq**](#function-density_eq-13) () const<br>_A method for accessing the m\_density\_eq member variable of the class._  |
|  host\_t&lt; ConstFieldSp&lt; double &gt; &gt; | [**density\_eq**](#function-density_eq-33) () const<br>_A method for accessing the m\_density\_eq member variable of the class._  |
|  host\_t&lt; DConstFieldSp &gt; | [**mass**](#function-mass) () const<br>_A method for accessing the m\_mass member variable of the class._  |
|  host\_t&lt; DConstFieldSp &gt; | [**mean\_velocity\_eq**](#function-mean_velocity_eq-13) () const<br>_A method for accessing the m\_mean\_velocity\_eq member variable of the class._  |
|  host\_t&lt; DConstFieldSp &gt; | [**mean\_velocity\_eq**](#function-mean_velocity_eq-13) () const<br>_A method for accessing the m\_mean\_velocity\_eq member variable of the class._  |
|  host\_t&lt; ConstFieldSp&lt; double &gt; &gt; | [**mean\_velocity\_eq**](#function-mean_velocity_eq-33) () const<br>_A method for accessing the m\_mean\_velocity\_eq member variable of the class._  |
| virtual DFieldSpVparMu | [**operator()**](#function-operator) (DFieldSpVparMu allfequilibrium) override const<br>_Initialises allfequilibrium as a Maxwellian._  |
| virtual DFieldSpVx | [**operator()**](#function-operator_1) (DFieldSpVx allfequilibrium) override const<br>_Initialises allfequilibrium as a Maxwellian._  |
| virtual DFieldSpVxVy | [**operator()**](#function-operator_2) (DFieldSpVxVy allfequilibrium) override const<br>_Initialises allfequilibrium as a Maxwellian._  |
|  host\_t&lt; DConstFieldSp &gt; | [**temperature\_eq**](#function-temperature_eq-13) () const<br>_A method for accessing the m\_temperature\_eq member variable of the class._  |
|  host\_t&lt; DConstFieldSp &gt; | [**temperature\_eq**](#function-temperature_eq-13) () const<br>_A method for accessing the m\_temperature\_eq member variable of the class._  |
|  host\_t&lt; ConstFieldSp&lt; double &gt; &gt; | [**temperature\_eq**](#function-temperature_eq-33) () const<br>_A method for accessing the m\_temperature\_eq member variable of the class._  |
|   | [**~MaxwellianEquilibrium**](#function-maxwellianequilibrium-13) () override<br> |
|   | [**~MaxwellianEquilibrium**](#function-maxwellianequilibrium-13) () override<br> |
|   | [**~MaxwellianEquilibrium**](#function-maxwellianequilibrium-13) () override<br> |


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
|  void | [**compute\_maxwellian**](#function-compute_maxwellian-13) (DFieldVparMu const fMaxwellian, double const mass, double const density, double const temperature, double const mean\_velocity, double const magnetic\_field) <br>_Compute a Maxwellian distribution function. The Maxwellian distribution function is defined as Compute $fM(v,mu) = (2\*PI\*T)\*\*1.5\*n\*exp(-E)$ with._  |
|  void | [**compute\_maxwellian**](#function-compute_maxwellian-23) (DFieldVx const fMaxwellian, double const density, double const temperature, double const mean\_velocity) <br>_Compute a Maxwellian distribution function. The Maxwellian distribution function is defined as $f\_M(v) = n/(sqrt(2\*PI\*T))\*exp(-(v-u)\*\*2/(2\*T))$ with $n$ the density, $T$ the temperature and $u$ is the mean velocity._  |
|  void | [**compute\_maxwellian**](#function-compute_maxwellian-33) (DFieldVxVy const fMaxwellian, double const density, double const temperature, double const mean\_velocity) <br>_Compute a Maxwellian distribution function. The Maxwellian distribution function is defined as $f\_M(v) = n/(sqrt(2\*PI\*T))\*exp(-(v-u)\*\*2/(2\*T))$ with $n$ the density, $T$ the temperature and $u$ is the mean velocity._  |
|  [**MaxwellianEquilibrium**](classMaxwellianEquilibrium.md) | [**init\_from\_input**](#function-init_from_input-13) (IdxRangeSp idx\_range\_kinsp, PC\_tree\_t const & yaml\_input\_file) <br>_Read the density, temperature and mean velocity required to initialise the Maxwellian in a YAML input file._  |
|  [**MaxwellianEquilibrium**](classMaxwellianEquilibrium.md) | [**init\_from\_input**](#function-init_from_input-23) (IdxRangeSp idx\_range\_kinsp, PC\_tree\_t const & yaml\_input\_file) <br>_Read the density, temperature and mean velocity required to initialise the Maxwellian in a YAML input file._  |
|  [**MaxwellianEquilibrium**](classMaxwellianEquilibrium.md) | [**init\_from\_input**](#function-init_from_input-23) (IdxRangeSp idx\_range\_kinsp, PC\_tree\_t const & yaml\_input\_file) <br>_Read the density, temperature and mean velocity required to initialise the Maxwellian in a YAML input file._  |








































































































## Detailed Description


A class that initialises the distribution function as a Maxwellian. 


    
## Public Functions Documentation




### function MaxwellianEquilibrium [1/3]

_The constructor for the_ [_**MaxwellianEquilibrium**_](classMaxwellianEquilibrium.md) _class._
```C++
MaxwellianEquilibrium::MaxwellianEquilibrium (
    host_t< DFieldMemSp > mass,
    host_t< DFieldMemSp > density_eq,
    host_t< DFieldMemSp > temperature_eq,
    host_t< DFieldMemSp > mean_velocity_eq,
    double magnetic_field
) 
```





**Parameters:**


* `mass` The mass of the species 
* `density_eq` The density of the Maxwellian 
* `temperature_eq` The temperature of the Maxwellian 
* `mean_velocity_eq` The mean velocity of the Maxwellian 
* `magnetic_field` The magnetic field 




        

<hr>



### function MaxwellianEquilibrium [2/3]

_The constructor for the_ [_**MaxwellianEquilibrium**_](classMaxwellianEquilibrium.md) _class._
```C++
MaxwellianEquilibrium::MaxwellianEquilibrium (
    host_t< DFieldMemSp > density_eq,
    host_t< DFieldMemSp > temperature_eq,
    host_t< DFieldMemSp > mean_velocity_eq
) 
```





**Parameters:**


* `density_eq` The density of the Maxwellian 
* `temperature_eq` The temperature of the Maxwellian 
* `mean_velocity_eq` The mean velocity of the Maxwellian 




        

<hr>



### function MaxwellianEquilibrium [2/3]

_The constructor for the_ [_**MaxwellianEquilibrium**_](classMaxwellianEquilibrium.md) _class._
```C++
MaxwellianEquilibrium::MaxwellianEquilibrium (
    host_t< DFieldMemSp > density_eq,
    host_t< DFieldMemSp > temperature_eq,
    host_t< DFieldMemSp > mean_velocity_eq
) 
```





**Parameters:**


* `density_eq` The density of the Maxwellian 
* `temperature_eq` The temperature of the Maxwellian 
* `mean_velocity_eq` The mean velocity of the Maxwellian 




        

<hr>



### function density\_eq [1/3]

_A method for accessing the m\_density\_eq member variable of the class._ 
```C++
inline host_t< DConstFieldSp > MaxwellianEquilibrium::density_eq () const
```





**Returns:**

A field containing the m\_density\_eq value. 





        

<hr>



### function density\_eq [1/3]

_A method for accessing the m\_density\_eq member variable of the class._ 
```C++
inline host_t< DConstFieldSp > MaxwellianEquilibrium::density_eq () const
```





**Returns:**

A view containing the m\_density\_eq value. 





        

<hr>



### function density\_eq [3/3]

_A method for accessing the m\_density\_eq member variable of the class._ 
```C++
inline host_t< ConstFieldSp< double > > MaxwellianEquilibrium::density_eq () const
```





**Returns:**

A field containing the m\_density\_eq value. 





        

<hr>



### function mass 

_A method for accessing the m\_mass member variable of the class._ 
```C++
inline host_t< DConstFieldSp > MaxwellianEquilibrium::mass () const
```





**Returns:**

A field containing the m\_mass value. 





        

<hr>



### function mean\_velocity\_eq [1/3]

_A method for accessing the m\_mean\_velocity\_eq member variable of the class._ 
```C++
inline host_t< DConstFieldSp > MaxwellianEquilibrium::mean_velocity_eq () const
```





**Returns:**

A field containing the m\_velocity\_eq value. 





        

<hr>



### function mean\_velocity\_eq [1/3]

_A method for accessing the m\_mean\_velocity\_eq member variable of the class._ 
```C++
inline host_t< DConstFieldSp > MaxwellianEquilibrium::mean_velocity_eq () const
```





**Returns:**

A view containing the m\_velocity\_eq value. 





        

<hr>



### function mean\_velocity\_eq [3/3]

_A method for accessing the m\_mean\_velocity\_eq member variable of the class._ 
```C++
inline host_t< ConstFieldSp< double > > MaxwellianEquilibrium::mean_velocity_eq () const
```





**Returns:**

A field containing the m\_velocity\_eq value. 





        

<hr>



### function operator() 

_Initialises allfequilibrium as a Maxwellian._ 
```C++
virtual DFieldSpVparMu MaxwellianEquilibrium::operator() (
    DFieldSpVparMu allfequilibrium
) override const
```





**Parameters:**


* `allfequilibrium` A Field containing a Maxwellian distribution function. 



**Returns:**

A Field containing a Maxwellian distribution function. 





        
Implements [*IEquilibrium::operator()*](classIEquilibrium.md#function-operator)


<hr>



### function operator() 

_Initialises allfequilibrium as a Maxwellian._ 
```C++
virtual DFieldSpVx MaxwellianEquilibrium::operator() (
    DFieldSpVx allfequilibrium
) override const
```





**Parameters:**


* `allfequilibrium` A Field containing a Maxwellian distribution function. 



**Returns:**

A Field containing a Maxwellian distribution function. 





        
Implements [*IEquilibrium::operator()*](classIEquilibrium.md#function-operator_1)


<hr>



### function operator() 

_Initialises allfequilibrium as a Maxwellian._ 
```C++
virtual DFieldSpVxVy MaxwellianEquilibrium::operator() (
    DFieldSpVxVy allfequilibrium
) override const
```





**Parameters:**


* `allfequilibrium` A Field containing a Maxwellian distribution function. 



**Returns:**

A Field containing a Maxwellian distribution function. 





        
Implements [*IEquilibrium::operator()*](classIEquilibrium.md#function-operator_2)


<hr>



### function temperature\_eq [1/3]

_A method for accessing the m\_temperature\_eq member variable of the class._ 
```C++
inline host_t< DConstFieldSp > MaxwellianEquilibrium::temperature_eq () const
```





**Returns:**

A field containing the m\_temperature\_eq value. 





        

<hr>



### function temperature\_eq [1/3]

_A method for accessing the m\_temperature\_eq member variable of the class._ 
```C++
inline host_t< DConstFieldSp > MaxwellianEquilibrium::temperature_eq () const
```





**Returns:**

A view containing the m\_temperature\_eq value. 





        

<hr>



### function temperature\_eq [3/3]

_A method for accessing the m\_temperature\_eq member variable of the class._ 
```C++
inline host_t< ConstFieldSp< double > > MaxwellianEquilibrium::temperature_eq () const
```





**Returns:**

A field containing the m\_temperature\_eq value. 





        

<hr>



### function ~MaxwellianEquilibrium [1/3]

```C++
MaxwellianEquilibrium::~MaxwellianEquilibrium () override
```




<hr>



### function ~MaxwellianEquilibrium [1/3]

```C++
MaxwellianEquilibrium::~MaxwellianEquilibrium () override
```




<hr>



### function ~MaxwellianEquilibrium [1/3]

```C++
MaxwellianEquilibrium::~MaxwellianEquilibrium () override
```




<hr>
## Public Static Functions Documentation




### function compute\_maxwellian [1/3]

_Compute a Maxwellian distribution function. The Maxwellian distribution function is defined as Compute $fM(v,mu) = (2\*PI\*T)\*\*1.5\*n\*exp(-E)$ with._ 
```C++
static void MaxwellianEquilibrium::compute_maxwellian (
    DFieldVparMu const fMaxwellian,
    double const mass,
    double const density,
    double const temperature,
    double const mean_velocity,
    double const magnetic_field
) 
```




* $n$ the density, $T$ the temperature and $u$ the mean velocity
* $B$ the magnetic field and
* $E$ the energy defined as $E = (0.5\*(v-u)\*\*2+mu\*B)/T$. 

**Parameters:**


  * `fMaxwellian` A Maxwellian distribution function. 
  * `mass` Mass of the species. 
  * `density` A parameter that represents the density of Maxwellian. 
  * `temperature` A parameter that represents the temperature of Maxwellian. 
  * `mean_velocity` A parameter that represents the mean velocity of Maxwellian. 
  * `magnetic_field` Magnetic field. 






        

<hr>



### function compute\_maxwellian [2/3]

_Compute a Maxwellian distribution function. The Maxwellian distribution function is defined as $f\_M(v) = n/(sqrt(2\*PI\*T))\*exp(-(v-u)\*\*2/(2\*T))$ with $n$ the density, $T$ the temperature and $u$ is the mean velocity._ 
```C++
static void MaxwellianEquilibrium::compute_maxwellian (
    DFieldVx const fMaxwellian,
    double const density,
    double const temperature,
    double const mean_velocity
) 
```





**Parameters:**


* `fMaxwellian` A Maxwellian distribution function. 
* `density` A parameter that represents the density of Maxwellian. 
* `temperature` A parameter that represents the temperature of Maxwellian. 
* `mean_velocity` A parameter that represents the mean velocity of Maxwellian. 




        

<hr>



### function compute\_maxwellian [3/3]

_Compute a Maxwellian distribution function. The Maxwellian distribution function is defined as $f\_M(v) = n/(sqrt(2\*PI\*T))\*exp(-(v-u)\*\*2/(2\*T))$ with $n$ the density, $T$ the temperature and $u$ is the mean velocity._ 
```C++
static void MaxwellianEquilibrium::compute_maxwellian (
    DFieldVxVy const fMaxwellian,
    double const density,
    double const temperature,
    double const mean_velocity
) 
```





**Parameters:**


* `fMaxwellian` A Maxwellian distribution function. 
* `density` A parameter that represents the density of Maxwellian. 
* `temperature` A parameter that represents the temperature of Maxwellian. 
* `mean_velocity` A parameter that represents the mean velocity of Maxwellian. 




        

<hr>



### function init\_from\_input [1/3]

_Read the density, temperature and mean velocity required to initialise the Maxwellian in a YAML input file._ 
```C++
static MaxwellianEquilibrium MaxwellianEquilibrium::init_from_input (
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



### function init\_from\_input [2/3]

_Read the density, temperature and mean velocity required to initialise the Maxwellian in a YAML input file._ 
```C++
static MaxwellianEquilibrium MaxwellianEquilibrium::init_from_input (
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



### function init\_from\_input [2/3]

_Read the density, temperature and mean velocity required to initialise the Maxwellian in a YAML input file._ 
```C++
static MaxwellianEquilibrium MaxwellianEquilibrium::init_from_input (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryVparMu/initialisation/maxwellianequilibrium.hpp`

