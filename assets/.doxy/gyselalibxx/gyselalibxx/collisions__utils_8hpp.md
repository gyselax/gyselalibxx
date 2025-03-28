

# File collisions\_utils.hpp



[**FileList**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**collisions\_utils.hpp**](collisions__utils_8hpp.md)

[Go to the source code of this file](collisions__utils_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "geometry.hpp"`
* `#include "quadrature.hpp"`
* `#include "trapezoid_quadrature.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**compute\_Dcoll**](#function-compute_dcoll) (DField&lt; IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), LocalGridVx &gt; &gt; Dcoll, DConstFieldSpX collfreq, DConstFieldSpX density, DConstFieldSpX temperature) <br>_Compute the intra species collision operator diffusion coefficient._  |
|  void | [**compute\_Nucoll**](#function-compute_nucoll) (DField&lt; IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), LocalGridVx &gt; &gt; Nucoll, DField&lt; IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), LocalGridVx &gt; &gt; Dcoll, DConstFieldSpX Vcoll, DConstFieldSpX Tcoll) <br>_Compute the intra species collision operator advection coefficient._  |
|  void | [**compute\_Vcoll\_Tcoll**](#function-compute_vcoll_tcoll) (DFieldSpX Vcoll, DFieldSpX Tcoll, DConstFieldSpXVx allfdistribu, DField&lt; IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), LocalGridVx &gt; &gt; Dcoll, DField&lt; IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), LocalGridVx &gt; &gt; dvDcoll) <br>_Compute the Vcoll and Tcoll coefficients, used for building the linear system._  |
|  void | [**compute\_collfreq**](#function-compute_collfreq) (DFieldSpX collfreq, DConstFieldSpX nustar\_profile, DConstFieldSpX density, DConstFieldSpX temperature) <br>_Compute the collision frequency for each species._  |
|  void | [**compute\_collfreq\_ab**](#function-compute_collfreq_ab) (DFieldSpX collfreq\_ab, DConstFieldSpX nustar\_profile, DConstFieldSpX density, DConstFieldSpX temperature) <br>_Compute the collision frequency between species a and b._  |
|  void | [**compute\_dvDcoll**](#function-compute_dvdcoll) (DField&lt; IdxRange&lt; [**Species**](structSpecies.md), [**GridX**](structGridX.md), LocalGridVx &gt; &gt; dvDcoll, DConstFieldSpX collfreq, DConstFieldSpX density, DConstFieldSpX temperature) <br>_Compute the velocity derivative of the collision operator diffusion coefficient._  |
|  void | [**compute\_momentum\_energy\_exchange**](#function-compute_momentum_energy_exchange) (DFieldSpX momentum\_exchange\_ab, DFieldSpX energy\_exchange\_ab, DConstFieldSpX collfreq\_ab, DConstFieldSpX density, DConstFieldSpX mean\_velocity, DConstFieldSpX temperature) <br>_Compute the momentum and energy exchange terms between species a and b._  |
|  void | [**compute\_nustar\_profile**](#function-compute_nustar_profile) (DFieldSpX nustar\_profile, double nustar0) <br>_Compute the collisionality spatial profile._  |




























## Public Functions Documentation




### function compute\_Dcoll 

_Compute the intra species collision operator diffusion coefficient._ 
```C++
template<class LocalGridVx>
void compute_Dcoll (
    DField< IdxRange< Species , GridX , LocalGridVx > > Dcoll,
    DConstFieldSpX collfreq,
    DConstFieldSpX density,
    DConstFieldSpX temperature
) 
```





**Parameters:**


* `Dcoll` A Field representing the diffusion coefficient. 
* `collfreq` The collision frequency for each species. 
* `density` The density of each species. 
* `temperature` The temperature of each species. 




        

<hr>



### function compute\_Nucoll 

_Compute the intra species collision operator advection coefficient._ 
```C++
template<class LocalGridVx>
void compute_Nucoll (
    DField< IdxRange< Species , GridX , LocalGridVx > > Nucoll,
    DField< IdxRange< Species , GridX , LocalGridVx > > Dcoll,
    DConstFieldSpX Vcoll,
    DConstFieldSpX Tcoll
) 
```





**Parameters:**


* `Nucoll` A Field representing the advection coefficient. 
* `Dcoll` A Field representing the diffusion coefficient. 
* `Vcoll` The Vcoll coefficient. 
* `Tcoll` The Tcoll coefficient. 




        

<hr>



### function compute\_Vcoll\_Tcoll 

_Compute the Vcoll and Tcoll coefficients, used for building the linear system._ 
```C++
template<class LocalGridVx>
void compute_Vcoll_Tcoll (
    DFieldSpX Vcoll,
    DFieldSpX Tcoll,
    DConstFieldSpXVx allfdistribu,
    DField< IdxRange< Species , GridX , LocalGridVx > > Dcoll,
    DField< IdxRange< Species , GridX , LocalGridVx > > dvDcoll
) 
```



Computation of Vcoll and Tcoll, which are the moments of the kernel maxwellian function of the intra species collision operator. Vcoll and Tcoll are defined as follows:
* 
* 
*  where the 5 integrals are defined as:   The brackets  represent the integral in velocity: 






**Parameters:**


* `Vcoll` The Vcoll coefficient. 
* `Tcoll` The Tcoll coefficient. 
* `allfdistribu` The distribution function. 
* `Dcoll` The collision operator diffusion coefficient. 
* `dvDcoll` The collision operator derivative of the diffusion coefficient. 




        

<hr>



### function compute\_collfreq 

_Compute the collision frequency for each species._ 
```C++
void compute_collfreq (
    DFieldSpX collfreq,
    DConstFieldSpX nustar_profile,
    DConstFieldSpX density,
    DConstFieldSpX temperature
) 
```





**Parameters:**


* `collfreq` A Field representing the collision frequency for each species. 
* `nustar_profile` The collisionality profile. 
* `density` The density of each species. 
* `temperature` The temperature of each species.

Computes the space and species dependent collision frequency collfreq. 


        

<hr>



### function compute\_collfreq\_ab 

_Compute the collision frequency between species a and b._ 
```C++
void compute_collfreq_ab (
    DFieldSpX collfreq_ab,
    DConstFieldSpX nustar_profile,
    DConstFieldSpX density,
    DConstFieldSpX temperature
) 
```





**Parameters:**


* `collfreq_ab` The collision frequency between species a and b. 
* `nustar_profile` The collisionality profile. 
* `density` The density of each species. 
* `temperature` The temperature of each species.

Computes the two species collision frequency collfreq\_ei 


        

<hr>



### function compute\_dvDcoll 

_Compute the velocity derivative of the collision operator diffusion coefficient._ 
```C++
template<class LocalGridVx>
void compute_dvDcoll (
    DField< IdxRange< Species , GridX , LocalGridVx > > dvDcoll,
    DConstFieldSpX collfreq,
    DConstFieldSpX density,
    DConstFieldSpX temperature
) 
```





**Parameters:**


* `dvDcoll` A Field representing the derivative of the diffusion coefficient. 
* `collfreq` The collision frequency for each species. 
* `density` The density of each species. 
* `temperature` The temperature of each species. 




        

<hr>



### function compute\_momentum\_energy\_exchange 

_Compute the momentum and energy exchange terms between species a and b._ 
```C++
void compute_momentum_energy_exchange (
    DFieldSpX momentum_exchange_ab,
    DFieldSpX energy_exchange_ab,
    DConstFieldSpX collfreq_ab,
    DConstFieldSpX density,
    DConstFieldSpX mean_velocity,
    DConstFieldSpX temperature
) 
```





**Parameters:**


* `momentum_exchange_ab` The momentum exchange term between species a and b. 
* `energy_exchange_ab` The energy exchange term between species a and b. 
* `collfreq_ab` The collision frequency between species a and b. 
* `density` The density of each species. 
* `mean_velocity` The mean velocity of each species. 
* `temperature` The temperature of each species.

Computes the momentum and energy exchange terms between ions and electrons 


        

<hr>



### function compute\_nustar\_profile 

_Compute the collisionality spatial profile._ 
```C++
void compute_nustar_profile (
    DFieldSpX nustar_profile,
    double nustar0
) 
```





**Parameters:**


* `nustar_profile` The collisionality profile. 
* `nustar0` normalised collisionality coefficient.

Computes the spatial profile of nustar, which is constant here, but could be space dependent (for instance to have no collisions in some specific parts of the simulation box). 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/collisions_utils.hpp`

