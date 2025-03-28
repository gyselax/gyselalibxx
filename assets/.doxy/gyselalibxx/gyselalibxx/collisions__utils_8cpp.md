

# File collisions\_utils.cpp



[**FileList**](files.md) **>** [**geometryXVx**](dir_e51b496b46dd687775e46e0826614574.md) **>** [**rhs**](dir_53474cb30a3389ee74cb3186cae99ac0.md) **>** [**collisions\_utils.cpp**](collisions__utils_8cpp.md)

[Go to the source code of this file](collisions__utils_8cpp_source.md)



* `#include <cassert>`
* `#include <cmath>`
* `#include <optional>`
* `#include "collisions_utils.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "quadrature.hpp"`
* `#include "species_info.hpp"`
* `#include "trapezoid_quadrature.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**compute\_collfreq**](#function-compute_collfreq) (DFieldSpX collfreq, DConstFieldSpX nustar\_profile, DConstFieldSpX density, DConstFieldSpX temperature) <br>_Compute the collision frequency for each species._  |
|  void | [**compute\_collfreq\_ab**](#function-compute_collfreq_ab) (DFieldSpX collfreq\_ab, DConstFieldSpX nustar\_profile, DConstFieldSpX density, DConstFieldSpX temperature) <br>_Compute the collision frequency between species a and b._  |
|  void | [**compute\_momentum\_energy\_exchange**](#function-compute_momentum_energy_exchange) (DFieldSpX momentum\_exchange\_ab, DFieldSpX energy\_exchange\_ab, DConstFieldSpX collfreq\_ab, DConstFieldSpX density, DConstFieldSpX mean\_velocity, DConstFieldSpX temperature) <br>_Compute the momentum and energy exchange terms between species a and b._  |
|  void | [**compute\_nustar\_profile**](#function-compute_nustar_profile) (DFieldSpX nustar\_profile, double nustar0) <br>_Compute the collisionality spatial profile._  |




























## Public Functions Documentation




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



Computes the two species collision frequency collfreq\_ei 


        

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



Computes the spatial profile of nustar, which is constant here, but could be space dependent (for instance to have no collisions in some specific parts of the simulation box). 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/rhs/collisions_utils.cpp`

