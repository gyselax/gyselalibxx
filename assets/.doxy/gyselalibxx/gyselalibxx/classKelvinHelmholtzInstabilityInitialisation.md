

# Class KelvinHelmholtzInstabilityInitialisation



[**ClassList**](annotated.md) **>** [**KelvinHelmholtzInstabilityInitialisation**](classKelvinHelmholtzInstabilityInitialisation.md)



_Initialise the allfdistribu function._ [More...](#detailed-description)

* `#include <initialisation_Kelvin_Helmholtz.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**KelvinHelmholtzInstabilityInitialisation**](#function-kelvinhelmholtzinstabilityinitialisation) (double const epsilon, double const mode\_k) <br>_Instantiate the initializer._  |
|  void | [**operator()**](#function-operator) (DFieldXY allfdistribu, DFieldXY allfdistribu\_equilibrium) <br>_Initialise_  _and_ _._ |
|   | [**~KelvinHelmholtzInstabilityInitialisation**](#function-kelvinhelmholtzinstabilityinitialisation) () = default<br> |




























## Detailed Description


Set the allfdistribu at 


and ,


with  an amplitude of perturbation and  mode equal to  divided by the length of the index range on . 


    
## Public Functions Documentation




### function KelvinHelmholtzInstabilityInitialisation 

_Instantiate the initializer._ 
```C++
inline KelvinHelmholtzInstabilityInitialisation::KelvinHelmholtzInstabilityInitialisation (
    double const epsilon,
    double const mode_k
) 
```





**Parameters:**


* `epsilon` , the amplitude of perturbation. 
* `mode_k` , the perturbation mode. 




        

<hr>



### function operator() 

_Initialise_  _and_ _._
```C++
inline void KelvinHelmholtzInstabilityInitialisation::operator() (
    DFieldXY allfdistribu,
    DFieldXY allfdistribu_equilibrium
) 
```





**Parameters:**


* `allfdistribu` Field referring to the  function. 
* `allfdistribu_equilibrium` Field referring to the  function. 




        

<hr>



### function ~KelvinHelmholtzInstabilityInitialisation 

```C++
KelvinHelmholtzInstabilityInitialisation::~KelvinHelmholtzInstabilityInitialisation () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXY/initialisation/initialisation_Kelvin_Helmholtz.hpp`

