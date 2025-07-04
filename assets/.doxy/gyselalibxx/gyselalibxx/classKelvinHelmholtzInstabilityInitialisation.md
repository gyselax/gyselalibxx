

# Class KelvinHelmholtzInstabilityInitialisation



[**ClassList**](annotated.md) **>** [**KelvinHelmholtzInstabilityInitialisation**](classKelvinHelmholtzInstabilityInitialisation.md)



_Initialise the allfdistribu function._ [More...](#detailed-description)

* `#include <initialisation_Kelvin_Helmholtz.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**KelvinHelmholtzInstabilityInitialisation**](#function-kelvinhelmholtzinstabilityinitialisation) (double const epsilon, double const mode\_k) <br>_Instantiate the initializer._  |
|  void | [**operator()**](#function-operator) (DFieldXY allfdistribu, DFieldXY allfdistribu\_equilibrium) <br>_Initialise_ \(f_{eq}\) _and_\(f\) _._ |
|   | [**~KelvinHelmholtzInstabilityInitialisation**](#function-kelvinhelmholtzinstabilityinitialisation) () = default<br> |




























## Detailed Description


Set the allfdistribu at 
\[f_{eq}(x, y) = \sin(y)\]



and \(f_0(x, y) = f_{eq}(x, y) + \varepsilon \cos(kx)\),


with \(\varepsilon\) an amplitude of perturbation and \(k\) mode equal to \(2\pi\) divided by the length of the index range on \(x\). 


    
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


* `epsilon` \(\varepsilon\), the amplitude of perturbation. 
* `mode_k` \(k\), the perturbation mode. 




        

<hr>



### function operator() 

_Initialise_ \(f_{eq}\) _and_\(f\) _._
```C++
inline void KelvinHelmholtzInstabilityInitialisation::operator() (
    DFieldXY allfdistribu,
    DFieldXY allfdistribu_equilibrium
) 
```





**Parameters:**


* `allfdistribu` Field referring to the \(f\) function. 
* `allfdistribu_equilibrium` Field referring to the \(f_{eq}\) function. 




        

<hr>



### function ~KelvinHelmholtzInstabilityInitialisation 

```C++
KelvinHelmholtzInstabilityInitialisation::~KelvinHelmholtzInstabilityInitialisation () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXY/initialisation/initialisation_Kelvin_Helmholtz.hpp`

