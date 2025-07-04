

# Class DiocotronDensitySolution



[**ClassList**](annotated.md) **>** [**DiocotronDensitySolution**](classDiocotronDensitySolution.md)



_The diocotron exact solution of the density_ \(\rho\) _._[More...](#detailed-description)

* `#include <diocotron_initialisation_equilibrium.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**DiocotronDensitySolution**](#function-diocotrondensitysolution) (double const W1, double const R1, double const R2, double const W2, double const Q, int const l, double const eps) <br>_Instantiate a_ [_**DiocotronDensitySolution**_](classDiocotronDensitySolution.md) _._ |
|  double | [**equilibrium**](#function-equilibrium) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Get the equilibrium of the density._  |
|  double | [**get\_frequency**](#function-get_frequency) () const<br>_Get the frequency of the perturbation._  |
|  double | [**get\_slope**](#function-get_slope) () const<br>_Get the slope of the perturbation._  |
|  double | [**initialisation**](#function-initialisation) (Coord&lt; [**R**](structR.md), [**Theta**](structTheta.md) &gt; const & coord) const<br>_Get the initial condition of the density._  |




























## Detailed Description


The equations of the diocotron simulation are the Vlasov-Poisson equations



* \(\partial_t \rho - E_y \partial_x \rho + E_x \partial_y \rho = 0\),
* \(L \phi = - \nabla \cdot (\alpha \nabla \phi) + \beta \phi = \rho\),
* \(E = -\nabla \phi\).




The [**DiocotronDensitySolution**](classDiocotronDensitySolution.md) provides an initial perturbed solution \(\rho_0\) for the density \(\rho\) and its associated equilibrium solution \(\rho_{eq}\). 


    
## Public Functions Documentation




### function DiocotronDensitySolution 

_Instantiate a_ [_**DiocotronDensitySolution**_](classDiocotronDensitySolution.md) _._
```C++
inline DiocotronDensitySolution::DiocotronDensitySolution (
    double const W1,
    double const R1,
    double const R2,
    double const W2,
    double const Q,
    int const l,
    double const eps
) 
```





**Parameters:**


* `W1` The value of the r-coordinate where the domain starts. 
* `R1` The value of the r-coordinate where the domain of non null initial condition starts. 
* `R2` The value of the r-coordinate where the domain of non null initial condition ends. 
* `W2` The value of the r-coordinate where the domain ends. 
* `Q` The charge carried by the inner conductor at \(r = W_1\). 
* `l` The mode of the perturbation. 
* `eps` The amplitude of the perturbation. 




        

<hr>



### function equilibrium 

_Get the equilibrium of the density._ 
```C++
double DiocotronDensitySolution::equilibrium (
    Coord< R , Theta > const & coord
) const
```



The equilibrium is given by



* if \(R_1 \leq r \leq R_2\), \(\rho_0(r, \theta) = \exp\left( - \left(\frac{r - \bar{r}}{d}\right)^p \right)\),
* otherwise, \(\rho_0(r, \theta) = 0\),




with \(p = 50\), \(\bar{r} = \frac{R_1 + R_2}{2}\) and \(d = \frac{R_2 - R_1}{2}\).




**Parameters:**


* `coord` The coordinate where we evaluate the equilibrium.



**Returns:**

The value of the equilibrium at the given coordinate. 





        

<hr>



### function get\_frequency 

_Get the frequency of the perturbation._ 
```C++
double DiocotronDensitySolution::get_frequency () const
```



The frequency of the perturbation is given by the real part of \(\omega\), where \(\omega\) is the solution of the dispersion relation.




**Returns:**

The frequency of the perturbation. 





        

<hr>



### function get\_slope 

_Get the slope of the perturbation._ 
```C++
double DiocotronDensitySolution::get_slope () const
```



The slope of the perturbation is given by the imaginary part of \(\omega\), where \(\omega\) is the solution of the dispersion relation.




**Returns:**

The slope of the perturbation. 





        

<hr>



### function initialisation 

_Get the initial condition of the density._ 
```C++
double DiocotronDensitySolution::initialisation (
    Coord< R , Theta > const & coord
) const
```



The initial condition is given by



* if \(R_1 \leq r \leq R_2\), \(\rho_0(r, \theta) = (1 + \varepsilon \cos(l\theta)) \exp\left( - \left(\frac{r - \bar{r}}{d}\right)^p \right)\),
* otherwise, \(\rho_0(r, \theta) = 0\),




with \(p = 50\), \(\bar{r} = \frac{R_1 + R_2}{2}\) and \(d = \frac{R_2 - R_1}{2}\).




**Parameters:**


* `coord` The coordinate where we evaluate the initial condition.



**Returns:**

The value of the initial condition at the given coordinate. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/initialisation/diocotron_initialisation_equilibrium.hpp`

