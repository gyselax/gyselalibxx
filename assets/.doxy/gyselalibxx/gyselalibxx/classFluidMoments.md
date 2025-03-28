

# Class FluidMoments



[**ClassList**](annotated.md) **>** [**FluidMoments**](classFluidMoments.md)



_A class that computes fluid moments of the distribution function._ [More...](#detailed-description)

* `#include <fluid_moments.hpp>`















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**MomentDensity**](structFluidMoments_1_1MomentDensity.md) <br> |
| struct | [**MomentTemperature**](structFluidMoments_1_1MomentTemperature.md) <br> |
| struct | [**MomentVelocity**](structFluidMoments_1_1MomentVelocity.md) <br> |








## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr [**MomentDensity**](structFluidMoments_1_1MomentDensity.md) | [**s\_density**](#variable-s_density)   = `[**MomentDensity**](structFluidMoments_1_1MomentDensity.md)()`<br> |
|  constexpr [**MomentTemperature**](structFluidMoments_1_1MomentTemperature.md) | [**s\_temperature**](#variable-s_temperature)   = `[**MomentTemperature**](structFluidMoments_1_1MomentTemperature.md)()`<br> |
|  constexpr [**MomentVelocity**](structFluidMoments_1_1MomentVelocity.md) | [**s\_velocity**](#variable-s_velocity)   = `[**MomentVelocity**](structFluidMoments_1_1MomentVelocity.md)()`<br> |














## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**FluidMoments**](#function-fluidmoments) ([**Quadrature**](classQuadrature.md)&lt; IdxRangeVx, IdxRangeSpXVx &gt; integrate\_v) <br> |
|  void | [**operator()**](#function-operator) (double & density, DConstFieldVx fdistribu, [**MomentDensity**](structFluidMoments_1_1MomentDensity.md) moment\_density) <br> |
|  void | [**operator()**](#function-operator_1) (DFieldSpX density, DConstFieldSpXVx allfdistribu, [**MomentDensity**](structFluidMoments_1_1MomentDensity.md) moment\_density) <br> |
|  void | [**operator()**](#function-operator_2) (double & mean\_velocity, DConstFieldVx fdistribu, double density, [**MomentVelocity**](structFluidMoments_1_1MomentVelocity.md) moment\_velocity) <br> |
|  void | [**operator()**](#function-operator_3) (DFieldSpX mean\_velocity, DConstFieldSpXVx allfdistribu, DConstFieldSpX density, [**MomentVelocity**](structFluidMoments_1_1MomentVelocity.md) moment\_velocity) <br> |
|  void | [**operator()**](#function-operator_4) (double & temperature, DConstFieldVx fdistribu, double density, double mean\_velocity, [**MomentTemperature**](structFluidMoments_1_1MomentTemperature.md) moment\_temperature) <br> |
|  void | [**operator()**](#function-operator_5) (DFieldSpX temperature, DConstFieldSpXVx allfdistribu, DConstFieldSpX density, DConstFieldSpX mean\_velocity, [**MomentTemperature**](structFluidMoments_1_1MomentTemperature.md) moment\_temperature) <br> |
|   | [**~FluidMoments**](#function-fluidmoments) () = default<br> |




























## Detailed Description


These fluid moments are the density, mean velocity and temperature of the distribution function. 


    
## Public Static Attributes Documentation




### variable s\_density 

```C++
constexpr MomentDensity FluidMoments::s_density;
```



A static instance of [**MomentDensity**](structFluidMoments_1_1MomentDensity.md) that can be used to indicated to the operator() that the density should be calculated. 


        

<hr>



### variable s\_temperature 

```C++
constexpr MomentTemperature FluidMoments::s_temperature;
```



A static instance of [**MomentTemperature**](structFluidMoments_1_1MomentTemperature.md) that can be used to indicated to the operator() that the temperature should be calculated. 


        

<hr>



### variable s\_velocity 

```C++
constexpr MomentVelocity FluidMoments::s_velocity;
```



A static instance of [**MomentVelocity**](structFluidMoments_1_1MomentVelocity.md) that can be used to indicated to the operator() that the velocity should be calculated. 


        

<hr>
## Public Functions Documentation




### function FluidMoments 

```C++
FluidMoments::FluidMoments (
    Quadrature < IdxRangeVx, IdxRangeSpXVx > integrate_v
) 
```



The constructor for the operator.




**Parameters:**


* `integrate_v` A quadrature method which integrates over the velocity space. 




        

<hr>



### function operator() 

```C++
void FluidMoments::operator() (
    double & density,
    DConstFieldVx fdistribu,
    MomentDensity moment_density
) 
```



Calculate the density at a specific point of the distribution function.




**Parameters:**


* `density` The density at the point. 
* `fdistribu` A slice in velocity space of the distribution function at the given point. 
* `moment_density` A tag to ensure that the correct operator is called. 




        

<hr>



### function operator() 

```C++
void FluidMoments::operator() (
    DFieldSpX density,
    DConstFieldSpXVx allfdistribu,
    MomentDensity moment_density
) 
```



Calculate the density of the distribution function.




**Parameters:**


* `density` The density at various points for different species. 
* `allfdistribu` The distribution function. 
* `moment_density` A tag to ensure that the correct operator is called. 




        

<hr>



### function operator() 

```C++
void FluidMoments::operator() (
    double & mean_velocity,
    DConstFieldVx fdistribu,
    double density,
    MomentVelocity moment_velocity
) 
```



Calculate the mean velocity at a specific point of the distribution function.




**Parameters:**


* `mean_velocity` The mean velocity at the point. 
* `fdistribu` A slice in velocity space of the distribution function at the given point. 
* `density` The density at the point. 
* `moment_velocity` A tag to ensure that the correct operator is called. 




        

<hr>



### function operator() 

```C++
void FluidMoments::operator() (
    DFieldSpX mean_velocity,
    DConstFieldSpXVx allfdistribu,
    DConstFieldSpX density,
    MomentVelocity moment_velocity
) 
```



Calculate the mean velocity of the distribution function.




**Parameters:**


* `mean_velocity` The mean velocity at various points for different species. 
* `allfdistribu` The distribution function. 
* `density` The density at various points for different species. 
* `moment_velocity` A tag to ensure that the correct operator is called. 




        

<hr>



### function operator() 

```C++
void FluidMoments::operator() (
    double & temperature,
    DConstFieldVx fdistribu,
    double density,
    double mean_velocity,
    MomentTemperature moment_temperature
) 
```



Calculate the temperature at a specific point of the distribution function.




**Parameters:**


* `temperature` The mean temperature at the point. 
* `fdistribu` A slice in velocity space of the distribution function at the given point. 
* `density` The density at the point. 
* `mean_velocity` The mean velocity at the point. 
* `moment_temperature` A tag to ensure that the correct operator is called. 




        

<hr>



### function operator() 

```C++
void FluidMoments::operator() (
    DFieldSpX temperature,
    DConstFieldSpXVx allfdistribu,
    DConstFieldSpX density,
    DConstFieldSpX mean_velocity,
    MomentTemperature moment_temperature
) 
```



Calculate the mean temperature of the distribution function.




**Parameters:**


* `temperature` The mean temperature at various points for different species. 
* `allfdistribu` The distribution function. 
* `density` The density at various points for different species. 
* `mean_velocity` The mean velocity at various points for different species. 
* `moment_temperature` A tag to ensure that the correct operator is called. 




        

<hr>



### function ~FluidMoments 

```C++
FluidMoments::~FluidMoments () = default
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/utils/fluid_moments.hpp`

