

# Class ITimeSolverRTheta



[**ClassList**](annotated.md) **>** [**ITimeSolverRTheta**](classITimeSolverRTheta.md)



_Base class for the time solvers._ 

* `#include <itimesolver.hpp>`





Inherited by the following classes: [BslExplicitPredCorrRTheta](classBslExplicitPredCorrRTheta.md),  [BslImplicitPredCorrRTheta](classBslImplicitPredCorrRTheta.md),  [BslPredCorrRTheta](classBslPredCorrRTheta.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual host\_t&lt; DFieldRTheta &gt; | [**operator()**](#function-operator) (host\_t&lt; DFieldRTheta &gt; allfdistribu, double const dt, int const steps=1) const = 0<br>_Solves on_  _the equations system._ |
| virtual  | [**~ITimeSolverRTheta**](#function-itimesolverrtheta) () = default<br> |
























## Protected Functions

| Type | Name |
| ---: | :--- |
|  void | [**display\_time\_difference**](#function-display_time_difference) (std::string const & title, std::chrono::time\_point&lt; std::chrono::system\_clock &gt; const & start\_time, std::chrono::time\_point&lt; std::chrono::system\_clock &gt; const & end\_time) const<br>_Displays the time difference between two given times and a title._  |




## Public Functions Documentation




### function operator() 

_Solves on_  _the equations system._
```C++
virtual host_t< DFieldRTheta > ITimeSolverRTheta::operator() (
    host_t< DFieldRTheta > allfdistribu,
    double const dt,
    int const steps=1
) const = 0
```





**Parameters:**


* `allfdistribu` On input: the initial condition. On output: the solution at . 
* `dt` The time step. 
* `steps` The number  of time interactions.



**Returns:**

A Field toward allfdistribu. 





        

<hr>



### function ~ITimeSolverRTheta 

```C++
virtual ITimeSolverRTheta::~ITimeSolverRTheta () = default
```




<hr>
## Protected Functions Documentation




### function display\_time\_difference 

_Displays the time difference between two given times and a title._ 
```C++
inline void ITimeSolverRTheta::display_time_difference (
    std::string const & title,
    std::chrono::time_point< std::chrono::system_clock > const & start_time,
    std::chrono::time_point< std::chrono::system_clock > const & end_time
) const
```



Useful to display the duration of a simulation.




**Parameters:**


* `title` A string printed in front of the computed duration. 
* `start_time` The start time point. 
* `end_time` The end time point. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryRTheta/time_solver/itimesolver.hpp`

