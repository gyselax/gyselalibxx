

# Class ExplicitTimeStepperBuilder

**template &lt;template&lt; class FieldMem, class [**DerivFieldMem**](classDerivFieldMem.md), class ExecSpace &gt; typename TimeStepper&gt;**



[**ClassList**](annotated.md) **>** [**ExplicitTimeStepperBuilder**](classExplicitTimeStepperBuilder.md)



_A class to indicate that an explicit time stepper should be constructed for use in other operators._ [More...](#detailed-description)

* `#include <itimestepper.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef TimeStepper&lt; FieldMem, [**DerivFieldMem**](classDerivFieldMem.md), ExecSpace &gt; | [**time\_stepper\_t**](#typedef-time_stepper_t)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ExplicitTimeStepperBuilder**](#function-explicittimestepperbuilder) () <br>_A constructor for the TimeStepperBuilder._  |
|  auto | [**preallocate**](#function-preallocate) (typename ChosenTimeStepper::IdxRange const idx\_range) const<br>_Allocate the TimeStepper object._  |




























## Detailed Description


This class is a time stepper builder. A time stepper builder is designed to construct a time stepper upon request. This allows the simulation to choose the method without needing to know the specifics of the types with which it should be initialised. This class should be specialised for the explicit time stepper builders. 


    
## Public Types Documentation




### typedef time\_stepper\_t 

```C++
using ExplicitTimeStepperBuilder< TimeStepper >::time_stepper_t =  TimeStepper<FieldMem, DerivFieldMem, ExecSpace>;
```



The type of the TimeStepper that will be constructed to solve an equation whose field and derivative(s) have the specified type. 


        

<hr>
## Public Functions Documentation




### function ExplicitTimeStepperBuilder 

_A constructor for the TimeStepperBuilder._ 
```C++
inline ExplicitTimeStepperBuilder::ExplicitTimeStepperBuilder () 
```




<hr>



### function preallocate 

_Allocate the TimeStepper object._ 
```C++
template<class ChosenTimeStepper>
inline auto ExplicitTimeStepperBuilder::preallocate (
    typename ChosenTimeStepper::IdxRange const idx_range
) const
```





**Template parameters:**


* `ChosenTimeStepper` The type of the TimeStepper to be constructed (obtained from time\_stepper\_t). 



**Parameters:**


* `idx_range` The index range on which the operator will act (and allocate memory). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/itimestepper.hpp`

