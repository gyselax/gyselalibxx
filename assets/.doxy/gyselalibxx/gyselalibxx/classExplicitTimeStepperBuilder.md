

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


## Public Static Functions

| Type | Name |
| ---: | :--- |
|  auto | [**preallocate**](#function-preallocate-12) (typename ChosenTimeStepper::IdxRange const idx\_range) <br>_Allocate the TimeStepper object for FieldLike types._  |
|  auto | [**preallocate**](#function-preallocate-22) () <br>_Allocate the TimeStepper object for scalar (non-FieldLike) types._  |


























## Detailed Description


This class is a time stepper builder. A time stepper builder is designed to construct a time stepper upon request. This allows the simulation to choose the method without needing to know the specifics of the types with which it should be initialised. This class should be specialised for the explicit time stepper builders. 


    
## Public Types Documentation




### typedef time\_stepper\_t 

```C++
using ExplicitTimeStepperBuilder< TimeStepper >::time_stepper_t =  TimeStepper<FieldMem, DerivFieldMem, ExecSpace>;
```



The type of the TimeStepper that will be constructed to solve an equation whose field and derivative(s) have the specified type. 

**Template parameters:**


* `FieldMem` The type of the data storage for the function. 
* [**DerivFieldMem**](classDerivFieldMem.md) The type of the data storage for the derivative of the function. 
* `ExecSpace` The space (CPU/GPU) where the calculations are carried out. This template parameter is ignored if the FieldMem is a scalar. 




        

<hr>
## Public Functions Documentation




### function ExplicitTimeStepperBuilder 

_A constructor for the TimeStepperBuilder._ 
```C++
inline ExplicitTimeStepperBuilder::ExplicitTimeStepperBuilder () 
```




<hr>
## Public Static Functions Documentation




### function preallocate [1/2]

_Allocate the TimeStepper object for FieldLike types._ 
```C++
template<class ChosenTimeStepper>
static inline auto ExplicitTimeStepperBuilder::preallocate (
    typename ChosenTimeStepper::IdxRange const idx_range
) 
```





**Template parameters:**


* `ChosenTimeStepper` The type of the TimeStepper to be constructed (obtained from time\_stepper\_t). 



**Parameters:**


* `idx_range` The index range on which the operator will act (and allocate memory). 




        

<hr>



### function preallocate [2/2]

_Allocate the TimeStepper object for scalar (non-FieldLike) types._ 
```C++
template<class ChosenTimeStepper>
static inline auto ExplicitTimeStepperBuilder::preallocate () 
```





**Template parameters:**


* `ChosenTimeStepper` The type of the TimeStepper to be constructed (obtained from time\_stepper\_t). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/itimestepper.hpp`

