

# Class CrankNicolsonBuilder



[**ClassList**](annotated.md) **>** [**CrankNicolsonBuilder**](classCrankNicolsonBuilder.md)



_A class to indicate that a Crank-Nicolson time stepper should be constructed for use in other operators._ [More...](#detailed-description)

* `#include <crank_nicolson.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**CrankNicolson**](classCrankNicolson.md)&lt; FieldMem, [**DerivFieldMem**](classDerivFieldMem.md), ExecSpace &gt; | [**time\_stepper\_t**](#typedef-time_stepper_t)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CrankNicolsonBuilder**](#function-cranknicolsonbuilder) (int const counter=20, double const epsilon=1e-12) <br>_A constructor for the_ [_**CrankNicolsonBuilder**_](classCrankNicolsonBuilder.md) _._ |
|  auto | [**preallocate**](#function-preallocate) (typename TimeStepper::IdxRange const idx\_range) const<br>_Allocate the TimeStepper object._  |




























## Detailed Description


This class is a time stepper builder. A time stepper builder is designed to construct a time stepper upon request. This allows the simulation to choose the method without needing to know the specifics of the types with which it should be initialised. 


    
## Public Types Documentation




### typedef time\_stepper\_t 

```C++
using CrankNicolsonBuilder::time_stepper_t =  CrankNicolson<FieldMem, DerivFieldMem, ExecSpace>;
```



The type of the TimeStepper that will be constructed to solve an equation whose field and derivative(s) have the specified type. 


        

<hr>
## Public Functions Documentation




### function CrankNicolsonBuilder 

_A constructor for the_ [_**CrankNicolsonBuilder**_](classCrankNicolsonBuilder.md) _._
```C++
inline explicit CrankNicolsonBuilder::CrankNicolsonBuilder (
    int const counter=20,
    double const epsilon=1e-12
) 
```





**Parameters:**


* `counter` The maximal number of loops for the implicit method. 
* `epsilon` The  upperbound of the difference of two steps in the implicit method: . 




        

<hr>



### function preallocate 

_Allocate the TimeStepper object._ 
```C++
template<class TimeStepper>
inline auto CrankNicolsonBuilder::preallocate (
    typename TimeStepper::IdxRange const idx_range
) const
```





**Template parameters:**


* `ChosenTimeStepper` The type of the TimeStepper to be constructed (obtained from time\_stepper\_t). 



**Parameters:**


* `idx_range` The index range on which the operator will act (and allocate memory). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/crank_nicolson.hpp`

