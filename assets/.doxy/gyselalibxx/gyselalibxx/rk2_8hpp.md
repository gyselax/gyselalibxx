

# File rk2.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**timestepper**](dir_ddbbe171637b3a2a6c78c931b02a7373.md) **>** [**rk2.hpp**](rk2_8hpp.md)

[Go to the source code of this file](rk2_8hpp_source.md)



* `#include <array>`
* `#include <type_traits>`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "itimestepper.hpp"`
* `#include "vector_field_common.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**RK2**](classRK2.md) &lt;class FieldMem, class [**DerivFieldMem**](classITimeStepper.md#typedef-derivfieldmem), class ExecSpace&gt;<br>_A class which provides an implementation of a second-order Runge-Kutta method._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ExplicitTimeStepperBuilder**](classExplicitTimeStepperBuilder.md)&lt; [**RK2**](classRK2.md) &gt; | [**RK2Builder**](#typedef-rk2builder)  <br> |
















































## Public Types Documentation




### typedef RK2Builder 

```C++
using RK2Builder =  ExplicitTimeStepperBuilder<RK2>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/rk2.hpp`

