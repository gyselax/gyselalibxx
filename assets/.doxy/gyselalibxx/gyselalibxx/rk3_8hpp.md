

# File rk3.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**timestepper**](dir_ddbbe171637b3a2a6c78c931b02a7373.md) **>** [**rk3.hpp**](rk3_8hpp.md)

[Go to the source code of this file](rk3_8hpp_source.md)



* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "itimestepper.hpp"`
* `#include "vector_field_common.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**RK3**](classRK3.md) &lt;class ValType, class DerivType, class ExecSpace&gt;<br>_A class which provides an implementation of a third-order Runge-Kutta method._  |
| class | [**RK3&lt; FieldMem, DerivFieldMem, ExecSpace &gt;**](classRK3_3_01FieldMem_00_01DerivFieldMem_00_01ExecSpace_01_4.md) &lt;FieldMem, DerivFieldMem, class ExecSpace&gt;<br>_A class which provides an implementation of a third-order Runge-Kutta method._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef [**ExplicitTimeStepperBuilder**](classExplicitTimeStepperBuilder.md)&lt; [**RK3**](classRK3.md) &gt; | [**RK3Builder**](#typedef-rk3builder)  <br> |
















































## Public Types Documentation




### typedef RK3Builder 

```C++
using RK3Builder =  ExplicitTimeStepperBuilder<RK3>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/rk3.hpp`

