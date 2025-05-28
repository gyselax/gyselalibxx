

# File itimestepper.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**timestepper**](dir_ddbbe171637b3a2a6c78c931b02a7373.md) **>** [**itimestepper.hpp**](itimestepper_8hpp.md)

[Go to the source code of this file](itimestepper_8hpp_source.md)



* `#include <array>`
* `#include <type_traits>`
* `#include <ddc/ddc.hpp>`
* `#include "multipatch_field.hpp"`
* `#include "multipatch_field_mem.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**ExplicitTimeStepperBuilder**](classExplicitTimeStepperBuilder.md) &lt;TimeStepper&gt;<br>_A class to indicate that an explicit time stepper should be constructed for use in other operators._  |
| class | [**ITimeStepper**](classITimeStepper.md) &lt;class FieldMem, class DerivFieldMemType, class ExecSpace&gt;<br>_The superclass from which all timestepping methods inherit._  |






## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**is\_timestepper\_builder\_v**](#variable-is_timestepper_builder_v)   = `detail::enable\_is\_timestepper\_builder&lt;std::remove\_const\_t&lt;std::remove\_reference\_t&lt;Type&gt;&gt;&gt;`<br> |












































## Public Attributes Documentation




### variable is\_timestepper\_builder\_v 

```C++
constexpr bool is_timestepper_builder_v;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/timestepper/itimestepper.hpp`

