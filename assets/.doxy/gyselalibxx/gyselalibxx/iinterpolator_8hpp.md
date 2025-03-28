

# File iinterpolator.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**iinterpolator.hpp**](iinterpolator_8hpp.md)

[Go to the source code of this file](iinterpolator_8hpp_source.md)



* `#include <memory>`
* `#include <ddc/ddc.hpp>`
* `#include <ddc/kernels/splines/deriv.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "ddc_helper.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**IInterpolator**](classIInterpolator.md) &lt;class GridInterp, Grid1D&gt;<br>_A class which provides an interpolating function._  |
| class | [**IPreallocatableInterpolator**](classIPreallocatableInterpolator.md) &lt;class GridInterp, Grid1D&gt;<br>_A class which provides access to an interpolating function which can be preallocated where useful._  |
| struct | [**interpolator\_on\_idx\_range**](structinterpolator__on__idx__range.md) &lt;Interp, class GridInterp, class IdxRange&gt;<br> |
| struct | [**interpolator\_on\_idx\_range&lt; Interp, GridInterp, IdxRange&lt; Grid1D... &gt; &gt;**](structinterpolator__on__idx__range_3_01Interp_00_01GridInterp_00_01IdxRange_3_01Grid1D_8_8_8_01_4_01_4.md) &lt;Interp, class GridInterp, Grid1D&gt;<br> |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef typename [**interpolator\_on\_idx\_range**](structinterpolator__on__idx__range.md)&lt; Interp, GridInterp, IdxRange &gt;::type | [**interpolator\_on\_idx\_range\_t**](#typedef-interpolator_on_idx_range_t)  <br> |
















































## Public Types Documentation




### typedef interpolator\_on\_idx\_range\_t 

```C++
using interpolator_on_idx_range_t =  typename interpolator_on_idx_range<Interp, GridInterp, IdxRange>::type;
```



A template function which returns an interpolation type.




**Template parameters:**


* `Interp` The interpolator class being built. 
* `GridInterp` The dimension along which the operator will interpolate. 
* `IdxRange` The index range on which the data being interpolated is defined. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/iinterpolator.hpp`

