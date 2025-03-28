

# Struct interpolator\_on\_idx\_range&lt; Interp, GridInterp, IdxRange&lt; Grid1D... &gt; &gt;

**template &lt;template&lt; class... &gt; class Interp, class GridInterp, class... Grid1D&gt;**



[**ClassList**](annotated.md) **>** [**interpolator\_on\_idx\_range&lt; Interp, GridInterp, IdxRange&lt; Grid1D... &gt; &gt;**](structinterpolator__on__idx__range_3_01Interp_00_01GridInterp_00_01IdxRange_3_01Grid1D_8_8_8_01_4_01_4.md)



[More...](#detailed-description)

* `#include <iinterpolator.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef Interp&lt; GridInterp, Grid1D... &gt; | [**type**](#typedef-type)  <br>_The type of the interpolator._  |
















































## Detailed Description


A structure which builds an interpolation type.




**Template parameters:**


* `Interp` The interpolator class being built. 
* `GridInterp` The dimension along which the operator will interpolate. 
* `Grid1D...` The dimensions on which the data being interpolated are defined. 




    
## Public Types Documentation




### typedef type 

_The type of the interpolator._ 
```C++
using interpolator_on_idx_range< Interp, GridInterp, IdxRange< Grid1D... > >::type =  Interp<GridInterp, Grid1D...>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/iinterpolator.hpp`

