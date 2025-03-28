

# File polar\_bsplines.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**polar\_splines**](dir_a6779ae02b71d57f488d261458bab1ce.md) **>** [**polar\_bsplines.hpp**](polar__bsplines_8hpp.md)

[Go to the source code of this file](polar__bsplines_8hpp_source.md)



* `#include <array>`
* `#include <vector>`
* `#include <ddc/ddc.hpp>`
* `#include "bernstein.hpp"`
* `#include "cartesian_to_barycentric.hpp"`
* `#include "ddc_helper.hpp"`
* `#include "discrete_to_cartesian.hpp"`
* `#include "mapping_tools.hpp"`
* `#include "polar_spline.hpp"`
* `#include "view.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**PolarBSplines**](classPolarBSplines.md) &lt;class [**BSplinesR**](structBSplinesR.md), class [**BSplinesTheta**](structBSplinesTheta.md), C&gt;<br> |
| class | [**Impl**](classPolarBSplines_1_1Impl.md) &lt;class DDim, class MemorySpace&gt;<br> |
| struct | [**Corner1Tag**](structPolarBSplines_1_1Impl_1_1Corner1Tag.md) <br>_The tag for the first corner of the Barycentric coordinates._  |
| struct | [**Corner2Tag**](structPolarBSplines_1_1Impl_1_1Corner2Tag.md) <br>_The tag for the second corner of the Barycentric coordinates._  |
| struct | [**Corner3Tag**](structPolarBSplines_1_1Impl_1_1Corner3Tag.md) <br>_The tag for the third corner of the Barycentric coordinates._  |
| struct | [**IntermediateBernsteinBasis**](structPolarBSplines_1_1Impl_1_1IntermediateBernsteinBasis.md) &lt;class DiscreteMapping&gt;<br> |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  [**PolarSpline**](structPolarSpline.md)&lt; DDim, MemorySpace &gt; | [**integrals**](#function-integrals) (ExecSpace const & execution\_space, [**PolarSpline**](structPolarSpline.md)&lt; DDim, MemorySpace &gt; int\_vals) <br> |




























## Public Functions Documentation




### function integrals 

```C++
template<class ExecSpace, class DDim, class MemorySpace>
PolarSpline < DDim, MemorySpace > integrals (
    ExecSpace const & execution_space,
    PolarSpline < DDim, MemorySpace > int_vals
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/polar_splines/polar_bsplines.hpp`

