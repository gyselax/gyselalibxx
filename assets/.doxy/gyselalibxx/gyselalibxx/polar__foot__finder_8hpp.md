

# File polar\_foot\_finder.hpp



[**FileList**](files.md) **>** [**advection**](dir_b90fde0f10c67a9aef841a6e6700f1f6.md) **>** [**polar\_foot\_finder.hpp**](polar__foot__finder_8hpp.md)

[Go to the source code of this file](polar__foot__finder_8hpp_source.md)



* `#include <source_location>`
* `#include <ddc/ddc.hpp>`
* `#include "polar_foot_finders/elementwise_choice.hpp"`
* `#include "polar_foot_finders/logical_advection_logical_foot_finder.hpp"`
* `#include "polar_foot_finders/logical_advection_pseudo_physical_foot_finder.hpp"`
* `#include "polar_foot_finders/physical_advection_physical_foot_finder.hpp"`
* `#include "polar_foot_finders/physical_advection_pseudo_physical_foot_finder.hpp"`
* `#include "ddc_alias_inline_functions.hpp"`
* `#include "ddc_aliases.hpp"`
* `#include "geometry_pseudo_cartesian.hpp"`
* `#include "i_interpolation.hpp"`
* `#include "i_interpolation_builder.hpp"`
* `#include "l_norm_tools.hpp"`
* `#include "vector_index_tools.hpp"`















## Classes

| Type | Name |
| ---: | :--- |
| class | [**PolarFootFinder**](classPolarFootFinder.md) &lt;FFSpace, AFSpace, LogicalToPhysicalMapping, class IdxRangeBatched, class TimeStepperBuilder, RThetaAdvectionInterpolator&gt;<br>_Operator for finding the feet of the characteristics on a polar slice._  |






















## Public Functions

| Type | Name |
| ---: | :--- |
|  auto | [**make\_polar\_foot\_finder**](#function-make_polar_foot_finder) (TimeStepperBuilder const & time\_stepper, LogicalToPhysicalMapping const & mapping, IdxRangeBatched const & idx\_range, RThetaAdvectionInterpolator const & interpolator\_advection\_field, Coord&lt; [**X\_pC**](structX__pC.md), [**Y\_pC**](structY__pC.md) &gt; coord\_centre=Coord&lt; [**X\_pC**](structX__pC.md), [**Y\_pC**](structY__pC.md) &gt;(0, 0), double epsilon=1e-12) <br>_Construct a_ [_**PolarFootFinder**_](classPolarFootFinder.md) _, deducing all type template parameters from the arguments._ |




























## Public Functions Documentation




### function make\_polar\_foot\_finder 

_Construct a_ [_**PolarFootFinder**_](classPolarFootFinder.md) _, deducing all type template parameters from the arguments._
```C++
template<FootFindingSpace FFSpace, AdvectionFieldSpace AFSpace, concepts::Mapping LogicalToPhysicalMapping, class IdxRangeBatched, class TimeStepperBuilder, concepts::Interpolation RThetaAdvectionInterpolator>
auto make_polar_foot_finder (
    TimeStepperBuilder const & time_stepper,
    LogicalToPhysicalMapping const & mapping,
    IdxRangeBatched const & idx_range,
    RThetaAdvectionInterpolator const & interpolator_advection_field,
    Coord< X_pC , Y_pC > coord_centre=Coord< X_pC , Y_pC >(0, 0),
    double epsilon=1e-12
) 
```



`FFSpace` and `AFSpace` must still be provided explicitly (they are enum values, not types). All remaining template parameters are deduced from the function arguments, eliminating the need for `decltype` at the call site.




**Parameters:**


* `time_stepper` A builder for the time integration method. 
* `mapping` The mapping from the logical domain to the physical domain. 
* `idx_range` The batched index range over which the operator works (used only for type deduction; its value is not forwarded to the constructor). 
* `interpolator_advection_field` An interpolator to build and evaluate an approximation of the advection field. 
* `coord_centre` The polar-centre coordinate in pseudo-Cartesian space. 
* `epsilon` Linearisation parameter near the O-point. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/advection/polar_foot_finder.hpp`

