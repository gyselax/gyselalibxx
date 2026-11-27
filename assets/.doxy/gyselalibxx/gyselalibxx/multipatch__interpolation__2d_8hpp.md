

# File multipatch\_interpolation\_2d.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_54175eb71bf3bdb1ca864f2a0d352472.md) **>** [**multipatch\_interpolation\_2d.hpp**](multipatch__interpolation__2d_8hpp.md)

[Go to the source code of this file](multipatch__interpolation__2d_8hpp_source.md)



* `#include <tuple>`
* `#include <type_traits>`
* `#include <ddc/ddc.hpp>`
* `#include "ddc_aliases.hpp"`
* `#include "i_interpolation.hpp"`
* `#include "interpolation_builder_2d_collection.hpp"`
* `#include "interpolation_evaluator_2d_collection.hpp"`
* `#include "multipatch_type.hpp"`
* `#include "type_seq_tools.hpp"`
* `#include "types.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef detail::MultipatchInterpolation2DImpl&lt; typename PatchLocator::PatchOrdering, type\_seq\_transform\_t&lt; InterpOnPatch, typename PatchLocator::PatchOrdering &gt;, PatchLocator, MultipatchExtrapolationRule, MultipatchBuilder &gt; | [**MultipatchInterpolation2D**](#typedef-multipatchinterpolation2d)  <br>_An owning object that bundles the 2D interpolations of every patch of a multipatch geometry._  |
















































## Public Types Documentation




### typedef MultipatchInterpolation2D 

_An owning object that bundles the 2D interpolations of every patch of a multipatch geometry._ 
```C++
using MultipatchInterpolation2D =  detail::MultipatchInterpolation2DImpl<
        typename PatchLocator::PatchOrdering,
        type_seq_transform_t<InterpOnPatch, typename PatchLocator::PatchOrdering>,
        PatchLocator,
        MultipatchExtrapolationRule,
        MultipatchBuilder>;
```



MultipatchInterpolation2D constructs and owns one 2D interpolation (see concepts::Interpolation, e.g. SplineInterpolator or LagrangeInterpolator) per patch, and exposes them through an [**InterpolationBuilder2DCollection**](classInterpolationBuilder2DCollection.md) (get\_builder()) and an [**InterpolationEvaluator2DCollection**](classInterpolationEvaluator2DCollection.md) (get\_evaluator()). It is the multipatch equivalent of the Interpolation classes and avoids building the builders, evaluators and extrapolation rules of each patch manually.


The patches are taken from the patch locator (PatchLocator::PatchOrdering).


By default the builder (get\_builder()) is the collection of the independent builders of each patch (get\_local\_builder()). A different MultipatchBuilder (e.g. one which also estimates the derivatives at the patch interfaces) can be provided. It is constructed from the collection and any additional arguments passed to the constructor.


The builders refer to the builders owned by this object, so this object must outlive any copy of the builders and can be neither copied nor moved.




**Template parameters:**


* `InterpOnPatch` A type template such that InterpOnPatch&lt;Patch&gt; is the 2D interpolation used on the given patch. Its interpolation grids must be the grids (Grid1, Grid2) of the patch. 
* `PatchLocator` The operator used to find the patch where a coordinate is physically located. 
* `MultipatchExtrapolationRule` The extrapolation rule used outside of the global domain (e.g. [**ConstantExtrapolationRuleOnion**](structConstantExtrapolationRuleOnion.md) or [**NullExtrapolationRule**](structNullExtrapolationRule.md)). 
* `MultipatchBuilder` The multipatch builder (default: the [**InterpolationBuilder2DCollection**](classInterpolationBuilder2DCollection.md) of the builders of each patch). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interpolation/multipatch_interpolation_2d.hpp`

