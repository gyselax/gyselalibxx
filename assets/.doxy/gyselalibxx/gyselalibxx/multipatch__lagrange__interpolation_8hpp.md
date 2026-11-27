

# File multipatch\_lagrange\_interpolation.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_54175eb71bf3bdb1ca864f2a0d352472.md) **>** [**multipatch\_lagrange\_interpolation.hpp**](multipatch__lagrange__interpolation_8hpp.md)

[Go to the source code of this file](multipatch__lagrange__interpolation_8hpp_source.md)



* `#include "ddc_aliases.hpp"`
* `#include "lagrange_interpolation.hpp"`
* `#include "multipatch_interpolation_2d.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef MultipatchInterpolation2D&lt; detail::LagrangeInterpolatorOnPatch&lt; typename PatchLocator::exec\_space, DataType, Basis1OnPatch, Basis2OnPatch &gt;::template type, PatchLocator, MultipatchExtrapolationRule &gt; | [**MultipatchLagrangeInterpolator2D**](#typedef-multipatchlagrangeinterpolator2d)  <br>_A MultipatchInterpolation2D which uses a 2D Lagrange interpolation on every patch._  |
















































## Public Types Documentation




### typedef MultipatchLagrangeInterpolator2D 

_A MultipatchInterpolation2D which uses a 2D Lagrange interpolation on every patch._ 
```C++
using MultipatchLagrangeInterpolator2D =  MultipatchInterpolation2D<
        detail::LagrangeInterpolatorOnPatch<
                typename PatchLocator::exec_space,
                DataType,
                Basis1OnPatch,
                Basis2OnPatch>::template type,
        PatchLocator,
        MultipatchExtrapolationRule>;
```



On each patch, the Lagrange interpolation is defined on the bases (Basis1OnPatch&lt;Patch&gt;, Basis2OnPatch&lt;Patch&gt;) and the grids (Grid1, Grid2) of the patch.


The extrapolation outside of the global domain is handled by MultipatchExtrapolationRule, so the extrapolation rules on each patch are null (or periodic along periodic dimensions).


The execution space is the one of the patch locator (PatchLocator::exec\_space).




**Template parameters:**


* `PatchLocator` The operator used to find the patch where a coordinate is physically located. 
* `MultipatchExtrapolationRule` The extrapolation rule used outside of the global domain (e.g. [**ConstantExtrapolationRuleOnion**](structConstantExtrapolationRuleOnion.md) or [**NullExtrapolationRule**](structNullExtrapolationRule.md)). 
* `Basis1OnPatch` A type template such that Basis1OnPatch&lt;Patch&gt; is the Lagrange basis ([**UniformLagrangeBasis**](classUniformLagrangeBasis.md) or [**NonUniformLagrangeBasis**](classNonUniformLagrangeBasis.md)) along the first dimension of the patch. 
* `Basis2OnPatch` A type template such that Basis2OnPatch&lt;Patch&gt; is the Lagrange basis along the second dimension of the patch. 
* `DataType` The type of the function values (default: double). 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interpolation/multipatch_lagrange_interpolation.hpp`

