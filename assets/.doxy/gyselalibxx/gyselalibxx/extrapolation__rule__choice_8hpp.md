

# File extrapolation\_rule\_choice.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**extrapolation\_rule\_choice.hpp**](extrapolation__rule__choice_8hpp.md)

[Go to the source code of this file](extrapolation__rule__choice_8hpp_source.md)



* `#include "constant_identity_interpolation_extrapolation_rule.hpp"`
* `#include "i_interpolation.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef details::GetExtrapolationRuleClass&lt; Rule, CoeffGrid, DataType &gt;::type | [**extrapolation\_rule\_t**](#typedef-extrapolation_rule_t)  <br> |




















## Public Functions

| Type | Name |
| ---: | :--- |
|  extrapolation\_rule\_t&lt; Rule, CoeffGrid, DataType &gt; | [**get\_extrapolation**](#function-get_extrapolation) (Extremity extremity) <br> |




























## Public Types Documentation




### typedef extrapolation\_rule\_t 

```C++
using extrapolation_rule_t =  details::GetExtrapolationRuleClass<Rule, CoeffGrid, DataType>::type;
```




<hr>
## Public Functions Documentation




### function get\_extrapolation 

```C++
template<ExtrapolationRule Rule, class CoeffGrid, class DataType, class Basis>
extrapolation_rule_t< Rule, CoeffGrid, DataType > get_extrapolation (
    Extremity extremity
) 
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/extrapolation_rule_choice.hpp`

