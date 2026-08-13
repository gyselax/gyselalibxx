

# File extrapolation\_rule\_choice.hpp



[**FileList**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**extrapolation\_rule\_choice.hpp**](extrapolation__rule__choice_8hpp.md)

[Go to the source code of this file](extrapolation__rule__choice_8hpp_source.md)



* `#include "constant_identity_interpolation_extrapolation_rule.hpp"`
* `#include "i_interpolation.hpp"`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::detail::TypeSeq&lt; typename details::ExtrapolationRuleResolver&lt; ddc::type\_seq\_element\_t&lt; 0, Rule &gt;, DataType, Basis, IdxRangeCoeff &gt;::type, typename details::ExtrapolationRuleResolver&lt; ddc::type\_seq\_element\_t&lt; 1, Rule &gt;, DataType, Basis, IdxRangeCoeff &gt;::type &gt; | [**extrapolation\_rule\_t**](#typedef-extrapolation_rule_t)  <br>_Resolve Rule (either a self-describing tag or a concrete extrapolation rule class) to the concrete extrapolation rule class to use for a given Basis/DataType._  |




## Public Attributes

| Type | Name |
| ---: | :--- |
|  constexpr bool | [**is\_extrapolation\_rule\_auto\_constructible\_v**](#variable-is_extrapolation_rule_auto_constructible_v)   = `/* multi line expression */`<br>_True if get\_extrapolation&lt;Rule, CoeffGrid, DataType, ...&gt; can build the rule with no extra information: the resolved rule type is default-constructible, or Rule is the_ [_**ExtrapolationRule::Constant**_](structExtrapolationRule_1_1Constant.md) _tag (which get\_extrapolation knows how to build from the boundary of the discrete space)._ |
















## Public Functions

| Type | Name |
| ---: | :--- |
|  Rule | [**get\_extrapolation**](#function-get_extrapolation) (Extremity extremity) <br>_Initialise the extrapolation rule._  |




























## Public Types Documentation




### typedef extrapolation\_rule\_t 

_Resolve Rule (either a self-describing tag or a concrete extrapolation rule class) to the concrete extrapolation rule class to use for a given Basis/DataType._ 
```C++
using extrapolation_rule_t =  ddc::detail::TypeSeq< typename details::ExtrapolationRuleResolver< ddc::type_seq_element_t<0, Rule>, DataType, Basis, IdxRangeCoeff>::type, typename details::ExtrapolationRuleResolver< ddc::type_seq_element_t<1, Rule>, DataType, Basis, IdxRangeCoeff>::type>;
```




<hr>
## Public Attributes Documentation




### variable is\_extrapolation\_rule\_auto\_constructible\_v 

_True if get\_extrapolation&lt;Rule, CoeffGrid, DataType, ...&gt; can build the rule with no extra information: the resolved rule type is default-constructible, or Rule is the_ [_**ExtrapolationRule::Constant**_](structExtrapolationRule_1_1Constant.md) _tag (which get\_extrapolation knows how to build from the boundary of the discrete space)._
```C++
constexpr bool is_extrapolation_rule_auto_constructible_v;
```




<hr>
## Public Functions Documentation




### function get\_extrapolation 

_Initialise the extrapolation rule._ 
```C++
template<class Rule, class DataType, class CDim, class IdxRangeCoeff, class IdxRangeBasis>
Rule get_extrapolation (
    Extremity extremity
) 
```



Initialise the extrapolation rule using the default constructor if available. For a constant extrapolation, initialise the extrapolation from the selected extremity. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/extrapolation_rule_choice.hpp`

