

# Namespace ExtrapolationRule



[**Namespace List**](namespaces.md) **>** [**ExtrapolationRule**](namespaceExtrapolationRule.md)



_A namespace containing tag types describing how a function is extrapolated outside the interpolation domain._ [More...](#detailed-description)
















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Constant**](structExtrapolationRule_1_1Constant.md) <br>_Tag selecting constant extrapolation._  |
| struct | [**NullValue**](structExtrapolationRule_1_1NullValue.md) <br>_Tag selecting null extrapolation: the function evaluates to zero outside the domain._  |
| struct | [**OneSidedPeriodic**](structExtrapolationRule_1_1OneSidedPeriodic.md) <br>_Tag selecting periodic extrapolation._  |


## Public Types

| Type | Name |
| ---: | :--- |
| typedef ddc::detail::TypeSeq&lt; [**Constant**](structExtrapolationRule_1_1Constant.md), [**Constant**](structExtrapolationRule_1_1Constant.md) &gt; | [**Constant\_Constant**](#typedef-constant_constant)  <br>_Convenience pairing of_ [_**ExtrapolationRule::Constant**_](structExtrapolationRule_1_1Constant.md) _for both boundaries._ |
| typedef ddc::detail::TypeSeq&lt; [**NullValue**](structExtrapolationRule_1_1NullValue.md), [**NullValue**](structExtrapolationRule_1_1NullValue.md) &gt; | [**Null\_Null**](#typedef-null_null)  <br>_Convenience pairing of_ [_**ExtrapolationRule::NullValue**_](structExtrapolationRule_1_1NullValue.md) _for both boundaries._ |
| typedef ddc::detail::TypeSeq&lt; [**OneSidedPeriodic**](structExtrapolationRule_1_1OneSidedPeriodic.md), [**OneSidedPeriodic**](structExtrapolationRule_1_1OneSidedPeriodic.md) &gt; | [**Periodic**](#typedef-periodic)  <br>_Convenience pairing of_ [_**ExtrapolationRule::OneSidedPeriodic**_](structExtrapolationRule_1_1OneSidedPeriodic.md) _for both boundaries._ |
















































## Detailed Description


Each tag exposes a type alias template that resolves, given the coefficient grid and data type of the interpolator, to the concrete extrapolation rule class to use. This keeps the set of extrapolation rules open for extension: user code may define its own tag following the same pattern, or bypass tags entirely and pass an already-concrete extrapolation rule class directly to SplineInterpolator or LagrangeInterpolator (see extrapolation\_rule\_choice.hpp for the generic resolution mechanism). 


    
## Public Types Documentation




### typedef Constant\_Constant 

_Convenience pairing of_ [_**ExtrapolationRule::Constant**_](structExtrapolationRule_1_1Constant.md) _for both boundaries._
```C++
using ExtrapolationRule::Constant_Constant = typedef ddc::detail::TypeSeq<Constant, Constant>;
```




<hr>



### typedef Null\_Null 

_Convenience pairing of_ [_**ExtrapolationRule::NullValue**_](structExtrapolationRule_1_1NullValue.md) _for both boundaries._
```C++
using ExtrapolationRule::Null_Null = typedef ddc::detail::TypeSeq<NullValue, NullValue>;
```




<hr>



### typedef Periodic 

_Convenience pairing of_ [_**ExtrapolationRule::OneSidedPeriodic**_](structExtrapolationRule_1_1OneSidedPeriodic.md) _for both boundaries._
```C++
using ExtrapolationRule::Periodic = typedef ddc::detail::TypeSeq<OneSidedPeriodic, OneSidedPeriodic>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation.hpp`

