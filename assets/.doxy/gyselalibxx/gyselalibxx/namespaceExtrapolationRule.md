

# Namespace ExtrapolationRule



[**Namespace List**](namespaces.md) **>** [**ExtrapolationRule**](namespaceExtrapolationRule.md)



_A namespace containing tag types describing how a function is extrapolated outside the interpolation domain._ [More...](#detailed-description)
















## Classes

| Type | Name |
| ---: | :--- |
| struct | [**Constant**](structExtrapolationRule_1_1Constant.md) <br>_Tag selecting constant extrapolation._  |
| struct | [**NullValue**](structExtrapolationRule_1_1NullValue.md) <br>_Tag selecting null extrapolation: the function evaluates to zero outside the domain._  |
| struct | [**Periodic**](structExtrapolationRule_1_1Periodic.md) <br>_Tag selecting periodic extrapolation._  |


















































## Detailed Description


Each tag exposes a type alias template that resolves, given the coefficient grid and data type of the interpolator, to the concrete extrapolation rule class to use. This keeps the set of extrapolation rules open for extension: user code may define its own tag following the same pattern, or bypass tags entirely and pass an already-concrete extrapolation rule class directly to [**SplineInterpolator**](classSplineInterpolator.md) or [**LagrangeInterpolator**](classLagrangeInterpolator.md) (see extrapolation\_rule\_choice.hpp for the generic resolution mechanism). 


    

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/i_interpolation.hpp`

