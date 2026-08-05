

# Struct SplineBoundaryClosures

**template &lt;ddc::SplineBuilderClosure MinClosure, ddc::SplineBuilderClosure MaxClosure&gt;**



[**ClassList**](annotated.md) **>** [**SplineBoundaryClosures**](structSplineBoundaryClosures.md)



_Groups the lower (min) and upper (max) ddc::SplineBuilderClosure of a spline builder into a single non-type template argument._ 

* `#include <spline_interpolation.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef std::integral\_constant&lt; ddc::SplineBuilderClosure, MaxClosure &gt; | [**max**](#typedef-max)  <br>_The ddc::SplineBuilderClosure at the upper boundary of the spline builder._  |
| typedef std::integral\_constant&lt; ddc::SplineBuilderClosure, MinClosure &gt; | [**min**](#typedef-min)  <br>_The ddc::SplineBuilderClosure at the lower boundary of the spline builder._  |
















































## Public Types Documentation




### typedef max 

_The ddc::SplineBuilderClosure at the upper boundary of the spline builder._ 
```C++
using SplineBoundaryClosures< MinClosure, MaxClosure >::max =  std::integral_constant<ddc::SplineBuilderClosure, MaxClosure>;
```




<hr>



### typedef min 

_The ddc::SplineBuilderClosure at the lower boundary of the spline builder._ 
```C++
using SplineBoundaryClosures< MinClosure, MaxClosure >::min =  std::integral_constant<ddc::SplineBuilderClosure, MinClosure>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/interpolation/spline_interpolation.hpp`

