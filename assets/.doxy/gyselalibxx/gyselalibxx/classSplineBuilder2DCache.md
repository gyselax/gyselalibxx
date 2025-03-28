

# Class SplineBuilder2DCache

**template &lt;class SplineBuilder2D&gt;**



[**ClassList**](annotated.md) **>** [**SplineBuilder2DCache**](classSplineBuilder2DCache.md)



_A class that stores spline builder coefficients and recomputes them when required._ 

* `#include <spline_builder_2d_cache.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SplineBuilder2DCache**](#function-splinebuilder2dcache) (SplineBuilder2D const & spline\_builder) <br>_Construct an instance of the class_ [_**SplineBuilder2DCache**_](classSplineBuilder2DCache.md) _._ |
|  DConstFieldSplineCoeffs | [**compute\_coeffs**](#function-compute_coeffs) (DConstField&lt; IdxRangeField &gt; field\_values) <br>_Compute the spline coefficients of the spline representation of field\_values if required, i.e. if not already computed by the other dimension than DimOfInterest._  |
|  DConstFieldSplineCoeffs | [**operator()**](#function-operator) () const<br>_Returns a constant field reference to the spline coefficients._  |




























## Public Functions Documentation




### function SplineBuilder2DCache 

_Construct an instance of the class_ [_**SplineBuilder2DCache**_](classSplineBuilder2DCache.md) _._
```C++
inline explicit SplineBuilder2DCache::SplineBuilder2DCache (
    SplineBuilder2D const & spline_builder
) 
```





**Parameters:**


* `spline_builder` A 2D spline builder. 




        

<hr>



### function compute\_coeffs 

_Compute the spline coefficients of the spline representation of field\_values if required, i.e. if not already computed by the other dimension than DimOfInterest._ 
```C++
template<class DimOfInterest>
inline DConstFieldSplineCoeffs SplineBuilder2DCache::compute_coeffs (
    DConstField< IdxRangeField > field_values
) 
```





**Parameters:**


* `field_values` The field to be used to compute spline coefficients.



**Returns:**

The spline coefficients updated if required. 





        

<hr>



### function operator() 

_Returns a constant field reference to the spline coefficients._ 
```C++
inline DConstFieldSplineCoeffs SplineBuilder2DCache::operator() () const
```





**Returns:**

A reference to a constant field that contains the spline coefficients. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/spline_builder_2d_cache.hpp`

