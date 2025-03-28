

# Class Spline2DPartialDerivativeCreator

**template &lt;class [**SplineBuilder2DCache**](classSplineBuilder2DCache.md), class SplineEvaluator2D, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**Spline2DPartialDerivativeCreator**](classSpline2DPartialDerivativeCreator.md)



_A class which stores information necessary to create a pointer to an instance of the_ [_**Spline2DPartialDerivative**_](classSpline2DPartialDerivative.md) _class._[More...](#detailed-description)

* `#include <spline_2d_partial_derivative.hpp>`



Inherits the following classes: [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Spline2DPartialDerivativeCreator**](#function-spline2dpartialderivativecreator) ([**SplineBuilder2DCache**](classSplineBuilder2DCache.md) & builder\_cache, SplineEvaluator2D const & evaluator) <br>_Construct an instance of the_ [_**Spline2DPartialDerivativeCreator**_](classSpline2DPartialDerivativeCreator.md) _class._ |
|  std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; typename SplineEvaluator2D::batched\_evaluation\_domain\_type, DerivativeDimension &gt; &gt; | [**create\_instance**](#function-create_instance) (DConstFieldType field) const<br> |


## Public Functions inherited from IPartialDerivativeCreator

See [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)

| Type | Name |
| ---: | :--- |
| virtual std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](classIPartialDerivativeCreator.md#function-create_instance) (DConstField&lt; IdxRangeFull &gt; field) const = 0<br>_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._ |






















































## Detailed Description


This class allows an instance of the [**Spline2DPartialDerivative**](classSpline2DPartialDerivative.md) class to be instantiated where necessary. Typically, the [**Spline2DPartialDerivativeCreator**](classSpline2DPartialDerivativeCreator.md) is instantiated in the initialisation of the simulation, and the corresponding [**Spline2DPartialDerivative**](classSpline2DPartialDerivative.md) object is instantiated where computing partial derivatives is required.




**Template parameters:**


* `SplineBuilder2D` A 2D spline builder. 
* `SplineEvaluator2D` A 2D spline evaluator. 




    
## Public Functions Documentation




### function Spline2DPartialDerivativeCreator 

_Construct an instance of the_ [_**Spline2DPartialDerivativeCreator**_](classSpline2DPartialDerivativeCreator.md) _class._
```C++
inline Spline2DPartialDerivativeCreator::Spline2DPartialDerivativeCreator (
    SplineBuilder2DCache & builder_cache,
    SplineEvaluator2D const & evaluator
) 
```





**Parameters:**


* `builder_cache` A 2d spline builder cache. 
* `evaluator` A 2d spline evaluator. 




        

<hr>



### function create\_instance 

```C++
inline std::unique_ptr< IPartialDerivative < typename SplineEvaluator2D::batched_evaluation_domain_type, DerivativeDimension > > Spline2DPartialDerivativeCreator::create_instance (
    DConstFieldType field
) const
```



Create a pointer to an instance of the abstract class [**IPartialDerivative**](classIPartialDerivative.md). The type of the returned object will be determined when the pointer is dereferenced.




**Parameters:**


* `field` A field to be differentiated.



**Returns:**

A pointer to an instance of the [**IPartialDerivative**](classIPartialDerivative.md) class. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/spline_2d_partial_derivative.hpp`

