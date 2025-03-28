

# Class Spline1DPartialDerivativeCreator

**template &lt;class Spline1DBuilder, class Spline1DEvaluator&gt;**



[**ClassList**](annotated.md) **>** [**Spline1DPartialDerivativeCreator**](classSpline1DPartialDerivativeCreator.md)



_A class which stores information necessary to create a pointer to an instance of the_ [_**Spline1DPartialDerivative**_](classSpline1DPartialDerivative.md) _class._[More...](#detailed-description)

* `#include <spline_1d_partial_derivative.hpp>`



Inherits the following classes: [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**Spline1DPartialDerivativeCreator**](#function-spline1dpartialderivativecreator) (Spline1DBuilder const & builder, Spline1DEvaluator const & evaluator) <br>_Construct an instance of the_ [_**Spline1DPartialDerivativeCreator**_](classSpline1DPartialDerivativeCreator.md) _class._ |
|  std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; typename Spline1DBuilder::batched\_interpolation\_domain\_type, typename Spline1DBuilder::continuous\_dimension\_type &gt; &gt; | [**create\_instance**](#function-create_instance) (DConstFieldType field) const<br> |


## Public Functions inherited from IPartialDerivativeCreator

See [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)

| Type | Name |
| ---: | :--- |
| virtual std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](classIPartialDerivativeCreator.md#function-create_instance) (DConstField&lt; IdxRangeFull &gt; field) const = 0<br>_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._ |






















































## Detailed Description


This class allows an instance of the [**Spline1DPartialDerivative**](classSpline1DPartialDerivative.md) class to be instantiated where necessary. Typically, the [**Spline1DPartialDerivativeCreator**](classSpline1DPartialDerivativeCreator.md) is instantiated in the initialisation of the simulation, and the corresponding [**Spline1DPartialDerivative**](classSpline1DPartialDerivative.md) object is instantiated where computing partial derivatives is required. 

**Template parameters:**


* `Spline1DBuilder` A 1D spline builder. 
* `Spline1DEvaluator` A 1D spline evaluator. 




    
## Public Functions Documentation




### function Spline1DPartialDerivativeCreator 

_Construct an instance of the_ [_**Spline1DPartialDerivativeCreator**_](classSpline1DPartialDerivativeCreator.md) _class._
```C++
inline Spline1DPartialDerivativeCreator::Spline1DPartialDerivativeCreator (
    Spline1DBuilder const & builder,
    Spline1DEvaluator const & evaluator
) 
```





**Parameters:**


* `builder` A 1d spline builder. 
* `evaluator` A 1d spline evaluator. 




        

<hr>



### function create\_instance 

```C++
inline std::unique_ptr< IPartialDerivative < typename Spline1DBuilder::batched_interpolation_domain_type, typename Spline1DBuilder::continuous_dimension_type > > Spline1DPartialDerivativeCreator::create_instance (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/spline_1d_partial_derivative.hpp`

