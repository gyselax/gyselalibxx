

# Class IPartialDerivativeCreator

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**IPartialDerivativeCreator**](classIPartialDerivativeCreator.md)



_An abstract class which provides a create\_instance function to instantiate an object of the_ [_**IPartialDerivative**_](classIPartialDerivative.md) _class where required._

* `#include <ipartial_derivative.hpp>`





Inherited by the following classes: [CentralFDMPartialDerivativeCreator](classCentralFDMPartialDerivativeCreator.md),  [CentralFDMPartialDerivativeWithBValueCreator](classCentralFDMPartialDerivativeWithBValueCreator.md)
































## Public Functions

| Type | Name |
| ---: | :--- |
| virtual std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](#function-create_instance) (DConstField&lt; IdxRangeFull &gt; field) const = 0<br>_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._ |




























## Public Functions Documentation




### function create\_instance 

_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._
```C++
virtual std::unique_ptr< IPartialDerivative < IdxRangeFull, DerivativeDimension > > IPartialDerivativeCreator::create_instance (
    DConstField< IdxRangeFull > field
) const = 0
```





**Parameters:**


* `field` A field to be passed to the constructor of [**IPartialDerivative**](classIPartialDerivative.md).



**Returns:**

A pointer to an [**IPartialDerivative**](classIPartialDerivative.md) object.




**See also:** [**IPartialDerivative**](classIPartialDerivative.md) 



        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/ipartial_derivative.hpp`

