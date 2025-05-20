

# Class ConstantPartialDerivativeCreator

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**ConstantPartialDerivativeCreator**](classConstantPartialDerivativeCreator.md)



_A class to create a_ [_**ConstantPartialDerivative**_](classConstantPartialDerivative.md) _via a create\_instance function._[More...](#detailed-description)

* `#include <constant_partial_derivatives.hpp>`



Inherits the following classes: [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**ConstantPartialDerivativeCreator**](#function-constantpartialderivativecreator) (double deriv\_value) <br>_Create an instance of_ [_**ConstantPartialDerivativeCreator**_](classConstantPartialDerivativeCreator.md) _._ |
| virtual std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](#function-create_instance) (DConstField&lt; IdxRangeFull &gt; field) const<br>_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._ |


## Public Functions inherited from IPartialDerivativeCreator

See [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)

| Type | Name |
| ---: | :--- |
| virtual std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](classIPartialDerivativeCreator.md#function-create_instance) (DConstField&lt; IdxRangeFull &gt; field) const = 0<br>_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._ |






















































## Detailed Description




**Template parameters:**


* `IdxRangeFull` The index range of the field on which the operator acts (with all dimensions, batched and dimension of interest, used for inheritance). 
* `DerivativeDimension` The dimension on which the partial derivative is calculated. (used for inheritance). 




    
## Public Functions Documentation




### function ConstantPartialDerivativeCreator 

_Create an instance of_ [_**ConstantPartialDerivativeCreator**_](classConstantPartialDerivativeCreator.md) _._
```C++
inline explicit ConstantPartialDerivativeCreator::ConstantPartialDerivativeCreator (
    double deriv_value
) 
```





**Parameters:**


* `deriv_value` The value that should be returned as the constant value of the derivative. 




        

<hr>



### function create\_instance 

_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._
```C++
inline virtual std::unique_ptr< IPartialDerivative < IdxRangeFull, DerivativeDimension > > ConstantPartialDerivativeCreator::create_instance (
    DConstField< IdxRangeFull > field
) const
```





**Parameters:**


* `field` The field whose derivative should be calculated.



**Returns:**

A pointer to an [**IPartialDerivative**](classIPartialDerivative.md) object. 





        
Implements [*IPartialDerivativeCreator::create\_instance*](classIPartialDerivativeCreator.md#function-create_instance)


<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/constant_partial_derivatives.hpp`

