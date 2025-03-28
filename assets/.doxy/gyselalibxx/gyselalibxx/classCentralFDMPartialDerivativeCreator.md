

# Class CentralFDMPartialDerivativeCreator

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**CentralFDMPartialDerivativeCreator**](classCentralFDMPartialDerivativeCreator.md)



_A class which stores information necessary to create a pointer to an instance of the_ [_**CentralFDMPartialDerivative**_](classCentralFDMPartialDerivative.md) _class._[More...](#detailed-description)

* `#include <central_fdm_partial_derivatives.hpp>`



Inherits the following classes: [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CentralFDMPartialDerivativeCreator**](#function-centralfdmpartialderivativecreator) () = default<br>_Construct an instance of the_ [_**CentralFDMPartialDerivativeCreator**_](classCentralFDMPartialDerivativeCreator.md) _class._ |
|  std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](#function-create_instance) (DConstFieldType field\_ref) const<br> |


## Public Functions inherited from IPartialDerivativeCreator

See [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)

| Type | Name |
| ---: | :--- |
| virtual std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](classIPartialDerivativeCreator.md#function-create_instance) (DConstField&lt; IdxRangeFull &gt; field) const = 0<br>_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._ |






















































## Detailed Description


This class allows an instance of the [**CentralFDMPartialDerivative**](classCentralFDMPartialDerivative.md) class to be instantiated where necessary. Typically, the [**CentralFDMPartialDerivativeCreator**](classCentralFDMPartialDerivativeCreator.md) is instantiated in the initialisation of the simulation, and the corresponding [**CentralFDMPartialDerivative**](classCentralFDMPartialDerivative.md) object is instantiated where computing partial derivatives is required. 

**Template parameters:**


* `IdxRangeFull` The index range of the field on which the operator acts (with all dimensions, batched and dimension of interest). 
* `DerivativeDimension` The dimension on which the partial derivative is calculated. 




    
## Public Functions Documentation




### function CentralFDMPartialDerivativeCreator 

_Construct an instance of the_ [_**CentralFDMPartialDerivativeCreator**_](classCentralFDMPartialDerivativeCreator.md) _class._
```C++
CentralFDMPartialDerivativeCreator::CentralFDMPartialDerivativeCreator () = default
```




<hr>



### function create\_instance 

```C++
inline std::unique_ptr< IPartialDerivative < IdxRangeFull, DerivativeDimension > > CentralFDMPartialDerivativeCreator::create_instance (
    DConstFieldType field_ref
) const
```



Create a pointer to an instance of the abstract class [**IPartialDerivative**](classIPartialDerivative.md). The type of the returned object will be determined when the pointer is dereferenced.




**Parameters:**


* `field_ref` A field to be differentiated.



**Returns:**

A pointer to an instance of the [**IPartialDerivative**](classIPartialDerivative.md) class. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/central_fdm_partial_derivatives.hpp`

