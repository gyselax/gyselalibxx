

# Class CentralFDMPartialDerivativeWithBValueCreator

**template &lt;class IdxRangeFull, class DerivativeDimension&gt;**



[**ClassList**](annotated.md) **>** [**CentralFDMPartialDerivativeWithBValueCreator**](classCentralFDMPartialDerivativeWithBValueCreator.md)



_A class which stores information necessary to create a pointer to an instance of the_ [_**CentralFDMPartialDerivativeWithBValue**_](classCentralFDMPartialDerivativeWithBValue.md) _class._[More...](#detailed-description)

* `#include <central_fdm_partial_derivatives_with_boundary_values.hpp>`



Inherits the following classes: [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)






















































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**CentralFDMPartialDerivativeWithBValueCreator**](#function-centralfdmpartialderivativewithbvaluecreator) (double bvalue\_left, double bvalue\_right) <br>_Construct an instance of the_ [_**CentralFDMPartialDerivativeWithBValueCreator**_](classCentralFDMPartialDerivativeWithBValueCreator.md) _class._ |
|  std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](#function-create_instance) (DConstFieldType field\_ref) const<br> |


## Public Functions inherited from IPartialDerivativeCreator

See [IPartialDerivativeCreator](classIPartialDerivativeCreator.md)

| Type | Name |
| ---: | :--- |
| virtual std::unique\_ptr&lt; [**IPartialDerivative**](classIPartialDerivative.md)&lt; IdxRangeFull, DerivativeDimension &gt; &gt; | [**create\_instance**](classIPartialDerivativeCreator.md#function-create_instance) (DConstField&lt; IdxRangeFull &gt; field) const = 0<br>_Create an instance of a pointer to an_ [_**IPartialDerivative**_](classIPartialDerivative.md) _object._ |






















































## Detailed Description


This class allows an instance of the [**CentralFDMPartialDerivativeWithBValue**](classCentralFDMPartialDerivativeWithBValue.md) class to be instantiated where necessary. Typically, the [**CentralFDMPartialDerivativeCreator**](classCentralFDMPartialDerivativeCreator.md) is instantiated in the initialisation of the simulation, and the corresponding [**CentralFDMPartialDerivativeWithBValue**](classCentralFDMPartialDerivativeWithBValue.md) object is instantiated where computing partial derivatives is required. 

**Template parameters:**


* `IdxRangeFull` The index range of the field on which the operator acts (with all dimensions, batched and dimension of interest). 
* `DerivativeDimension` The dimension on which the partial derivative is calculated. 




    
## Public Functions Documentation




### function CentralFDMPartialDerivativeWithBValueCreator 

_Construct an instance of the_ [_**CentralFDMPartialDerivativeWithBValueCreator**_](classCentralFDMPartialDerivativeWithBValueCreator.md) _class._
```C++
inline CentralFDMPartialDerivativeWithBValueCreator::CentralFDMPartialDerivativeWithBValueCreator (
    double bvalue_left,
    double bvalue_right
) 
```





**Parameters:**


* `bvalue_left` The left boundary value. 
* `bvalue_right` The right boundary value. 




        

<hr>



### function create\_instance 

```C++
inline std::unique_ptr< IPartialDerivative < IdxRangeFull, DerivativeDimension > > CentralFDMPartialDerivativeWithBValueCreator::create_instance (
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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/math_tools/central_fdm_partial_derivatives_with_boundary_values.hpp`

