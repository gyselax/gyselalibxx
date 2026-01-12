

# Class SingleInterfaceDerivativesCalculatorCollection

**template &lt;class... Interfaces&gt;**



[**ClassList**](annotated.md) **>** [**SingleInterfaceDerivativesCalculatorCollection**](classSingleInterfaceDerivativesCalculatorCollection.md)



_A class to store a collection of interface derivative calculators templated on the interfaces._ [More...](#detailed-description)

* `#include <single_interface_derivatives_calculator_collection.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**SingleInterfaceDerivativesCalculatorCollection**](#function-singleinterfacederivativescalculatorcollection) ([**SingleInterfaceDerivativesCalculator**](classSingleInterfaceDerivativesCalculator.md)&lt; Interfaces &gt; const &... derivative\_calculators) <br>_Instantiate a_ [_**SingleInterfaceDerivativesCalculatorCollection**_](classSingleInterfaceDerivativesCalculatorCollection.md) _from a list of interface derivative calculators._ |
|  [**SingleInterfaceDerivativesCalculator**](classSingleInterfaceDerivativesCalculator.md)&lt; [**Interface**](structInterface.md) &gt; const & | [**get**](#function-get) () const<br>_Get a derivative calculator of the collection. The output cannot be copied. This operator only allows to get a temporary reference to call one of the operators of the_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _class._ |




























## Detailed Description


The class stores a constant reference of interface derivative calculators. It should not be use to copy the elements outside of the class but it should be use to access the operators of the stored interface derivative calculators.




**Template parameters:**


* `Interfaces` Types of interface that defined the interface derivative calculators.



**Warning:**

For each interface, only one interface derivative calculator should be defined.




**See also:** [**SingleInterfaceDerivativesCalculator**](classSingleInterfaceDerivativesCalculator.md). 



    
## Public Functions Documentation




### function SingleInterfaceDerivativesCalculatorCollection 

_Instantiate a_ [_**SingleInterfaceDerivativesCalculatorCollection**_](classSingleInterfaceDerivativesCalculatorCollection.md) _from a list of interface derivative calculators._
```C++
inline explicit SingleInterfaceDerivativesCalculatorCollection::SingleInterfaceDerivativesCalculatorCollection (
    SingleInterfaceDerivativesCalculator < Interfaces > const &... derivative_calculators
) 
```





**Parameters:**


* `derivative_calculators` [**Interface**](structInterface.md) derivative calculators. 




        

<hr>



### function get 

_Get a derivative calculator of the collection. The output cannot be copied. This operator only allows to get a temporary reference to call one of the operators of the_ [_**SingleInterfaceDerivativesCalculator**_](classSingleInterfaceDerivativesCalculator.md) _class._
```C++
template<class Interface>
inline SingleInterfaceDerivativesCalculator < Interface > const & SingleInterfaceDerivativesCalculatorCollection::get () const
```





**Template parameters:**


* [**Interface**](structInterface.md) The interface where the required interface derivative calculator is defined.



**Returns:**

The required interface derivative calculator as a constant reference. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/multipatch/interface_derivatives/single_interface_derivatives_calculator_collection.hpp`

