

# Struct tensor\_tools::VectorIndexIdMap

**template &lt;char Id, class AssociatedVectorIndexSet&gt;**



[**ClassList**](annotated.md) **>** [**tensor\_tools**](namespacetensor__tools.md) **>** [**VectorIndexIdMap**](structtensor__tools_1_1VectorIndexIdMap.md)



_A class representing a vector index identifier._ [More...](#detailed-description)

* `#include <vector_index_tools.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef AssociatedVectorIndexSet | [**possible\_idx\_values**](#typedef-possible_idx_values)  <br>_The VectorIndexSet describing valid indices for this component._  |






## Public Static Attributes

| Type | Name |
| ---: | :--- |
|  constexpr char | [**id**](#variable-id)   = `Id`<br>_The character identifying the index._  |










































## Detailed Description


A vector is indexed at a certain position using an identifier (a character) which can take one of multiple possible values (types).




**Template parameters:**


* `Id` The character identifying the index. 
* `AssociatedVectorIndexSet` The VectorIndexSet describing the indices associated with this identifier. 




    
## Public Types Documentation




### typedef possible\_idx\_values 

_The VectorIndexSet describing valid indices for this component._ 
```C++
using tensor_tools::VectorIndexIdMap< Id, AssociatedVectorIndexSet >::possible_idx_values =  AssociatedVectorIndexSet;
```




<hr>
## Public Static Attributes Documentation




### variable id 

_The character identifying the index._ 
```C++
constexpr char tensor_tools::VectorIndexIdMap< Id, AssociatedVectorIndexSet >::id;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/vector_index_tools.hpp`

