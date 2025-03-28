

# Struct tensor\_tools::GetCovariantDims&lt; VectorIndexSet&lt; Dims... &gt; &gt;

**template &lt;class... Dims&gt;**



[**ClassList**](annotated.md) **>** [**tensor\_tools**](namespacetensor__tools.md) **>** [**GetCovariantDims&lt; VectorIndexSet&lt; Dims... &gt; &gt;**](structtensor__tools_1_1GetCovariantDims_3_01VectorIndexSet_3_01Dims_8_8_8_01_4_01_4.md)



_A class to get a VectorIndexSet containing only covariant dimensions._ [More...](#detailed-description)

* `#include <vector_index_tools.hpp>`

















## Public Types

| Type | Name |
| ---: | :--- |
| typedef VectorIndexSet&lt; std::conditional\_t&lt; Dims::IS\_COVARIANT, Dims, typename Dims::Dual &gt;... &gt; | [**type**](#typedef-type)  <br>_The type of the VectorIndexSet containing only covariant dimensions._  |
















































## Detailed Description




**Template parameters:**


* `AnyVectorIndexSet` The original VectorIndexSet. 




    
## Public Types Documentation




### typedef type 

_The type of the VectorIndexSet containing only covariant dimensions._ 
```C++
using tensor_tools::GetCovariantDims< VectorIndexSet< Dims... > >::type =  VectorIndexSet<std::conditional_t<Dims::IS_COVARIANT, Dims, typename Dims::Dual>...>;
```




<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/data_types/vector_index_tools.hpp`

