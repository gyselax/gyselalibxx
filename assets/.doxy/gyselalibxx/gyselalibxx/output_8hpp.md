

# File output.hpp



[**FileList**](files.md) **>** [**io**](dir_c184e51c84f2c3f0345bbc8a0d75d3e1.md) **>** [**output.hpp**](output_8hpp.md)

[Go to the source code of this file](output_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/pdi.hpp>`
* `#include "ddc_aliases.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**expose\_mesh\_to\_pdi**](#function-expose_mesh_to_pdi) (std::string pdi\_name, IdxRange&lt; Mesh &gt; idx\_range) <br>_Expose a IdxRange to PDI._  |




























## Public Functions Documentation




### function expose\_mesh\_to\_pdi 

_Expose a IdxRange to PDI._ 
```C++
template<class Mesh>
void expose_mesh_to_pdi (
    std::string pdi_name,
    IdxRange< Mesh > idx_range
) 
```





**Parameters:**


* `pdi_name` The name given in PDI. 
* `idx_range` IdxRange that is exposed 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/io/output.hpp`

