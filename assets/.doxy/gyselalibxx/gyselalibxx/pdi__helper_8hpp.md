

# File pdi\_helper.hpp



[**FileList**](files.md) **>** [**io**](dir_c184e51c84f2c3f0345bbc8a0d75d3e1.md) **>** [**pdi\_helper.hpp**](pdi__helper_8hpp.md)

[Go to the source code of this file](pdi__helper_8hpp_source.md)



* `#include <ddc/ddc.hpp>`
* `#include <ddc/pdi.hpp>`
* `#include <pdi.h>`
* `#include "vector_field.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**PDI\_expose\_idx\_range**](#function-pdi_expose_idx_range) (IdxRange&lt; Grids... &gt; index\_range, std::string name) <br> |
|  void | [**PDI\_expose\_vector\_field**](#function-pdi_expose_vector_field) (std::string const & name\_stem, [**VectorConstField**](classVectorField.md)&lt; ElementType, IdxRangeType, VectorIndexSet&lt; IndexTag... &gt;, Kokkos::HostSpace &gt; out\_vector, StringType const &... name\_suffixes) <br>_A helper function to expose a vector field to PDI ready for output to a file._  |
|  void | [**PDI\_get\_arrays**](#function-pdi_get_arrays) (std::string const & event\_name, std::string const & name, std::vector&lt; T &gt; & out\_vector, Args &... input\_args) <br> |




























## Public Functions Documentation




### function PDI\_expose\_idx\_range 

```C++
template<class... Grids>
void PDI_expose_idx_range (
    IdxRange< Grids... > index_range,
    std::string name
) 
```




<hr>



### function PDI\_expose\_vector\_field 

_A helper function to expose a vector field to PDI ready for output to a file._ 
```C++
template<class ElementType, class IdxRangeType, class... IndexTag, class... StringType>
void PDI_expose_vector_field (
    std::string const & name_stem,
    VectorConstField < ElementType, IdxRangeType, VectorIndexSet< IndexTag... >, Kokkos::HostSpace > out_vector,
    StringType const &... name_suffixes
) 
```





**Parameters:**


* `name_stem` The prefix for the names of the elements of the vector field. 
* `out_vector` The vector field to be exposed. 
* `name_suffixes` The suffixes for the names of the elements of the vector field. There must be the same number of suffixes as there are dimensions in the vector field. 




        

<hr>



### function PDI\_get\_arrays 

```C++
template<class T, class... Args>
void PDI_get_arrays (
    std::string const & event_name,
    std::string const & name,
    std::vector< T > & out_vector,
    Args &... input_args
) 
```



A helper function to read an unknown number of arrays from a file using PDI. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/io/pdi_helper.hpp`

