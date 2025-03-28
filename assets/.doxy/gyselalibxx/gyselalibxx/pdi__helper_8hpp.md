

# File pdi\_helper.hpp



[**FileList**](files.md) **>** [**io**](dir_c184e51c84f2c3f0345bbc8a0d75d3e1.md) **>** [**pdi\_helper.hpp**](pdi__helper_8hpp.md)

[Go to the source code of this file](pdi__helper_8hpp_source.md)



* `#include <pdi.h>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**PDI\_expose\_idx\_range**](#function-pdi_expose_idx_range) (IdxRange&lt; Grids... &gt; index\_range, std::string name) <br> |
|  void | [**PDI\_get\_arrays**](#function-pdi_get_arrays) (std::string const & event\_name, std::string const & name, std::vector&lt; [**T**](structT.md) &gt; & out\_vector, Args &... input\_args) <br> |




























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

