

# File input.hpp



[**FileList**](files.md) **>** [**io**](dir_c184e51c84f2c3f0345bbc8a0d75d3e1.md) **>** [**input.hpp**](input_8hpp.md)

[Go to the source code of this file](input_8hpp_source.md)



* `#include <sstream>`
* `#include <ddc/ddc.hpp>`
* `#include <paraconf.h>`
* `#include <pdi.h>`
* `#include "ddc_aliases.hpp"`
* `#include "mesh_builder.hpp"`
* `#include "non_uniform_interpolation_points.hpp"`
* `#include "paraconfpp.hpp"`
* `#include "pdi_helper.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  IdxRange&lt; Grid1D &gt; | [**init\_pseudo\_uniform\_spline\_dependent\_idx\_range**](#function-init_pseudo_uniform_spline_dependent_idx_range) (PC\_tree\_t const & conf\_gyselalibxx, std::string const & mesh\_identifier) <br> |
|  IdxRange&lt; Grid1D &gt; | [**init\_spline\_dependent\_idx\_range**](#function-init_spline_dependent_idx_range) (PC\_tree\_t const & conf\_gyselalibxx, std::string const & mesh\_identifier) <br> |
|  void | [**parse\_executable\_arguments**](#function-parse_executable_arguments) (PC\_tree\_t & conf\_gyselalibxx, long int & iter\_start, int argc, char \*\* argv, char const \*const params\_yaml) <br>_Extract the paraconf configuration and the restart iteration from the executable arguments._  |
|  PC\_tree\_t | [**parse\_executable\_arguments**](#function-parse_executable_arguments) (int argc, char \*\* argv, char const \*const params\_yaml) <br>_Extract the paraconf configuration from the executable arguments._  |




























## Public Functions Documentation




### function init\_pseudo\_uniform\_spline\_dependent\_idx\_range 

```C++
template<class Grid1D, class BSplines, class InterpPointInitMethod>
inline IdxRange< Grid1D > init_pseudo_uniform_spline_dependent_idx_range (
    PC_tree_t const & conf_gyselalibxx,
    std::string const & mesh_identifier
) 
```



Initialise an index range which will serve as an interpolation index range for splines.


The index range is initialised using information from an input yaml file. This function should be used for non-uniform B-splines, but it initialises the break points uniformly. Such splines are referred to as pseudo-uniform as the cells on which the polynomials are defined are uniform. However they are not strictly uniform as multiple knots will be found at the same position at the boundary.


The information to be read from the file is:
* .SplineMesh.&lt;mesh\_identifier&gt;\_min
* .SplineMesh.&lt;mesh\_identifier&gt;\_max
* .SplineMesh.&lt;mesh\_identifier&gt;\_ncells




The interpolation index range is then created using the specified method. 


        

<hr>



### function init\_spline\_dependent\_idx\_range 

```C++
template<class Grid1D, class BSplines, class InterpPointInitMethod>
inline IdxRange< Grid1D > init_spline_dependent_idx_range (
    PC_tree_t const & conf_gyselalibxx,
    std::string const & mesh_identifier
) 
```



Initialise an index range which will serve as an interpolation index range for splines.


The index range is initialised using information from an input yaml file. If the B-splines are uniform then the information to be read is:
* .SplineMesh.&lt;mesh\_identifier&gt;\_min
* .SplineMesh.&lt;mesh\_identifier&gt;\_max
* .SplineMesh.&lt;mesh\_identifier&gt;\_ncells




If the B-splines are non-uniform then the information to be read is:
* .SplineMesh.&lt;mesh\_identifier&gt;\_MeshFile




This string indicates the name of a file which contains the knots of the bspline.


This information is used to initialise the B-splines. The interpolation index range is then created using the specified method. 


        

<hr>



### function parse\_executable\_arguments 

_Extract the paraconf configuration and the restart iteration from the executable arguments._ 
```C++
void parse_executable_arguments (
    PC_tree_t & conf_gyselalibxx,
    long int & iter_start,
    int argc,
    char ** argv,
    char const *const params_yaml
) 
```





**Parameters:**


* `conf_gyselalibxx` The paraconf configuration describing the simulation. 
* `iter_start` The index of the iteration from which the simulation should restart. 
* `argc` The number of arguments passed to the executable. 
* `argv` The arguments passed to the executable. 
* `params_yaml` The default parameters for the yaml file. 




        

<hr>



### function parse\_executable\_arguments 

_Extract the paraconf configuration from the executable arguments._ 
```C++
PC_tree_t parse_executable_arguments (
    int argc,
    char ** argv,
    char const *const params_yaml
) 
```





**Parameters:**


* `argc` The number of arguments passed to the executable. 
* `argv` The arguments passed to the executable. 
* `params_yaml` The default parameters for the yaml file.



**Returns:**

The paraconf configuration describing the simulation. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/io/input.hpp`

