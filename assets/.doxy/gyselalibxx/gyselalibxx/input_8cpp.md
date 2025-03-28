

# File input.cpp



[**FileList**](files.md) **>** [**io**](dir_c184e51c84f2c3f0345bbc8a0d75d3e1.md) **>** [**input.cpp**](input_8cpp.md)

[Go to the source code of this file](input_8cpp_source.md)



* `#include <cstdlib>`
* `#include <filesystem>`
* `#include <fstream>`
* `#include <iostream>`
* `#include "input.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**parse\_executable\_arguments**](#function-parse_executable_arguments) (PC\_tree\_t & conf\_gyselalibxx, long int & iter\_start, int argc, char \*\* argv, char const \*const params\_yaml) <br>_Extract the paraconf configuration and the restart iteration from the executable arguments._  |
|  PC\_tree\_t | [**parse\_executable\_arguments**](#function-parse_executable_arguments) (int argc, char \*\* argv, char const \*const params\_yaml) <br>_Extract the paraconf configuration from the executable arguments._  |




























## Public Functions Documentation




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
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/io/input.cpp`

