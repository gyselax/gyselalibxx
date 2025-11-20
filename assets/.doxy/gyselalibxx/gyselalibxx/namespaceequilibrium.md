

# Namespace equilibrium



[**Namespace List**](namespaces.md) **>** [**equilibrium**](namespaceequilibrium.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  std::unique\_ptr&lt; [**IEquilibrium**](classIEquilibrium.md) &gt; | [**init\_from\_input**](#function-init_from_input) (IdxRangeSp idx\_range\_kinsp, PC\_tree\_t const & yaml\_input\_file) <br>_Determine the chosen equilibrium method from a YAML input file._  |




























## Public Functions Documentation




### function init\_from\_input 

_Determine the chosen equilibrium method from a YAML input file._ 
```C++
std::unique_ptr< IEquilibrium > equilibrium::init_from_input (
    IdxRangeSp idx_range_kinsp,
    PC_tree_t const & yaml_input_file
) 
```



Use the Algorithm.equilibrium key in a YAML input file to select the chosen equilibrium method and call the init\_from\_input method for that method.




**Parameters:**


* `idx_range_kinsp` Index range for the kinetic species 
* `yaml_input_file` YAML input file 



**Returns:**

A pointer to an equilibrium operator. A pointer is used for OOP. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/initialisation/iequilibrium.hpp`

