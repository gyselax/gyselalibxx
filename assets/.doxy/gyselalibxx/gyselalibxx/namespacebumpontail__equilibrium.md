

# Namespace bumpontail\_equilibrium



[**Namespace List**](namespaces.md) **>** [**bumpontail\_equilibrium**](namespacebumpontail__equilibrium.md)










































## Public Functions

| Type | Name |
| ---: | :--- |
|  [**BumpontailEquilibrium**](classBumpontailEquilibrium.md) | [**init\_from\_input**](#function-init_from_input) (IdxRangeSp idx\_range\_kinsp, PC\_tree\_t const & yaml\_input\_file) <br>_Read the density, temperature and mean velocity required to initialise the bump-on-tail Maxwellian in a YAML input file._  |




























## Public Functions Documentation




### function init\_from\_input 

_Read the density, temperature and mean velocity required to initialise the bump-on-tail Maxwellian in a YAML input file._ 
```C++
BumpontailEquilibrium bumpontail_equilibrium::init_from_input (
    IdxRangeSp idx_range_kinsp,
    PC_tree_t const & yaml_input_file
) 
```





**Parameters:**


* `idx_range_kinsp` Index range for the kinetic species 
* `yaml_input_file` YAML input file 



**Returns:**

an instance of Maxwellian distribution function. 





        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/geometryXVx/initialisation/bumpontailequilibrium.hpp`

