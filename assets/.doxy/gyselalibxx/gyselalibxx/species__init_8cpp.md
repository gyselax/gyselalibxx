

# File species\_init.cpp



[**FileList**](files.md) **>** [**speciesinfo**](dir_661be8452a62f1b4720eb6eb57123ae7.md) **>** [**species\_init.cpp**](species__init_8cpp.md)

[Go to the source code of this file](species__init_8cpp_source.md)



* `#include "species_init.hpp"`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  void | [**init\_all\_species**](#function-init_all_species) (IdxRangeSp & idx\_range\_kinsp, IdxRangeSp & idx\_range\_fluidsp, PC\_tree\_t conf\_gyselalibxx, int nb\_kinspecies, int nb\_fluidspecies) <br>_Initialise the species domain._  |
|  IdxRangeSp | [**init\_kinetic\_species**](#function-init_kinetic_species) () <br>_Initialise the kinetic species domain._  |
|  IdxRangeSp | [**init\_species**](#function-init_species) (PC\_tree\_t conf\_gyselalibxx) <br>_Initialise the species domain in the case of adiabatic and kinetic species._  |
|  void | [**init\_species\_withfluid**](#function-init_species_withfluid) (IdxRangeSp & idx\_range\_kinsp, IdxRangeSp & idx\_range\_fluidsp, PC\_tree\_t conf\_gyselalibxx) <br>_Initialise the species domain in the specific case of fluid species added to adiabatic and kinetic species._  |




























## Public Functions Documentation




### function init\_all\_species 

_Initialise the species domain._ 
```C++
void init_all_species (
    IdxRangeSp & idx_range_kinsp,
    IdxRangeSp & idx_range_fluidsp,
    PC_tree_t conf_gyselalibxx,
    int nb_kinspecies,
    int nb_fluidspecies
) 
```





**Parameters:**


* `idx_range_kinsp` kinetic species domain 
* `idx_range_fluidsp` fluid species domain 
* `conf_gyselalibxx` is the YAML input file 
* `nb_kinspecies` number of kinetic species 
* `nb_fluidspecies` number of fluid species 




        

<hr>



### function init\_kinetic\_species 

_Initialise the kinetic species domain._ 
```C++
IdxRangeSp init_kinetic_species () 
```





**Returns:**

the kinetic species domain 





        

<hr>



### function init\_species 

_Initialise the species domain in the case of adiabatic and kinetic species._ 
```C++
IdxRangeSp init_species (
    PC_tree_t conf_gyselalibxx
) 
```





**Parameters:**


* `conf_gyselalibxx` is the YAML input file 



**Returns:**

the kinetic species domain 





        

<hr>



### function init\_species\_withfluid 

_Initialise the species domain in the specific case of fluid species added to adiabatic and kinetic species._ 
```C++
void init_species_withfluid (
    IdxRangeSp & idx_range_kinsp,
    IdxRangeSp & idx_range_fluidsp,
    PC_tree_t conf_gyselalibxx
) 
```





**Parameters:**


* `idx_range_kinsp` kinetic species domain 
* `idx_range_fluidsp` fluid species domain 
* `conf_gyselalibxx` is the YAML input file 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/speciesinfo/species_init.cpp`

