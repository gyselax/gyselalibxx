

# File paraconfpp.hpp



[**FileList**](files.md) **>** [**paraconfpp**](dir_7700c957c7ac062c1c9c3e42e00e7e24.md) **>** [**paraconfpp.hpp**](paraconfpp_8hpp.md)

[Go to the source code of this file](paraconfpp_8hpp_source.md)



* `#include <string_view>`
* `#include <paraconf.h>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|  bool | [**PCpp\_bool**](#function-pcpp_bool) (PC\_tree\_t tree, std::string\_view key) <br>_Retrieves a configuration value as a boolean._  |
|  double | [**PCpp\_double**](#function-pcpp_double) (PC\_tree\_t tree, std::string\_view key) <br>_Retrieves a configuration value as a double._  |
|  PC\_tree\_t | [**PCpp\_get**](#function-pcpp_get) (PC\_tree\_t tree, std::string\_view key) <br>_Retrieves a configuration entry._  |
|  long | [**PCpp\_int**](#function-pcpp_int) (PC\_tree\_t tree, std::string\_view key) <br>_Retrieves a configuration value as an integer._  |
|  int | [**PCpp\_len**](#function-pcpp_len) (PC\_tree\_t tree, std::string\_view key) <br>_Retrieves the length of a configuration entry._  |
|  std::string | [**PCpp\_string**](#function-pcpp_string) (PC\_tree\_t tree, std::string\_view key) <br>_Retrieves a configuration value as a string._  |




























## Public Functions Documentation




### function PCpp\_bool 

_Retrieves a configuration value as a boolean._ 
```C++
bool PCpp_bool (
    PC_tree_t tree,
    std::string_view key
) 
```





**Parameters:**


* `tree` The configuration tree. 
* `key` The key identifying the configuration entry.



**Returns:**

The configuration value as a boolean.




**Exception:**


* `std::runtime_error` If the value cannot be retrieved or converted. 




        

<hr>



### function PCpp\_double 

_Retrieves a configuration value as a double._ 
```C++
double PCpp_double (
    PC_tree_t tree,
    std::string_view key
) 
```





**Parameters:**


* `tree` The configuration tree. 
* `key` The key identifying the configuration entry.



**Returns:**

The configuration value as a double.




**Exception:**


* `std::runtime_error` If the value cannot be retrieved or converted. 




        

<hr>



### function PCpp\_get 

_Retrieves a configuration entry._ 
```C++
PC_tree_t PCpp_get (
    PC_tree_t tree,
    std::string_view key
) 
```





**Parameters:**


* `tree` The configuration tree. 
* `key` The key identifying the configuration entry.



**Returns:**

The configuration entry, or an invalid handle if the entry cannot be found. 





        

<hr>



### function PCpp\_int 

_Retrieves a configuration value as an integer._ 
```C++
long PCpp_int (
    PC_tree_t tree,
    std::string_view key
) 
```





**Parameters:**


* `tree` The configuration tree. 
* `key` The key identifying the configuration entry.



**Returns:**

The configuration value as a long integer.




**Exception:**


* `std::runtime_error` If the value cannot be retrieved or converted. 




        

<hr>



### function PCpp\_len 

_Retrieves the length of a configuration entry._ 
```C++
int PCpp_len (
    PC_tree_t tree,
    std::string_view key
) 
```





**Parameters:**


* `tree` The configuration tree. 
* `key` The key identifying the configuration entry.



**Returns:**

The length of the configuration entry.




**Exception:**


* `std::runtime_error` If the length cannot be retrieved. 




        

<hr>



### function PCpp\_string 

_Retrieves a configuration value as a string._ 
```C++
std::string PCpp_string (
    PC_tree_t tree,
    std::string_view key
) 
```





**Parameters:**


* `tree` The configuration tree. 
* `key` The key identifying the configuration entry.



**Returns:**

The configuration value as a string.




**Exception:**


* `std::runtime_error` If the value cannot be retrieved or converted. 




        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/paraconfpp/paraconfpp.hpp`

