

# Class MpiScopeGuard



[**ClassList**](annotated.md) **>** [**MpiScopeGuard**](classMpiScopeGuard.md)



_RAII wrapper for MPI initialisation and finalisation._ [More...](#detailed-description)

* `#include <mpi_scope_guard.hpp>`





































## Public Functions

| Type | Name |
| ---: | :--- |
|   | [**MpiScopeGuard**](#function-mpiscopeguard-13) (int & argc, char \*\*& argv) noexcept<br>_Initialises the MPI environment using command-line arguments._  |
|   | [**MpiScopeGuard**](#function-mpiscopeguard-23) ([**MpiScopeGuard**](classMpiScopeGuard.md) const & rhs) = delete<br>_Deleted copy constructor._  |
|   | [**MpiScopeGuard**](#function-mpiscopeguard-33) ([**MpiScopeGuard**](classMpiScopeGuard.md) && rhs) noexcept<br>_Deleted move constructor._  |
|  [**MpiScopeGuard**](classMpiScopeGuard.md) & | [**operator=**](#function-operator) ([**MpiScopeGuard**](classMpiScopeGuard.md) const & rhs) = delete<br>_Deleted copy assignment operator._  |
|  [**MpiScopeGuard**](classMpiScopeGuard.md) & | [**operator=**](#function-operator_1) ([**MpiScopeGuard**](classMpiScopeGuard.md) && rhs) noexcept<br>_Deleted move assignment operator._  |
|   | [**~MpiScopeGuard**](#function-mpiscopeguard) () noexcept<br>_Finalises the MPI environment._  |




























## Detailed Description


This class manages the lifetime of the MPI environment using the Resource Acquisition Is Initialisation (RAII) idiom. Constructing an instance initialises MPI, and destroying the instance finalises MPI.


The class is non-copyable and non-movable to ensure unique ownership of the MPI runtime state within a process. 


    
## Public Functions Documentation




### function MpiScopeGuard [1/3]

_Initialises the MPI environment using command-line arguments._ 
```C++
MpiScopeGuard::MpiScopeGuard (
    int & argc,
    char **& argv
) noexcept
```



Calls `MPI_Init`, forwarding the provided command-line arguments to the MPI implementation.




**Parameters:**


* `argc` Number of command-line arguments. 
* `argv` Array of command-line argument strings.



**Note:**

The MPI implementation may modify `argc` and `argv`. 




**Note:**

This constructor assumes that MPI has not already been initialised elsewhere. 





        

<hr>



### function MpiScopeGuard [2/3]

_Deleted copy constructor._ 
```C++
MpiScopeGuard::MpiScopeGuard (
    MpiScopeGuard const & rhs
) = delete
```



Copying an `MpiScopeGuard` is disallowed because MPI lifetime ownership must remain unique. 


        

<hr>



### function MpiScopeGuard [3/3]

_Deleted move constructor._ 
```C++
MpiScopeGuard::MpiScopeGuard (
    MpiScopeGuard && rhs
) noexcept
```



Moving an `MpiScopeGuard` is disallowed because MPI lifetime ownership cannot be transferred safely. 


        

<hr>



### function operator= 

_Deleted copy assignment operator._ 
```C++
MpiScopeGuard & MpiScopeGuard::operator= (
    MpiScopeGuard const & rhs
) = delete
```




<hr>



### function operator= 

_Deleted move assignment operator._ 
```C++
MpiScopeGuard & MpiScopeGuard::operator= (
    MpiScopeGuard && rhs
) noexcept
```




<hr>



### function ~MpiScopeGuard 

_Finalises the MPI environment._ 
```C++
MpiScopeGuard::~MpiScopeGuard () noexcept
```



Calls `MPI_Finalize` if MPI was successfully initialised by this guard instance. 


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/mpi_parallelisation/mpi_scope_guard.hpp`

