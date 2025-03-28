

# File preprocessor.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**preprocessor.hpp**](preprocessor_8hpp.md)

[Go to the source code of this file](preprocessor_8hpp_source.md)



































































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**GSLX\_UNUSED**](preprocessor_8hpp.md#define-gslx_unused) (an\_expression) `static\_cast&lt;void&gt;(an\_expression)`<br> |
| define  | [**GSLX\_UTILITY\_STRINGIFY**](preprocessor_8hpp.md#define-gslx_utility_stringify) (an\_expression) `\_GSLX\_UTILITY\_STRINGIFY(an\_expression)`<br> |
| define  | [**\_GSLX\_UTILITY\_STRINGIFY**](preprocessor_8hpp.md#define-_gslx_utility_stringify) (an\_expression) `#an\_expression`<br> |

## Macro Definition Documentation





### define GSLX\_UNUSED 

```C++
#define GSLX_UNUSED (
    an_expression
) `static_cast<void>(an_expression)`
```



GSLX\_UNUSED


Prevent a compiler warning for an unused variable.


Example usage:  


        

<hr>



### define GSLX\_UTILITY\_STRINGIFY 

```C++
#define GSLX_UTILITY_STRINGIFY (
    an_expression
) `_GSLX_UTILITY_STRINGIFY(an_expression)`
```




<hr>



### define \_GSLX\_UTILITY\_STRINGIFY 

```C++
#define _GSLX_UTILITY_STRINGIFY (
    an_expression
) `#an_expression`
```



GSLX\_UTILITY\_STRINGIFY


This macro will convert the text as argument to a string literal.


Example usage:  Expands to :  


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/preprocessor.hpp`

