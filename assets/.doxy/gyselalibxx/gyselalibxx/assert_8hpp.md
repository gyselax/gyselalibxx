

# File assert.hpp



[**FileList**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**assert.hpp**](assert_8hpp.md)

[Go to the source code of this file](assert_8hpp_source.md)



* `#include "preprocessor.hpp"`
* `#include <cstdlib>`













## Namespaces

| Type | Name |
| ---: | :--- |
| namespace | [**gslx**](namespacegslx.md) <br> |
| namespace | [**error**](namespacegslx_1_1error.md) <br> |



















































## Macros

| Type | Name |
| ---: | :--- |
| define  | [**GSLX\_ASSERT**](assert_8hpp.md#define-gslx_assert) (an\_expression) `/* multi line expression */`<br> |
| define  | [**GSLX\_ASSERT\_ENABLED**](assert_8hpp.md#define-gslx_assert_enabled)  `1`<br> |
| define  | [**GSLX\_DEBUG\_BREAK**](assert_8hpp.md#define-gslx_debug_break) () `std::exit(EXIT\_FAILURE)`<br> |

## Macro Definition Documentation





### define GSLX\_ASSERT 

```C++
#define GSLX_ASSERT (
    an_expression
) `/* multi line expression */`
```




<hr>



### define GSLX\_ASSERT\_ENABLED 

```C++
#define GSLX_ASSERT_ENABLED `1`
```




<hr>



### define GSLX\_DEBUG\_BREAK 

```C++
#define GSLX_DEBUG_BREAK (
    
) `std::exit(EXIT_FAILURE)`
```



GSLX\_DEBUG\_BREAK


This macro will stop the program.


Example usage:  


        

<hr>

------------------------------
The documentation for this class was generated from the following file `/home/runner/work/gyselalibxx/gyselalibxx/code_branch/src/utils/assert.hpp`

