

# File preprocessor.hpp

[**File List**](files.md) **>** [**src**](dir_68267d1309a1af8e8297ef4c3efbcdba.md) **>** [**utils**](dir_313caf1132e152dd9b58bea13a4052ca.md) **>** [**preprocessor.hpp**](preprocessor_8hpp.md)

[Go to the documentation of this file](preprocessor_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once


#define GSLX_UNUSED(an_expression) static_cast<void>(an_expression)



#define _GSLX_UTILITY_STRINGIFY(an_expression) #an_expression
#define GSLX_UTILITY_STRINGIFY(an_expression) _GSLX_UTILITY_STRINGIFY(an_expression)
```


