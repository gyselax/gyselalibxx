

# File paraconfpp.hpp

[**File List**](files.md) **>** [**paraconfpp**](dir_7700c957c7ac062c1c9c3e42e00e7e24.md) **>** [**paraconfpp.hpp**](paraconfpp_8hpp.md)

[Go to the documentation of this file](paraconfpp_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>

#include <paraconf.h>

double PCpp_double(PC_tree_t tree, std::string_view key);

long PCpp_int(PC_tree_t tree, std::string_view key);

int PCpp_len(PC_tree_t tree, std::string_view key);

bool PCpp_bool(PC_tree_t tree, std::string_view key);

std::string PCpp_string(PC_tree_t tree, std::string_view key);

PC_tree_t PCpp_get(PC_tree_t tree, std::string_view key);
```


