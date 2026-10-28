

# File paraconfpp.hpp

[**File List**](files.md) **>** [**paraconfpp**](dir_7700c957c7ac062c1c9c3e42e00e7e24.md) **>** [**paraconfpp.hpp**](paraconfpp_8hpp.md)

[Go to the documentation of this file](paraconfpp_8hpp.md)


```C++
// SPDX-License-Identifier: MIT

#pragma once

#include <filesystem>
#include <stdexcept>
#include <string>
#include <utility>

#include <paraconf.h>

template <class... Args>
double PCpp_double(PC_tree_t tree, std::string const& str, Args&&... args)
{
    double data = 0.;
    if (PC_status_t s = PC_double(PC_get(tree, str.c_str(), std::forward<Args>(args)...), &data);
        s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

template <class... Args>
long PCpp_int(PC_tree_t tree, std::string const& str, Args&&... args)
{
    long data = 0;
    if (PC_status_t s = PC_int(PC_get(tree, str.c_str(), std::forward<Args>(args)...), &data);
        s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

template <class... Args>
int PCpp_len(PC_tree_t tree, std::string const& str, Args&&... args)
{
    int data = 0;
    if (PC_status_t s = PC_len(PC_get(tree, str.c_str(), std::forward<Args>(args)...), &data);
        s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

template <class... Args>
bool PCpp_bool(PC_tree_t tree, std::string const& str, Args&&... args)
{
    int data = 0;
    if (PC_status_t s = PC_bool(PC_get(tree, str.c_str(), std::forward<Args>(args)...), &data);
        s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

template <class... Args>
std::string PCpp_string(PC_tree_t tree, std::string const& str, Args&&... args)
{
    char* c_str = nullptr;
    if (PC_status_t s = PC_string(PC_get(tree, str.c_str(), std::forward<Args>(args)...), &c_str);
        s != PC_OK) {
        if (c_str != nullptr) {
            free(c_str);
        }
        throw std::runtime_error(PC_errmsg());
    }
    std::string const str_out(c_str);
    free(c_str);
    return str_out;
}

template <class... Args>
PC_tree_t PCpp_get(PC_tree_t tree, std::string const& str, Args&&... args)
{
    return PC_get(tree, str.c_str(), std::forward<Args>(args)...);
}
```


