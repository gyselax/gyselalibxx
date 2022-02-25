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
PC_tree_t PCpp_get(PC_tree_t tree, std::string const& str, Args&&... args)
{
    return PC_get(tree, str.c_str(), std::forward<Args>(args)...);
}
