// SPDX-License-Identifier: MIT

#include <cstdlib>
#include <stdexcept>
#include <string>

#include <paraconf.h>

#include "paraconfpp.hpp"

double PCpp_double(PC_tree_t tree, std::string const& str)
{
    double data = 0.;
    if (PC_status_t s = PC_double(PC_get(tree, "%s", str.c_str()), &data); s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

long PCpp_int(PC_tree_t tree, std::string const& str)
{
    long data = 0;
    if (PC_status_t s = PC_int(PC_get(tree, "%s", str.c_str()), &data); s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

int PCpp_len(PC_tree_t tree, std::string const& str)
{
    int data = 0;
    if (PC_status_t s = PC_len(PC_get(tree, "%s", str.c_str()), &data); s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

bool PCpp_bool(PC_tree_t tree, std::string const& str)
{
    int data = 0;
    if (PC_status_t s = PC_bool(PC_get(tree, "%s", str.c_str()), &data); s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

std::string PCpp_string(PC_tree_t tree, std::string const& str)
{
    char* c_str = nullptr;
    if (PC_status_t s = PC_string(PC_get(tree, "%s", str.c_str()), &c_str); s != PC_OK) {
        if (c_str != nullptr) {
            std::free(c_str);
        }
        throw std::runtime_error(PC_errmsg());
    }
    std::string const str_out(c_str);
    std::free(c_str);
    return str_out;
}

PC_tree_t PCpp_get(PC_tree_t tree, std::string const& str)
{
    return PC_get(tree, "%s", str.c_str());
}
