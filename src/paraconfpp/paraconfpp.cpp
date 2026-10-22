// SPDX-License-Identifier: MIT

#include <cstdlib>
#include <stdexcept>
#include <string_view>

#include <paraconf.h>

#include "paraconfpp.hpp"

double PCpp_double(PC_tree_t tree, std::string_view const sv)
{
    double data = 0.;
    if (PC_status_t s = PC_double(PCpp_get(tree, sv), &data); s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

long PCpp_int(PC_tree_t tree, std::string_view const sv)
{
    long data = 0;
    if (PC_status_t s = PC_int(PCpp_get(tree, sv), &data); s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

int PCpp_len(PC_tree_t tree, std::string_view const sv)
{
    int data = 0;
    if (PC_status_t s = PC_len(PCpp_get(tree, sv), &data); s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

bool PCpp_bool(PC_tree_t tree, std::string_view const sv)
{
    int data = 0;
    if (PC_status_t s = PC_bool(PCpp_get(tree, sv), &data); s != PC_OK) {
        throw std::runtime_error(PC_errmsg());
    }
    return data;
}

std::string PCpp_string(PC_tree_t tree, std::string_view const sv)
{
    char* c_str = nullptr;
    if (PC_status_t s = PC_string(PCpp_get(tree, sv), &c_str); s != PC_OK) {
        if (c_str != nullptr) {
            std::free(c_str);
        }
        throw std::runtime_error(PC_errmsg());
    }
    std::string const str_out(c_str);
    std::free(c_str);
    return str_out;
}

PC_tree_t PCpp_get(PC_tree_t tree, std::string_view const sv)
{
    int const size = sv.size();
    char const* const ptr = sv.data();
    return PC_get(tree, "%.*s", size, ptr);
}
