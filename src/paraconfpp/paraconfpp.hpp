// SPDX-License-Identifier: MIT

#pragma once

#include <string>

#include <paraconf.h>

double PCpp_double(PC_tree_t tree, std::string const& str);

long PCpp_int(PC_tree_t tree, std::string const& str);

int PCpp_len(PC_tree_t tree, std::string const& str);

bool PCpp_bool(PC_tree_t tree, std::string const& str);

std::string PCpp_string(PC_tree_t tree, std::string const& str);

PC_tree_t PCpp_get(PC_tree_t tree, std::string const& str);
