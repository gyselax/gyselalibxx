// SPDX-License-Identifier: MIT

#pragma once

#include <string_view>

#include <paraconf.h>

/**
 * @brief Retrieves a configuration value as a double.
 *
 * @param tree The configuration tree.
 * @param sv The key identifying the configuration entry.
 *
 * @return The configuration value as a double.
 *
 * @throws std::runtime_error If the value cannot be retrieved or converted.
 */
double PCpp_double(PC_tree_t tree, std::string_view sv);

/**
 * @brief Retrieves a configuration value as an integer.
 *
 * @param tree The configuration tree.
 * @param sv The key identifying the configuration entry.
 *
 * @return The configuration value as a long integer.
 *
 * @throws std::runtime_error If the value cannot be retrieved or converted.
 */
long PCpp_int(PC_tree_t tree, std::string_view sv);

/**
 * @brief Retrieves the length of a configuration entry.
 *
 * @param tree The configuration tree.
 * @param sv The key identifying the configuration entry.
 *
 * @return The length of the configuration entry.
 *
 * @throws std::runtime_error If the length cannot be retrieved.
 */
int PCpp_len(PC_tree_t tree, std::string_view sv);

/**
 * @brief Retrieves a configuration value as a boolean.
 *
 * @param tree The configuration tree.
 * @param sv The key identifying the configuration entry.
 *
 * @return The configuration value as a boolean.
 *
 * @throws std::runtime_error If the value cannot be retrieved or converted.
 */
bool PCpp_bool(PC_tree_t tree, std::string_view sv);

/**
 * @brief Retrieves a configuration value as a string.
 *
 * @param tree The configuration tree.
 * @param sv The key identifying the configuration entry.
 *
 * @return The configuration value as a string.
 *
 * @throws std::runtime_error If the value cannot be retrieved or converted.
 */
std::string PCpp_string(PC_tree_t tree, std::string_view sv);

/**
 * @brief Retrieves a configuration entry.
 *
 * @param tree The configuration tree.
 * @param sv The key identifying the configuration entry.
 *
 * @return The configuration entry, or an invalid handle if the entry
 *         cannot be found.
 */
PC_tree_t PCpp_get(PC_tree_t tree, std::string_view sv);
