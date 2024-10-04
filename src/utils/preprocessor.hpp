// SPDX-License-Identifier: MIT

#pragma once

////////////////////////////////////////
/// GSLX_UNUSED
///
/// Prevent a compiler warning for an unused variable.
///
/// Example usage:
///
///     int a;
///     GSLX_UNUSED(a);
///     return;
///
////////////////////////////////////////

#define GSLX_UNUSED(an_expression) static_cast<void>(an_expression)


////////////////////////////////////////
/// GSLX_UTILITY_STRINGIFY
///
/// This macro will convert the text as argument to a string literal.
///
/// Example usage:
///
///     GSLX_UTILITY_STRINGIFY(THIS IS A TEST);
///
/// Expands to :
///
///     "THIS IS A TEST";
///
////////////////////////////////////////

#define _GSLX_UTILITY_STRINGIFY(an_expression) #an_expression
#define GSLX_UTILITY_STRINGIFY(an_expression) _GSLX_UTILITY_STRINGIFY(an_expression)
