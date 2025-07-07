// SPDX-License-Identifier: MIT
#pragma once

/**
 * A basic implementation of std::cyl_bessel_j. This implementation
 * is provided for testing as std::cyl_bessel_j is not available on
 * libc++. It is not robust or optimised in any way so it should not
 * be used for production code or for tight error bounds tests.
 */
double cyl_bessel_j(double nu, double x);
