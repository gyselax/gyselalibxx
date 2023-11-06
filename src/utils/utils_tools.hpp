/**
 * @file utils_tools.hpp
 * File Describing useful functions.
 */

#pragma once

#include <ddc/ddc.hpp>


/**
 * @brief Compute the infinity norm.
 *
 * For a given vector @f$ x @f$ , compute
 * @f$|Vert x |Vert_{\infty} = \sup_n |x_n| @f$.
 *
 * @param[in] coord
 *      The given vector.
 *
 * @return A double containing the value of the infinty norm.
 */
template <class... Tags>
double norm_inf(ddc::Coordinate<Tags...> coord)
{
    double result = 0.0;
    ((result = std::max(result, fabs(coord.template get<Tags>()))), ...);
    return result;
};


/**
 * @brief Compute the infinity norm.
 *
 * In case of scalar, the infinity norm
 * returns the scalar.
 *
 * @param[in] coord
 *      The given double.
 *
 * @return A double containing the value of the infinty norm.
 */
inline double norm_inf(double const coord)
{
    return coord;
};
