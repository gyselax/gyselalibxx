// SPDX-License-Identifier: MIT
#pragma once

/**
 * @brief A class for describing curvilinear 2D mappings from the logical domain to the physical domain.
 * */
template <class X, class Y, class R, class Theta>
class Curvilinear2DToCartesian
{
public:
    /// @brief Indicate the first physical coordinate.
    using cartesian_tag_x = X;
    /// @brief Indicate the second physical coordinate.
    using cartesian_tag_y = Y;
    /// @brief Indicate the first logical coordinate.
    using curvilinear_tag_r = R;
    /// @brief Indicate the second logical coordinate.
    using curvilinear_tag_theta = Theta;
};
