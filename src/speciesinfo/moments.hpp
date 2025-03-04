// SPDX-License-Identifier: MIT

#pragma once

/// @brief Moments discrete dimension to access constant attributes related to fluid moments.
class Moments
{
public:
    /// alias of the discrete dimension
    using discrete_dimension_type = Moments;

public:
    /// @brief Impl object storing attributes in `MemorySpace`.
    template <class Grid1D, class MemorySpace>
    class Impl
    {
        template <class ODDim, class OMemorySpace>
        friend class Impl;

    public:
        /// alias of the discrete dimension
        using discrete_dimension_type = Moments;

        /**
         * @brief Conversion constructor between different memory spaces.
         * @param[in] impl object from `OMemorySpace` that will be used to initialise this object on `MemorySpace`
         */
        template <class OMemorySpace>
        explicit Impl(Impl<Grid1D, OMemorySpace> const& impl)
        {
        }

        /**
         * @brief Main constructor
         */
        Impl() {}
    };
};
