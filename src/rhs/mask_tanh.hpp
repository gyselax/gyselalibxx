#pragma once

#include <geometry.hpp>

enum class MaskType { Normal, Inverted };

DFieldX mask_tanh(
        IDomainX const& gridx,
        double const extent,
        double const stiffness,
        MaskType type,
        bool const normalized);
