#pragma once

#include <geometry.hpp>

enum class MaskType { Normal, Inverted };

DFieldX mask_tanh(
        IDomainX const& gridx,
        double extent,
        double stiffness,
        MaskType const type,
        bool normalized);
