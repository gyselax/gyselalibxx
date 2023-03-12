// SPDX-License-Identifier: MIT

#pragma once

template <class Tag>
struct Fourier
{
    using base_tag_type = Tag;
    static bool constexpr PERIODIC = Tag::PERIODIC;
};
