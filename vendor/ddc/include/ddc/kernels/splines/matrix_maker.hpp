// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#pragma once

#include <memory>
#include <optional>

#include "matrix_sparse.hpp"

namespace ddc::detail {

class MatrixMaker
{
public:
    template <typename ExecSpace>
    static std::unique_ptr<Matrix> make_new_sparse(
            int const n,
            std::optional<int> cols_per_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
    {
        return std::make_unique<
                Matrix_Sparse<ExecSpace>>(n, cols_per_chunk, preconditionner_max_block_size);
    }
};

} // namespace ddc::detail
