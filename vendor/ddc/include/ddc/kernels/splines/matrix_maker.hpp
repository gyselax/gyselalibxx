#pragma once
#include <memory>

#include "matrix_sparse.hpp"


namespace ddc::detail {
class MatrixMaker
{
public:
    template <typename ExecSpace>
    static std::unique_ptr<Matrix> make_new_sparse(
            int const n,
            std::optional<int> cols_per_par_chunk = std::nullopt,
            std::optional<int> par_chunks_per_seq_chunk = std::nullopt,
            std::optional<unsigned int> preconditionner_max_block_size = std::nullopt)
    {
        return std::make_unique<Matrix_Sparse<ExecSpace>>(
                n,
                cols_per_par_chunk,
                par_chunks_per_seq_chunk,
                preconditionner_max_block_size);
    }
};
} // namespace ddc::detail
