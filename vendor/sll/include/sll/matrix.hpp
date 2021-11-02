#pragma once

#include <iosfwd>
#include <memory>

#include "sll/view.hpp"

class Matrix
{
public:
    Matrix(int mat_size) : n(mat_size) {}
    virtual ~Matrix() = default;
    virtual double get_element(int i, int j) const = 0;
    virtual void set_element(int i, int j, double aij) = 0;
    virtual void factorize();
    virtual DSpan1D solve_inplace(DSpan1D b) const;
    virtual DSpan1D solve_transpose_inplace(DSpan1D b) const;
    virtual DSpan2D solve_inplace_matrix(DSpan2D bx) const;
    int get_size() const
    {
        return n;
    }
    static std::unique_ptr<Matrix> make_new_banded(int n, int kl, int ku, bool pds);
    static std::unique_ptr<Matrix> make_new_periodic_banded(int n, int kl, int ku, bool pds);
    static std::unique_ptr<Matrix> make_new_block_with_banded_region(
            int n,
            int kl,
            int ku,
            bool pds,
            int block1_size,
            int block2_size = 0);

protected:
    virtual int factorize_method() = 0;
    virtual int solve_inplace_method(double* b, char transpose, int nrows, int ncols) const = 0;
    int const n; // matrix size
};

std::ostream& operator<<(std::ostream& o, Matrix const& m);
