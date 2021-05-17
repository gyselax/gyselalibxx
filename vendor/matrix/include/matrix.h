#pragma once

#include <memory>
#include <string>

#include <view.h>

class Matrix
{
public:
    Matrix(int mat_size) : n(mat_size) {}
    virtual ~Matrix() = default;
    virtual double get_element(int i, int j) const = 0;
    virtual void set_element(int i, int j, double aij) = 0;
    virtual void factorize();
    virtual void solve_inplace(mdspan_1d& b) const;
    virtual void solve_transpose_inplace(mdspan_1d& b) const;
    virtual void solve_inplace_matrix(mdspan_2d& bx) const;
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
    virtual int solve_inplace_method(const char transpose, double* b, int nrows, int ncols)
            const = 0;
    const int n; // matrix size
};

std::ostream& operator<<(std::ostream& o, const Matrix& m);
