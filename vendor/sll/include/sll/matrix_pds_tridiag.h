#ifndef MATRIX_PDS_BANDED_H
#define MATRIX_PDS_BANDED_H

#include <memory>

#include "sll/matrix.h"

class Matrix_PDS_Tridiag : public Matrix
{
    /*
     * Represents a real symmetric positive definite matrix
     * stored in a block format
     * */
public:
    Matrix_PDS_Tridiag(int n);
    virtual double get_element(int i, int j) const override;
    virtual void set_element(int i, int j, double a_ij) override;

protected:
    virtual int factorize_method() override;
    virtual int solve_inplace_method(const char transpose, double* b, int nrows, int ncols)
            const override;
    std::unique_ptr<double[]> d; // diagonal
    std::unique_ptr<double[]> l; // lower diagonal
};

#endif // MATRIX_SYMMETRIC_BANDED_H
