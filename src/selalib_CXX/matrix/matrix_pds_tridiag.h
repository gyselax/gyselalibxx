#ifndef MATRIX_PDS_BANDED_H
#define MATRIX_PDS_BANDED_H
#include "matrix_banded.h"

class Matrix_PDS_Tridiag: public Matrix {
    /*
     * Represents a real symmetric positive definite matrix
     * stored in a block format
     * */
    public:
        Matrix_PDS_Tridiag(int n);
        virtual ~Matrix_PDS_Tridiag();
        virtual double get_element(int i, int j) const override;
        virtual void set_element(int i, int j, double a_ij) override;
    protected:
        virtual int factorize_method() override;
        virtual int solve_inplace_method(const char transpose, double* b, int nrows, int ncols) const override;
        double* d; //diagonal
        double* l; // lower diagonal
};

#endif // MATRIX_SYMMETRIC_BANDED_H

