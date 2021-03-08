#ifndef MATRIX_BANDED_H
#define MATRIX_BANDED_H
#include "matrix.h"

class Matrix_Banded: public Matrix {
    public:
        Matrix_Banded(int n, int kl, int ku);
        virtual ~Matrix_Banded();
        virtual double get_element(int i, int j) const override;
        virtual void set_element(int i, int j, double a_ij) override;
    protected:
        virtual int factorize_method() override;
        virtual int solve_inplace_method(const char transpose, double* b, int nrows, int ncols) const override;
        const int kl; // no. of subdiagonals
        const int ku; // no. of superdiagonals
        const int c;   // no. of columns in q
        int* ipiv;    // pivot indices
        double* q;    // banded matrix representation
};

#endif // MATRIX_BANDED_H
