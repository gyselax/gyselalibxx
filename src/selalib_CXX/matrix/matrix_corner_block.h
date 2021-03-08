#ifndef MATRIX_CORNER_BLOCK_H
#define MATRIX_CORNER_BLOCK_H
#include <memory>
#include "matrix.h"
#include "matrix_dense.h"

class Matrix_Corner_Block: public Matrix {
    public:
        Matrix_Corner_Block(int n, int k, std::unique_ptr<Matrix> q);
        virtual ~Matrix_Corner_Block();
        virtual double get_element(int i, int j) const override;
        virtual void set_element(int i, int j, double a_ij) override;
        virtual void factorize() override;
        virtual void solve_inplace(mdspan_1d& bx) const override;
        virtual void solve_transpose_inplace(mdspan_1d& bx) const override;
        virtual void solve_inplace_matrix(mdspan_2d& bx) const override;
    protected:
        Matrix_Corner_Block(int n, int k, std::unique_ptr<Matrix> q, bool is_virtual);
        virtual void allocate_lambda();
        virtual void calculate_delta_to_factorize();
        virtual void solve_lambda_section(mdspan_1d& v, mdspan_1d const& u) const;
        virtual void solve_lambda_section(mdspan_2d& v, mdspan_2d const& u) const;
        virtual void solve_lambda_section_transpose(mdspan_1d& v, mdspan_1d const& u) const;
        const int k;  // small block size
        const int nb; // main block matrix size
        double* Abm_1_gamma_ptr;
        double* lambda_ptr;
        //-------------------------------------
        //
        //    q = | q_block | gamma |
        //        |  lambda | delta |
        //
        //-------------------------------------
        std::unique_ptr<Matrix> q_block;
        Matrix_Dense            delta;
        mdspan_2d               Abm_1_gamma;
        mdspan_2d*              lambda;
    private:
        virtual int factorize_method() override {return 0;}
        virtual int solve_inplace_method(const char transpose, double* b,
                                        int nrows, int ncols) const override
             {return 0;}
};

#endif // MATRIX_CORNER_BLOCK_H
