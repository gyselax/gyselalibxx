#include "exact_splitting_solver.h"
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include <iostream>
#include <cmath>

using Matrix3 = Eigen::Matrix3d;
using Vector3 = Eigen::Vector3d;
using Vector3i = Eigen::Vector3i;
using RowVector3 = Eigen::RowVector3d;

// Kronecker product for 3x1 vector v and 1x3 row vector w -> 3x3 matrix
static Matrix3 kron(const Vector3& v, const RowVector3& w)
{
    Matrix3 result;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            result(i, j) = v(i) * w(j);
    return result;
}

// Exponential of 3x3 antisymmetric matrix using Rodrigues formula
static Matrix3 expAntiSym(const Vector3& Bvec)
{
    double theta = std::sqrt(Bvec[0] * Bvec[0] + Bvec[1] * Bvec[1] + Bvec[2] * Bvec[2]);
    
    Matrix3 Bmat;
    Bmat <<      0,  Bvec(2), -Bvec(1),
            -Bvec(2),      0,  Bvec(0),
             Bvec(1), -Bvec(0),      0;
    if (theta < 1e-14)
        return Matrix3::Identity() + Bmat;
    return Matrix3::Identity()
         + (std::sin(theta) / theta) * Bmat
         + ((1.0 - std::cos(theta)) / (theta * theta)) * (Bmat * Bmat);
}

// Exact splitting solver
void exact_splitting_solver(
    const std::array<double, 3>& magnetic_field,
    double t,
    std::array<double, 3>& yl,
    std::array<double, 3>& y2,
    std::array<double, 3>& y3,
    std::array<double, 3>& yr,
    const std::array<int, 3>& order_arr)
{
    // Convert arrays to Eigen
    Vector3 Bvec(magnetic_field[0], magnetic_field[1], magnetic_field[2]);
    Vector3i order(order_arr[0], order_arr[1], order_arr[2]);

    // Construct 3x3 antisymmetric matrix B
    Matrix3 B = Matrix3::Zero();
    B(0,1) =  Bvec(2); B(0,2) = -Bvec(1);
    B(1,0) = -Bvec(2); B(1,2) =  Bvec(0);
    B(2,0) =  Bvec(1); B(2,1) = -Bvec(0);

    Matrix3 I = Matrix3::Identity();

    // Initial vectors
    Vector3 ym = B.col(order(0));
    Vector3 y2_e = B.col(order(1));
    Vector3 y3_e = B.col(order(2));

    // --- Compute yr like Python version ---
    Matrix3 ym_kron = kron(ym, I.row(order(0)));
    Matrix3 y2_kron = kron(y2_e, I.row(order(1)));
    Matrix3 y3_kron = kron(y3_e, I.row(order(2)));

    Vector3 yr_e = -0.5 * ( (ym_kron * B - B * ym_kron + y2_kron * y3_kron - y3_kron * y2_kron).diagonal() );

    yr_e(order(0)) = 0.0;
    yr_e(order(1)) /= -B(order(0), order(1));
    yr_e(order(2)) /= -B(order(0), order(2));

    Vector3 yl_e;

    constexpr int Num = 500;
    //Eigen::VectorXd error_vec(Num);

    // Iterative exact splitting
    for (int i = 0; i < Num; ++i) {
        // Construct P
        Matrix3 P =
            (I + t * kron(ym - yr_e, I.row(order(0)))) *
            (I + t * kron(y2_e,      I.row(order(1)))) *
            (I + t * kron(y3_e,      I.row(order(2)))) *
            (I + t * kron(yr_e,      I.row(order(0))));

        // Matrix logarithm
        Matrix3 g = (1.0 / t) * P.log();

        Vector3 vec = Vector3::Zero();
        vec(order(1)) = g(order(1), order(1)) / (-B(order(0), order(1)));
        vec(order(2)) = g(order(2), order(2)) / (-B(order(0), order(2)));

        Matrix3 II_tau = B * kron(vec, I.row(order(0))) - kron(vec, I.row(order(0))) * B;
        Matrix3 diff = g - II_tau;

        Vector3 ymi = diff.col(order(0));
        Vector3 y2i = diff.col(order(1));
        Vector3 y3i = diff.col(order(2));

        ym   += B.col(order(0)) - ymi;
        y2_e += B.col(order(1)) - y2i;
        y3_e += B.col(order(2)) - y3i;
        yr_e -= 1.0/t * vec;
        yl_e  = ym - yr_e;

        // Optional: track error per iteration
        Matrix3 expB = expAntiSym(Bvec * t);
        Matrix3 approx =
            (I + t * kron(yl_e, I.row(order(0)))) *
            (I + t * kron(y2_e, I.row(order(1)))) *
            (I + t * kron(y3_e, I.row(order(2)))) *
            (I + t * kron(yr_e, I.row(order(0))));
        //error_vec(i) = (expB - approx).norm();
        double error_iter = (expB - approx).norm();
        if (error_iter < 1e-14) {
        //std::cout << "Converged at iteration " << i << ", error = " << error_iter << std::endl;
        break;
        }
    }

    // Copy results back to std::array
    for (int i = 0; i < 3; ++i) {
        yl[i] = yl_e(i);
        y2[i] = y2_e(i);
        y3[i] = y3_e(i);
        yr[i] = yr_e(i);
    }
}
