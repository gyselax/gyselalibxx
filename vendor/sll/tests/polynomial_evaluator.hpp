#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <random>

#include <ddc/ddc.hpp>

#include <sll/math_tools.hpp>

#include "ddc_aliases.hpp"

/**
 * A wrapper around a polynomial Evaluator that can be used in tests.
 */
struct PolynomialEvaluator
{
    /**
     * @brief An analytical evaluator to be used for exact comparisons in tests.
     */
    template <class Grid1D, std::size_t Degree>
    class Evaluator
    {
    public:
        /// The grid dimension.
        using Dim = Grid1D;

    private:
        std::array<double, Degree + 1> m_coeffs;
        int const m_degree;
        double const m_xN;

    public:
        /**
         * @brief A constructor taking the range over which the evaluator will be applied.
         * The polynomial is initialised with random values between 0.0 and 1.0
         *
         * @param[in] idx_range The range of the grid over which the evaluator will be applied.
         */
        template <class IdxRange>
        Evaluator(IdxRange idx_range)
            : m_degree(Degree)
            , m_xN(std::max(std::abs(rmin(idx_range)), std::abs(rmax(idx_range))))
        {
            for (int i(0); i < m_degree + 1; ++i) {
                m_coeffs[i] = double(rand() % 100) / 100.0;
            }
        }

        /**
         * Evaluate the equation at x
         * @param[in] x The position where the equation is evaluated.
         * @returns The result of the equation.
         */
        double operator()(double const x) const noexcept
        {
            return eval(x, 0);
        }

        /**
         * Evaluate the equation at the positions on which chunk is defined.
         * @param[out] chunk The result of the equation.
         */
        void operator()(ddc::ChunkSpan<double, ddc::DiscreteDomain<Grid1D>> chunk) const
        {
            auto const& idx_range = chunk.idx_range();

            for (ddc::DiscreteElement<Grid1D> const i : idx_range) {
                chunk(i) = eval(ddc::coordinate(i), 0);
            }
        }

        /**
         * Evaluate the derivative of the equation at x
         * @param[in] x The position where the equation is evaluated.
         * @param[in] derivative The order of the derivative.
         * @returns The result of the equation.
         */
        double deriv(double const x, int const derivative) const noexcept
        {
            return eval(x, derivative);
        }

        /**
         * Evaluate the derivative of the equation at the positions on which chunk is defined.
         * @param[out] chunk The result of the equation.
         * @param[in] derivative The order of the derivative.
         */
        void deriv(ddc::ChunkSpan<double, ddc::DiscreteDomain<Grid1D>> chunk, int const derivative)
                const
        {
            auto const& idx_range = chunk.idx_range();

            for (ddc::DiscreteElement<Grid1D> const i : idx_range) {
                chunk(i) = eval(ddc::coordinate(i), derivative);
            }
        }

        /**
         * The maximum norm of this equation. Used for normalisation.
         * @param[in] diff The derivative whose max norm should be returned.
         * @returns The maximum value of the infinite norm.
         */
        double max_norm(int diff = 0) const
        {
            return std::abs(deriv(m_xN, diff));
        }

    private:
        double eval(double const x, int const derivative) const
        {
            double result(0.0);
            int start = derivative < 0 ? 0 : derivative;
            for (int i(start); i < m_degree + 1; ++i) {
                double v = double(falling_factorial(i, derivative)) * std::pow(x, i - derivative);
                result += m_coeffs[i] * v;
            }
            return result;
        }

        double falling_factorial(int i, int d) const
        {
            double c = 1.0;
            if (d >= 0) {
                for (int k(0); k < d; ++k) {
                    c *= (i - k);
                }
            } else {
                for (int k(-1); k > d - 1; --k) {
                    c /= (i - k);
                }
            }
            return c;
        }
    };
};
