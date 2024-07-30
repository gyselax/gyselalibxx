#pragma once
#include <cmath>

#include <ddc/ddc.hpp>

#include <sll/math_tools.hpp>

#include "ddc_aliases.hpp"

/**
 * A wrapper around a cosine Evaluator that can be used in tests.
 */
struct CosineEvaluator
{
    /**
     * @brief An analytical evaluator to be used for exact comparisons in tests.
     */
    template <class Grid1D>
    class Evaluator
    {
    public:
        /// The grid dimension.
        using Dim = Grid1D;

    private:
        static inline constexpr double s_2_pi = 2. * M_PI;

    private:
        double m_c0;

        double m_c1;

    public:
        /**
         * @brief A constructor taking the range over which the evaluator will be applied.
         * This sets the equation to:
         * @f$ cos(2 \pi x) @f$
         *
         * @param[in] idx_range The range of the grid over which the evaluator will be applied.
         */
        template <class IdxRange>
        Evaluator(IdxRange idx_range) : m_c0(1.0)
                                      , m_c1(0.0)
        {
        }

        /**
         * @brief A constructor parametrising the equation.
         * This sets the equation to:
         * @f$ cos(2 \pi (c_0 x + c_1)) @f$
         *
         * @param[in] c0 The first parameter of the equation.
         * @param[in] c1 The second parameter of the equation.
         */
        Evaluator(double c0, double c1) : m_c0(c0), m_c1(c1) {}

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
            return ipow(s_2_pi * m_c0, diff);
        }

    private:
        double eval(double const x, int const derivative) const noexcept
        {
            return ipow(s_2_pi * m_c0, derivative)
                   * std::cos(M_PI_2 * derivative + s_2_pi * (m_c0 * x + m_c1));
        }
    };
};
