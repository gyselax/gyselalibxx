#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"

/**
 * A wrapper around a 2D Evaluator that can be used in tests.
 */
struct Evaluator2D
{
    /**
     * @brief An analytical 2D evaluator combining 2 1D evaluators, to be used for exact comparisons in tests.
     */
    template <class Eval1, class Eval2>
    class Evaluator
    {
    private:
        Eval1 eval_func1;
        Eval2 eval_func2;

    public:
        /**
         * @brief A constructor taking the range over which the evaluator will be applied.
         * The polynomial is initialised with random values between 0.0 and 1.0
         *
         * @param[in] idx_range The range of the grid over which the evaluator will be applied.
         */
        template <class IdxRange>
        Evaluator(IdxRange idx_range)
            : eval_func1(ddc::select<typename Eval1::Dim>(idx_range))
            , eval_func2(ddc::select<typename Eval2::Dim>(idx_range))
        {
        }

        /**
         * Evaluate the equation at (x,y)
         * @param[in] x The x-position where the equation is evaluated.
         * @param[in] y The y-position where the equation is evaluated.
         * @returns The result of the equation.
         */
        double operator()(double const x, double const y) const noexcept
        {
            return eval_func1(x) * eval_func2(y);
        }

        /**
         * Evaluate the equation at @f$(x_1,x_2)@f$.
         * @param[in] x The position where the equation is evaluated.
         * @returns The result of the equation.
         */
        template <class Grid1, class Grid2>
        double operator()(ddc::Coordinate<Grid1, Grid2> const x) const noexcept
        {
            return eval_func1(ddc::get<Grid1>(x)) * eval_func2(ddc::get<Grid2>(x));
        }

        /**
         * Evaluate the equation at the positions on which chunk is defined.
         * @param[out] chunk The result of the equation.
         */
        template <class Grid1, class Grid2>
        void operator()(ddc::ChunkSpan<double, ddc::DiscreteDomain<Grid1, Grid2>> chunk) const
        {
            auto const& idx_range = chunk.idx_range();

            ddc::for_each(idx_range, [=](ddc::DiscreteElement<Grid1, Grid2> const i) {
                chunk(i) = eval_func1(ddc::coordinate(ddc::select<Grid1>(i)))
                           * eval_func2(ddc::coordinate(ddc::select<Grid2>(i)));
            });
        }

        /**
         * Evaluate the derivative of the equation at (x,y)
         * @param[in] x The x-position where the equation is evaluated.
         * @param[in] y The y-position where the equation is evaluated.
         * @param[in] derivative_x The order of the x-derivative.
         * @param[in] derivative_y The order of the x-derivative.
         * @returns The result of the equation.
         */
        double deriv(double const x, double const y, int const derivative_x, int const derivative_y)
                const noexcept
        {
            return eval_func1.deriv(x, derivative_x) * eval_func2.deriv(y, derivative_y);
        }

        /**
         * Evaluate the derivative of the equation at @f$(x_1,x_2)@f$.
         * @param[in] x The position where the equation is evaluated.
         * @param[in] derivative_x The order of the x-derivative.
         * @param[in] derivative_y The order of the x-derivative.
         * @returns The result of the equation.
         */
        template <class Grid1, class Grid2>
        double deriv(
                ddc::Coordinate<Grid1, Grid2> const x,
                int const derivative_x,
                int const derivative_y) const noexcept
        {
            return eval_func1.deriv(ddc::get<Grid1>(x), derivative_x)
                   * eval_func2.deriv(ddc::get<Grid2>(x), derivative_y);
        }

        /**
         * Evaluate the derivative of the equation at the positions on which chunk is defined.
         * @param[out] chunk The result of the equation.
         * @param[in] derivative_x The order of the x-derivative.
         * @param[in] derivative_y The order of the x-derivative.
         */
        template <class Grid1, class Grid2>
        void deriv(
                ddc::ChunkSpan<double, ddc::DiscreteDomain<Grid1, Grid2>> chunk,
                int const derivative_x,
                int const derivative_y) const
        {
            auto const& idx_range = chunk.idx_range();

            ddc::for_each(idx_range, [=](ddc::DiscreteElement<Grid1, Grid2> const i) {
                chunk(i) = eval_func1.deriv(ddc::coordinate(ddc::select<Grid1>(i)), derivative_x)
                           * eval_func2.deriv(ddc::coordinate(ddc::select<Grid2>(i)), derivative_y);
            });
        }

        /**
         * The maximum norm of this equation. Used for normalisation.
         * @param[in] diff1 The x-derivative whose max norm should be returned.
         * @param[in] diff2 The y-derivative whose max norm should be returned.
         * @returns The maximum value of the infinite norm.
         */
        double max_norm(int diff1 = 0, int diff2 = 0) const
        {
            return eval_func1.max_norm(diff1) * eval_func2.max_norm(diff2);
        }
    };
};
