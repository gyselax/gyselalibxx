#pragma once

#include <ddc/ddc.hpp>

struct Evaluator2D
{
    template <class Eval1, class Eval2>
    class Evaluator
    {
    private:
        Eval1 eval_func1;
        Eval2 eval_func2;

    public:
        Evaluator() = default;

        double operator()(double const x, double const y) const noexcept
        {
            return eval_func1(x) * eval_func2(y);
        }

        template <class DDim1, class DDim2>
        double operator()(Coordinate<DDim1, DDim2> const x) const noexcept
        {
            return eval_func1(get<DDim1>(x)) * eval_func2(get<DDim2>(x));
        }

        template <class DDim1, class DDim2>
        void operator()(ChunkSpan<double, DiscreteDomain<DDim1, DDim2>> chunk) const
        {
            auto const& domain = chunk.domain();

            for_each(domain, [=](DiscreteElement<DDim1, DDim2> const i) {
                chunk(i) = eval_func1(coordinate(select<DDim1>(i)))
                           * eval_func2(coordinate(select<DDim2>(i)));
            });
        }

        double deriv(double const x, double const y, int const derivative_x, int const derivative_y)
                const noexcept
        {
            return eval_func1.deriv(x, derivative_x) * eval_func2.deriv(y, derivative_y);
        }

        template <class DDim1, class DDim2>
        double deriv(
                Coordinate<DDim1, DDim2> const x,
                int const derivative_x,
                int const derivative_y) const noexcept
        {
            return eval_func1.deriv(get<DDim1>(x), derivative_x)
                   * eval_func2.deriv(get<DDim2>(x), derivative_y);
        }

        template <class DDim1, class DDim2>
        void deriv(
                ChunkSpan<double, DiscreteDomain<DDim1, DDim2>> chunk,
                int const derivative_x,
                int const derivative_y) const
        {
            auto const& domain = chunk.domain();

            for_each(domain, [=](DiscreteElement<DDim1, DDim2> const i) {
                chunk(i) = eval_func1.deriv(coordinate(select<DDim1>(i)), derivative_x)
                           * eval_func2.deriv(coordinate(select<DDim2>(i)), derivative_y);
            });
        }
    };
};
