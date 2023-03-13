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
        template <class Domain>
        Evaluator(Domain domain)
            : eval_func1(ddc::select<typename Eval1::Dim>(domain))
            , eval_func2(ddc::select<typename Eval2::Dim>(domain))
        {
        }

        double operator()(double const x, double const y) const noexcept
        {
            return eval_func1(x) * eval_func2(y);
        }

        template <class DDim1, class DDim2>
        double operator()(ddc::Coordinate<DDim1, DDim2> const x) const noexcept
        {
            return eval_func1(ddc::get<DDim1>(x)) * eval_func2(ddc::get<DDim2>(x));
        }

        template <class DDim1, class DDim2>
        void operator()(ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim1, DDim2>> chunk) const
        {
            auto const& domain = chunk.domain();

            ddc::for_each(domain, [=](ddc::DiscreteElement<DDim1, DDim2> const i) {
                chunk(i) = eval_func1(ddc::coordinate(ddc::select<DDim1>(i)))
                           * eval_func2(ddc::coordinate(ddc::select<DDim2>(i)));
            });
        }

        double deriv(double const x, double const y, int const derivative_x, int const derivative_y)
                const noexcept
        {
            return eval_func1.deriv(x, derivative_x) * eval_func2.deriv(y, derivative_y);
        }

        template <class DDim1, class DDim2>
        double deriv(
                ddc::Coordinate<DDim1, DDim2> const x,
                int const derivative_x,
                int const derivative_y) const noexcept
        {
            return eval_func1.deriv(ddc::get<DDim1>(x), derivative_x)
                   * eval_func2.deriv(ddc::get<DDim2>(x), derivative_y);
        }

        template <class DDim1, class DDim2>
        void deriv(
                ddc::ChunkSpan<double, ddc::DiscreteDomain<DDim1, DDim2>> chunk,
                int const derivative_x,
                int const derivative_y) const
        {
            auto const& domain = chunk.domain();

            ddc::for_each(domain, [=](ddc::DiscreteElement<DDim1, DDim2> const i) {
                chunk(i) = eval_func1.deriv(ddc::coordinate(ddc::select<DDim1>(i)), derivative_x)
                           * eval_func2.deriv(ddc::coordinate(ddc::select<DDim2>(i)), derivative_y);
            });
        }

        double max_norm(int diff1 = 0, int diff2 = 0) const
        {
            return eval_func1.max_norm(diff1) * eval_func2.max_norm(diff2);
        }
    };
};
