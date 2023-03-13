#pragma once

#include <sll/polar_spline.hpp>

template <class PolarBSplinesType>
class PolarSplineEvaluator
{
private:
    // Tags to determine what to evaluate
    struct eval_type
    {
    };

    struct eval_deriv_r_type
    {
    };

    struct eval_deriv_p_type
    {
    };

    struct eval_deriv_r_p_type
    {
    };

public:
    using bsplines_type = PolarBSplinesType;
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplinesP = typename PolarBSplinesType::BSplinesP_tag;
    using DimR = typename BSplinesR::tag_type;
    using DimP = typename BSplinesP::tag_type;

public:
    static int constexpr continuity = PolarBSplinesType::continuity;

private:
    PolarSplineBoundaryValue2D<PolarBSplinesType> const& m_outer_bc;

public:
    PolarSplineEvaluator() = delete;

    explicit PolarSplineEvaluator(PolarSplineBoundaryValue2D<PolarBSplinesType> const& outer_bc)
        : m_outer_bc(outer_bc)
    {
    }

    PolarSplineEvaluator(PolarSplineEvaluator const& x) = default;

    PolarSplineEvaluator(PolarSplineEvaluator&& x) = default;

    ~PolarSplineEvaluator() = default;

    PolarSplineEvaluator& operator=(PolarSplineEvaluator const& x) = default;

    PolarSplineEvaluator& operator=(PolarSplineEvaluator&& x) = default;

    double operator()(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        return eval(coord_eval, spline_coef);
    }

    template <class Domain>
    void operator()(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<DimR, DimP> const, Domain> const coords_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval(coords_eval(i), spline_coef);
        });
    }

    double deriv_dim_1(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_r_type());
    }

    double deriv_dim_2(
            ddc::Coordinate<DimR, DimP> coord_eval,
            double const coord_eval1,
            double const coord_eval2,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_p_type());
    }

    double deriv_1_and_2(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_r_p_type());
    }

    template <class Domain>
    void deriv_dim_1(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<DimR, DimP> const, Domain> const coords_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_r_type());
        });
    }

    template <class Domain>
    void deriv_dim_2(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<DimR, DimP> const, Domain> const coords_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_p_type());
        });
    }

    template <class Domain>
    void deriv_dim_1_and_2(
            ddc::ChunkSpan<double, Domain> const spline_eval,
            ddc::ChunkSpan<ddc::Coordinate<DimR, DimP> const, Domain> const coords_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        ddc::for_each(coords_eval.domain(), [=](auto i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_r_p_type());
        });
    }

    template <class Mapping>
    double integrate(PolarSplineView<PolarBSplinesType> const spline_coef, Mapping const mapping)
            const
    {
        int constexpr nr = ddc::discrete_space<BSplinesR>().ncells() + BSplinesR::degree() - 2;
        int constexpr np = ddc::discrete_space<BSplinesP>().ncells() + BSplinesP::degree();
        std::array<double, PolarBSplinesType::eval_size()> singular_values;
        DSpan1D singular_vals(singular_values.data(), PolarBSplinesType::n_singular_basis());
        std::array<double, nr * np> values;
        DSpan2D vals(values.data(), nr, np);

        ddc::discrete_space<PolarBSplinesType>().integrals(singular_vals, vals);

        double y = 0.;
        ddc::for_each(
                spline_coef.singular_spline_coef.domain(),
                [=](ddc::DiscreteElement<PolarBSplinesType> const i) {
                    y += spline_coef.singular_spline_coef(i) * singular_vals(i)
                         * mapping.determinant(i);
                });
        ddc::for_each(
                spline_coef.spline_coef.domain(),
                [=](ddc::DiscreteElement<BSplinesR, BSplinesP> const i) {
                    y += spline_coef.spline_coef(i)
                         * vals(ddc::select<BSplinesR>(i), ddc::select<BSplinesP>(i))
                         * mapping.determinant(i);
                });
        return y;
    }

private:
    double eval(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef) const
    {
        const double coord_eval1 = ddc::get<DimR>(coord_eval);
        double coord_eval2 = ddc::get<DimP>(coord_eval);
        if (coord_eval1 > ddc::discrete_space<BSplinesR>().rmax()) {
            return m_outer_bc(coord_eval1, coord_eval2, spline_coef);
        }
        if (coord_eval2 < ddc::discrete_space<BSplinesP>().rmin()
            || coord_eval2 > ddc::discrete_space<BSplinesP>().rmax()) {
            coord_eval2 -= std::floor(
                                   (coord_eval2 - ddc::discrete_space<BSplinesP>().rmin())
                                   / ddc::discrete_space<BSplinesP>().length())
                           * ddc::discrete_space<BSplinesP>().length();
        }
        return eval_no_bc(
                ddc::Coordinate<DimR, DimP>(coord_eval1, coord_eval2),
                spline_coef,
                eval_type());
    }

    template <class EvalType>
    double eval_no_bc(
            ddc::Coordinate<DimR, DimP> coord_eval,
            PolarSplineView<PolarBSplinesType> const spline_coef,
            EvalType const) const
    {
        static_assert(
                std::is_same_v<
                        EvalType,
                        eval_type> || std::is_same_v<EvalType, eval_deriv_r_type> || std::is_same_v<EvalType, eval_deriv_p_type> || std::is_same_v<EvalType, eval_deriv_r_p_type>);

        std::array<double, PolarBSplinesType::n_singular_basis()> singular_data;
        DSpan1D singular_vals(singular_data.data(), PolarBSplinesType::n_singular_basis());
        std::array<double, (BSplinesR::degree() + 1) * (BSplinesP::degree() + 1)> data;
        DSpan2D vals(data.data(), BSplinesR::degree() + 1, BSplinesP::degree() + 1);

        ddc::DiscreteElement<BSplinesR, BSplinesP> jmin;

        if constexpr (std::is_same_v<EvalType, eval_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_basis(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_r_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_r(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_p_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_p(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_r_p_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_r_and_p(singular_vals, vals, coord_eval);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < PolarBSplinesType::n_singular_basis(); ++i) {
            y += spline_coef.singular_spline_coef(ddc::DiscreteElement<PolarBSplinesType>(i))
                 * singular_vals(i);
        }
        ddc::DiscreteElement<BSplinesR> jmin_r = ddc::select<BSplinesR>(jmin);
        ddc::DiscreteElement<BSplinesP> jmin_p = ddc::select<BSplinesP>(jmin);
        int nr = BSplinesR::degree() + 1;
        if (jmin_r.uid() < continuity + 1) {
            nr = nr - (continuity + 1 - jmin_r.uid());
            jmin_r = ddc::DiscreteElement<BSplinesR>(continuity + 1);
        }
        for (int i = 0; i < nr; ++i) {
            for (std::size_t j = 0; j < BSplinesP::degree() + 1; ++j) {
                y += spline_coef.spline_coef(jmin_r + i, jmin_p + j) * vals(i, j);
            }
        }
        return y;
    }
};
