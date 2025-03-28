

# File polar\_spline\_evaluator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**polar\_splines**](dir_a6779ae02b71d57f488d261458bab1ce.md) **>** [**polar\_spline\_evaluator.hpp**](polar__spline__evaluator_8hpp.md)

[Go to the documentation of this file](polar__spline__evaluator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "polar_spline.hpp"
#include "view.hpp"

template <class PolarBSplinesType, class OuterExtrapolationRule>
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

    struct eval_deriv_theta_type
    {
    };

    struct eval_deriv_r_theta_type
    {
    };

public:
    using bsplines_type = PolarBSplinesType;
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;
    using DimR = typename BSplinesR::continuous_dimension_type;
    using DimTheta = typename BSplinesTheta::continuous_dimension_type;

public:
    static int constexpr continuity = PolarBSplinesType::continuity;

private:
    OuterExtrapolationRule m_outer_bc;

public:
    PolarSplineEvaluator() = delete;

    explicit PolarSplineEvaluator(OuterExtrapolationRule const& outer_bc) : m_outer_bc(outer_bc) {}

    PolarSplineEvaluator(PolarSplineEvaluator const& x) = default;

    PolarSplineEvaluator(PolarSplineEvaluator&& x) = default;

    ~PolarSplineEvaluator() = default;

    PolarSplineEvaluator& operator=(PolarSplineEvaluator const& x) = default;

    PolarSplineEvaluator& operator=(PolarSplineEvaluator&& x) = default;

    double operator()(
            ddc::Coordinate<DimR, DimTheta> coord_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        return eval(coord_eval, spline_coef);
    }

    template <class Domain>
    void operator()(
            DField<Domain, Kokkos::HostSpace> const spline_eval,
            ConstField<ddc::Coordinate<DimR, DimTheta>, Domain, Kokkos::HostSpace> const
                    coords_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::for_each(get_idx_range(coords_eval), [=](IdxEval i) {
            spline_eval(i) = eval(coords_eval(i), spline_coef);
        });
    }

    double deriv_dim_1(
            ddc::Coordinate<DimR, DimTheta> coord_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_r_type());
    }

    double deriv_dim_2(
            ddc::Coordinate<DimR, DimTheta> coord_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_theta_type());
    }


    double deriv_1_and_2(
            ddc::Coordinate<DimR, DimTheta> coord_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, eval_deriv_r_theta_type());
    }

    template <class Domain>
    void deriv_dim_1(
            DField<Domain, Kokkos::HostSpace> const spline_eval,
            ConstField<ddc::Coordinate<DimR, DimTheta>, Domain, Kokkos::HostSpace> const
                    coords_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::for_each(get_idx_range(coords_eval), [=](IdxEval i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_r_type());
        });
    }

    template <class Domain>
    void deriv_dim_2(
            DField<Domain, Kokkos::HostSpace> const spline_eval,
            ConstField<ddc::Coordinate<DimR, DimTheta>, Domain, Kokkos::HostSpace> const
                    coords_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::for_each(get_idx_range(coords_eval), [=](IdxEval i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_theta_type());
        });
    }

    template <class Domain>
    void deriv_dim_1_and_2(
            DField<Domain, Kokkos::HostSpace> const spline_eval,
            ConstField<ddc::Coordinate<DimR, DimTheta>, Domain, Kokkos::HostSpace> const
                    coords_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::for_each(get_idx_range(coords_eval), [=](IdxEval i) {
            spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, eval_deriv_r_theta_type());
        });
    }

    template <class Mapping>
    double integrate(
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef,
            Mapping const mapping) const
    {
        int constexpr nr = ddc::discrete_space<BSplinesR>().ncells() + BSplinesR::degree() - 2;
        int constexpr ntheta
                = ddc::discrete_space<BSplinesTheta>().ncells() + BSplinesTheta::degree();
        std::array<double, PolarBSplinesType::eval_size()> singular_values;
        DSpan1D singular_vals(singular_values.data(), PolarBSplinesType::n_singular_basis());
        std::array<double, nr * ntheta> values;
        DSpan2D vals(values.data(), nr, ntheta);

        ddc::discrete_space<PolarBSplinesType>().integrals(singular_vals, vals);

        double y = 0.;
        ddc::for_each(
                get_idx_range(spline_coef.singular_spline_coef),
                [=](Idx<PolarBSplinesType> const i) {
                    y += spline_coef.singular_spline_coef(i) * singular_vals(i)
                         * mapping.determinant(i);
                });
        ddc::for_each(
                get_idx_range(spline_coef.spline_coef),
                [=](Idx<BSplinesR, BSplinesTheta> const i) {
                    y += spline_coef.spline_coef(i)
                         * vals(ddc::select<BSplinesR>(i), ddc::select<BSplinesTheta>(i))
                         * mapping.determinant(i);
                });
        return y;
    }

private:
    double eval(
            ddc::Coordinate<DimR, DimTheta> coord_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef) const
    {
        const double coord_eval1 = ddc::get<DimR>(coord_eval);
        double coord_eval2 = ddc::get<DimTheta>(coord_eval);
        if (coord_eval1 > ddc::discrete_space<BSplinesR>().rmax()) {
            return m_outer_bc(coord_eval, spline_coef);
        }
        if (coord_eval2 < ddc::discrete_space<BSplinesTheta>().rmin()
            || coord_eval2 > ddc::discrete_space<BSplinesTheta>().rmax()) {
            coord_eval2 -= std::floor(
                                   (coord_eval2 - ddc::discrete_space<BSplinesTheta>().rmin())
                                   / ddc::discrete_space<BSplinesTheta>().length())
                           * ddc::discrete_space<BSplinesTheta>().length();
        }
        ddc::Coordinate<DimR, DimTheta> coord_eval_new(coord_eval1, coord_eval2);
        return eval_no_bc(coord_eval_new, spline_coef, eval_type());
    }

    template <class EvalType>
    double eval_no_bc(
            ddc::Coordinate<DimR, DimTheta> coord_eval,
            host_t<ConstPolarSpline<PolarBSplinesType>> const spline_coef,
            EvalType const) const
    {
        static_assert(
                std::is_same_v<
                        EvalType,
                        eval_type> || std::is_same_v<EvalType, eval_deriv_r_type> || std::is_same_v<EvalType, eval_deriv_theta_type> || std::is_same_v<EvalType, eval_deriv_r_theta_type>);

        std::array<double, PolarBSplinesType::n_singular_basis()> singular_data;
        DSpan1D singular_vals(singular_data.data(), PolarBSplinesType::n_singular_basis());
        std::array<double, (BSplinesR::degree() + 1) * (BSplinesTheta::degree() + 1)> data;
        DSpan2D vals(data.data(), BSplinesR::degree() + 1, BSplinesTheta::degree() + 1);

        Idx<BSplinesR, BSplinesTheta> jmin;

        if constexpr (std::is_same_v<EvalType, eval_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_basis(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_r_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_r(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_theta_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_theta(singular_vals, vals, coord_eval);
        } else if constexpr (std::is_same_v<EvalType, eval_deriv_r_theta_type>) {
            jmin = ddc::discrete_space<PolarBSplinesType>()
                           .eval_deriv_r_and_theta(singular_vals, vals, coord_eval);
        }

        double y = 0.0;
        for (std::size_t i = 0; i < PolarBSplinesType::n_singular_basis(); ++i) {
            y += spline_coef.singular_spline_coef(Idx<PolarBSplinesType>(i)) * singular_vals(i);
        }
        Idx<BSplinesR> jmin_r = ddc::select<BSplinesR>(jmin);
        Idx<BSplinesTheta> jmin_theta = ddc::select<BSplinesTheta>(jmin);
        int nr = BSplinesR::degree() + 1;
        if (jmin_r < Idx<BSplinesR>(continuity + 1)) {
            nr = nr - (Idx<BSplinesR>(continuity + 1) - jmin_r);
            jmin_r = Idx<BSplinesR>(continuity + 1);
        }
        for (int i = 0; i < nr; ++i) {
            for (std::size_t j = 0; j < BSplinesTheta::degree() + 1; ++j) {
                y += spline_coef.spline_coef(jmin_r + i, jmin_theta + j) * vals(i, j);
            }
        }
        return y;
    }
};
```


