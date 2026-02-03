

# File polar\_spline\_evaluator.hpp

[**File List**](files.md) **>** [**interpolation**](dir_264890e5c091f8c8d7fe1f842870c25e.md) **>** [**polar\_splines**](dir_a6779ae02b71d57f488d261458bab1ce.md) **>** [**polar\_spline\_evaluator.hpp**](polar__spline__evaluator_8hpp.md)

[Go to the documentation of this file](polar__spline__evaluator_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once

#include "view.hpp"

template <class ExecSpace, class MemorySpace, class PolarBSplinesType, class OuterExtrapolationRule>
class PolarSplineEvaluator
{
public:
    using bsplines_type = PolarBSplinesType;
    using BSplinesR = typename PolarBSplinesType::BSplinesR_tag;
    using BSplinesTheta = typename PolarBSplinesType::BSplinesTheta_tag;
    using R = typename BSplinesR::continuous_dimension_type;
    using Theta = typename BSplinesTheta::continuous_dimension_type;

    using exec_space = ExecSpace;

    using memory_space = MemorySpace;

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

    KOKKOS_FUNCTION double operator()(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval(coord_eval, spline_coef);
    }

    template <class Domain>
    void operator()(
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval(coords_eval(i), spline_coef);
                });
    }

    template <class Domain>
    void operator()(
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval(ddc::coordinate(i), spline_coef);
                });
    }

    template <class DerivDim>
    KOKKOS_FUNCTION double deriv(
            Idx<ddc::Deriv<DerivDim>> deriv_order,
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, deriv_order);
    }

    [[deprecated("Use deriv(Idx<ddc::Deriv<R>>(1), ...) instead")]] KOKKOS_FUNCTION double
    deriv_dim_1(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, Idx<ddc::Deriv<R>>(1));
    }

    [[deprecated("Use deriv(Idx<ddc::Deriv<Theta>>(1), ...) instead")]] KOKKOS_FUNCTION double
    deriv_dim_2(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, Idx<ddc::Deriv<Theta>>(1));
    }

    [[deprecated("Use eval_deriv(..., Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1)) "
                 "instead")]] KOKKOS_FUNCTION double
    deriv_1_and_2(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        return eval_no_bc(coord_eval, spline_coef, Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1));
    }

    template <class Domain, class... DerivDims>
    void deriv(
            Idx<DerivDims...> const deriv_order,
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, deriv_order);
                });
    }

    template <class Domain, class... DerivDims>
    void deriv(
            Idx<DerivDims...> const deriv_order,
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(ddc::coordinate(i), spline_coef, deriv_order);
                });
    }

    template <class Domain>
    [[deprecated("Use deriv(Idx<ddc::Deriv<R>>(1), ...) instead")]] void deriv_dim_1(
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(coords_eval(i), spline_coef, Idx<ddc::Deriv<R>>(1));
                });
    }

    template <class Domain>
    [[deprecated("Use deriv(Idx<ddc::Deriv<R>>(1), ...) instead")]] void deriv_dim_1(
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i)
                            = eval_no_bc(ddc::coordinate(i), spline_coef, Idx<ddc::Deriv<R>>(1));
                });
    }

    template <class Domain>
    [[deprecated("Use deriv(Idx<ddc::Deriv<Theta>>(1), ...) instead")]] void deriv_dim_2(
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i)
                            = eval_no_bc(coords_eval(i), spline_coef, Idx<ddc::Deriv<Theta>>(1));
                });
    }

    template <class Domain>
    [[deprecated("Use deriv(Idx<ddc::Deriv<Theta>>(1), ...) instead")]] void deriv_dim_2(
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(
                            ddc::coordinate(i),
                            spline_coef,
                            Idx<ddc::Deriv<Theta>>(1));
                });
    }

    template <class Domain>
    [[deprecated("Use eval_deriv(..., Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1)) instead")]] void
    deriv_dim_1_and_2(
            DField<Domain, MemorySpace> const spline_eval,
            ConstField<Coord<R, Theta>, Domain, MemorySpace> const coords_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(coords_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(
                            coords_eval(i),
                            spline_coef,
                            Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1));
                });
    }

    template <class Domain>
    [[deprecated("Use eval_deriv(..., Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1)) instead")]] void
    deriv_dim_1_and_2(
            DField<Domain, MemorySpace> const spline_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        using IdxEval = typename Domain::discrete_element_type;
        ddc::parallel_for_each(
                exec_space(),
                get_idx_range(spline_eval),
                KOKKOS_CLASS_LAMBDA(IdxEval i) {
                    spline_eval(i) = eval_no_bc(
                            ddc::coordinate(i),
                            spline_coef,
                            Idx<ddc::Deriv<R>, ddc::Deriv<Theta>>(1, 1));
                });
    }

private:
    KOKKOS_FUNCTION double eval(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef) const
    {
        const double coord_eval1 = ddc::get<R>(coord_eval);
        double coord_eval2 = ddc::get<Theta>(coord_eval);
        if (coord_eval1 > ddc::discrete_space<BSplinesR>().rmax()) {
            return m_outer_bc(coord_eval, spline_coef);
        }
        if (coord_eval2 < ddc::discrete_space<BSplinesTheta>().rmin()
            || coord_eval2 > ddc::discrete_space<BSplinesTheta>().rmax()) {
            coord_eval2 -= Kokkos::floor(
                                   (coord_eval2 - ddc::discrete_space<BSplinesTheta>().rmin())
                                   / ddc::discrete_space<BSplinesTheta>().length())
                           * ddc::discrete_space<BSplinesTheta>().length();
        }
        Coord<R, Theta> coord_eval_new(coord_eval1, coord_eval2);
        return eval_no_bc(coord_eval_new, spline_coef, Idx<>());
    }

    template <class... DerivDims>
    KOKKOS_FUNCTION double eval_no_bc(
            Coord<R, Theta> coord_eval,
            DConstField<IdxRange<PolarBSplinesType>, MemorySpace> const spline_coef,
            Idx<DerivDims...> const deriv_order) const
    {
        std::array<double, PolarBSplinesType::n_singular_basis()> singular_data;
        DSpan1D singular_vals(singular_data.data(), PolarBSplinesType::n_singular_basis());
        std::array<double, (BSplinesR::degree() + 1) * (BSplinesTheta::degree() + 1)> data;
        DSpan2D vals(data.data(), BSplinesR::degree() + 1, BSplinesTheta::degree() + 1);

        Idx<BSplinesR, BSplinesTheta> jmin
                = ddc::discrete_space<PolarBSplinesType>()
                          .eval_deriv(singular_vals, vals, coord_eval, deriv_order);

        double y = 0.0;
        for (std::size_t i = 0; i < PolarBSplinesType::n_singular_basis(); ++i) {
            y += spline_coef(Idx<PolarBSplinesType>(i)) * singular_vals(i);
        }
        Idx<BSplinesR> jmin_r(jmin);
        Idx<BSplinesTheta> jmin_theta(jmin);

        DConstField<IdxRange<BSplinesR, BSplinesTheta>, MemorySpace> spline_coef_2d
                = PolarBSplinesType::get_tensor_product_subset(spline_coef);
        IdxRange<BSplinesR, BSplinesTheta> tensor_prod_idx_range = get_idx_range(spline_coef_2d);
        IdxRange<BSplinesTheta> tensor_prod_idx_range_theta(tensor_prod_idx_range);
        Idx<BSplinesTheta> idx_theta_max = tensor_prod_idx_range_theta.back();
        IdxStep<BSplinesTheta> n_idx_theta = tensor_prod_idx_range_theta.extents();
        for (int i = 0; i < BSplinesR::degree() + 1; ++i) {
            for (std::size_t j = 0; j < BSplinesTheta::degree() + 1; ++j) {
                Idx<BSplinesTheta> idx_theta = jmin_theta + j;
                if (idx_theta > idx_theta_max) {
                    idx_theta -= n_idx_theta;
                }
                y += spline_coef_2d(jmin_r + i, idx_theta) * vals(i, j);
            }
        }
        return y;
    }
};
```


