

# File discrete\_to\_cartesian.hpp

[**File List**](files.md) **>** [**mapping**](dir_5300298560c4bf255ab9f36681603d89.md) **>** [**discrete\_to\_cartesian.hpp**](discrete__to__cartesian_8hpp.md)

[Go to the documentation of this file](discrete__to__cartesian_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <iostream>

#include <ddc/ddc.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "mapping_tools.hpp"
#include "math_tools.hpp"
#include "tensor.hpp"
#include "view.hpp"

template <
        class X,
        class Y,
        class SplineEvaluator,
        class R = typename SplineEvaluator::continuous_dimension_type1,
        class Theta = typename SplineEvaluator::continuous_dimension_type2,
        class MemorySpace = typename SplineEvaluator::memory_space>
class DiscreteToCartesian
{
    static_assert(std::is_same_v<MemorySpace, typename SplineEvaluator::memory_space>);

public:
    using BSplineR = typename SplineEvaluator::bsplines_type1;
    using BSplineTheta = typename SplineEvaluator::bsplines_type2;

    using cartesian_tag_x = X;
    using cartesian_tag_y = Y;
    using curvilinear_tag_r = R;
    using curvilinear_tag_theta = Theta;

    using CoordArg = Coord<R, Theta>;
    using CoordResult = Coord<X, Y>;

    using X_cov = typename X::Dual;
    using Y_cov = typename Y::Dual;
    using R_cov = typename R::Dual;
    using Theta_cov = typename Theta::Dual;

private:
    using spline_idx_range = IdxRange<BSplineR, BSplineTheta>;

    using SplineType = DConstField<spline_idx_range, MemorySpace>;

    using IdxRangeRTheta = typename SplineEvaluator::evaluation_domain_type;
    using IdxRangeTheta = typename SplineEvaluator::evaluation_domain_type2;
    using IdxTheta = typename IdxRangeTheta::discrete_element_type;

private:
    SplineType m_x_spline_representation;
    SplineType m_y_spline_representation;
    SplineEvaluator m_spline_evaluator;
    IdxRangeRTheta m_idx_range_singular_point;

public:
    KOKKOS_FUNCTION DiscreteToCartesian(
            SplineType curvilinear_to_x,
            SplineType curvilinear_to_y,
            SplineEvaluator const& evaluator,
            IdxRangeRTheta idx_range_singular_point)
        : m_x_spline_representation(curvilinear_to_x)
        , m_y_spline_representation(curvilinear_to_y)
        , m_spline_evaluator(evaluator)
        , m_idx_range_singular_point(idx_range_singular_point)
    {
    }

    KOKKOS_FUNCTION Coord<X, Y> operator()(
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        const double x = m_spline_evaluator(coord, get_const_field(m_x_spline_representation));
        const double y = m_spline_evaluator(coord, get_const_field(m_y_spline_representation));
        return Coord<X, Y>(x, y);
    }

    KOKKOS_FUNCTION DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix(
            Coord<R, Theta> const& coord) const
    {
        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> jacobian_matrix;
        ddcHelper::get<X, R_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_x_spline_representation));
        ddcHelper::get<X, Theta_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_x_spline_representation));
        ddcHelper::get<Y, R_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_y_spline_representation));
        ddcHelper::get<Y, Theta_cov>(jacobian_matrix)
                = m_spline_evaluator.deriv_dim_2(coord, get_const_field(m_y_spline_representation));
        return jacobian_matrix;
    }

    template <class IndexTag1, class IndexTag2>
    KOKKOS_INLINE_FUNCTION double jacobian_component(Coord<R, Theta> coord) const
    {
        static_assert(ddc::in_tags_v<IndexTag1, VectorIndexSet<X, Y>>);
        static_assert(ddc::in_tags_v<IndexTag2, VectorIndexSet<R_cov, Theta_cov>>);

        if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (1,1), i.e dx/dr
            return m_spline_evaluator
                    .deriv_dim_1(coord, get_const_field(m_x_spline_representation));
        } else if constexpr (std::is_same_v<IndexTag1, X> && std::is_same_v<IndexTag2, Theta_cov>) {
            // Component (1,2), i.e dx/dtheta
            return m_spline_evaluator
                    .deriv_dim_2(coord, get_const_field(m_x_spline_representation));
        } else if constexpr (std::is_same_v<IndexTag1, Y> && std::is_same_v<IndexTag2, R_cov>) {
            // Component (2,1), i.e dy/dr
            return m_spline_evaluator
                    .deriv_dim_1(coord, get_const_field(m_y_spline_representation));
        } else {
            // Component (2,2), i.e dy/dtheta
            return m_spline_evaluator
                    .deriv_dim_2(coord, get_const_field(m_y_spline_representation));
        }
    }

    KOKKOS_FUNCTION double jacobian(
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        Tensor J = jacobian_matrix(coord);
        return ddcHelper::get<X, R_cov>(J) * ddcHelper::get<Y, Theta_cov>(J)
               - ddcHelper::get<Y, R_cov>(J) * ddcHelper::get<X, Theta_cov>(J);
    }

    KOKKOS_INLINE_FUNCTION DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>>
    first_order_jacobian_matrix_r_rtheta(
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        DTensor<VectorIndexSet<X, Y>, VectorIndexSet<R_cov, Theta_cov>> J;
        ddcHelper::get<X, R_cov>(J)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_x_spline_representation));
        ddcHelper::get<X, Theta_cov>(J)
                = m_spline_evaluator
                          .deriv_1_and_2(coord, get_const_field(m_x_spline_representation));
        ddcHelper::get<Y, R_cov>(J)
                = m_spline_evaluator.deriv_dim_1(coord, get_const_field(m_y_spline_representation));
        ddcHelper::get<Y, Theta_cov>(J)
                = m_spline_evaluator
                          .deriv_1_and_2(coord, get_const_field(m_y_spline_representation));
        return J;
    }

    KOKKOS_INLINE_FUNCTION IdxRangeRTheta idx_range_singular_point() const
    {
        return m_idx_range_singular_point;
    }

    KOKKOS_INLINE_FUNCTION const Coord<X, Y> control_point(
            Idx<BSplineR, BSplineTheta> const& el) const
    {
        return Coord<X, Y>(m_x_spline_representation(el), m_y_spline_representation(el));
    }
};


namespace mapping_detail {
template <
        class X,
        class Y,
        class SplineEvaluator,
        class R,
        class Theta,
        class MemorySpace,
        class ExecSpace>
struct MappingAccessibility<
        ExecSpace,
        DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>
{
    static constexpr bool value = Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible;
};

template <class X, class Y, class SplineEvaluator, class R, class Theta, class MemorySpace>
struct IsCurvilinear2DMapping<DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>
    : std::true_type
{
};

template <class X, class Y, class SplineEvaluator, class R, class Theta, class MemorySpace>
struct SingularOPointInvJacobian<DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>
    : std::true_type
{
};

} // namespace mapping_detail
```


