

# File discrete\_to\_cartesian.hpp

[**File List**](files.md) **>** [**coord\_transformations**](dir_67161c4ffadea73fddf46ea451c2f62c.md) **>** [**discrete\_to\_cartesian.hpp**](discrete__to__cartesian_8hpp.md)

[Go to the documentation of this file](discrete__to__cartesian_8hpp.md)


```C++
// SPDX-License-Identifier: MIT
#pragma once
#include <iostream>

#include <ddc/ddc.hpp>

#include "coord_transformation_tools.hpp"
#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
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
    using CoordJacobian = CoordArg;

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

    using IdxRTheta = typename IdxRangeRTheta::discrete_element_type;
    using IdxStepRTheta = typename IdxRangeRTheta::discrete_vector_type;

private:
    SplineType m_x_spline_representation;
    SplineType m_y_spline_representation;
    SplineEvaluator m_spline_evaluator;
    IdxRangeRTheta m_idx_range_singular_point;
    Coord<X, Y> m_o_point;

public:
    DiscreteToCartesian(
            SplineType curvilinear_to_x,
            SplineType curvilinear_to_y,
            SplineEvaluator const& evaluator,
            IdxRangeRTheta idx_range_singular_point)
        : m_x_spline_representation(curvilinear_to_x)
        , m_y_spline_representation(curvilinear_to_y)
        , m_spline_evaluator(evaluator)
        , m_idx_range_singular_point(idx_range_singular_point)
    {
        if constexpr (Kokkos::SpaceAccessibility<Kokkos::DefaultHostExecutionSpace, MemorySpace>::
                              accessible) {
            Coord<R, Theta> centre_coord(0.0, 0.0);
            m_o_point = (*this)(centre_coord);
        } else {
            IdxRangeRTheta idx_range_o_point
                    = m_idx_range_singular_point.take_first(IdxStepRTheta(1, 1));
            FieldMem<Coord<X, Y>, IdxRangeRTheta> coord_centre_field_alloc(idx_range_o_point);
            Field<Coord<X, Y>, IdxRangeRTheta> coord_centre_field
                    = get_field(coord_centre_field_alloc);
            using ExecSpace = typename SplineEvaluator::exec_space;
            (*this)(ExecSpace(), coord_centre_field);
            auto coord_centre_field_host = ddc::create_mirror_view_and_copy(coord_centre_field);
            m_o_point = coord_centre_field_host(idx_range_o_point.front());
        }
    }

    KOKKOS_INLINE_FUNCTION Coord<X, Y> o_point() const
    {
        return m_o_point;
    }

    KOKKOS_FUNCTION Coord<X, Y> operator()(
            Coord<curvilinear_tag_r, curvilinear_tag_theta> const& coord) const
    {
        const double x = m_spline_evaluator(coord, get_const_field(m_x_spline_representation));
        const double y = m_spline_evaluator(coord, get_const_field(m_y_spline_representation));
        return Coord<X, Y>(x, y);
    }

    template <class ExecSpace, class GridR, class GridTheta>
    void operator()(ExecSpace exec_space, Field<Coord<X, Y>, IdxRange<GridR, GridTheta>> coords)
    {
        static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible);
        static_assert(std::is_same_v<R, typename GridR::continuous_dimension_type>);
        static_assert(std::is_same_v<Theta, typename GridTheta::continuous_dimension_type>);
        ddc::parallel_for_each(
                ExecSpace(),
                get_idx_range(coords),
                KOKKOS_CLASS_LAMBDA(IdxRTheta idx) {
                    coords(idx) = (*this)(ddc::coordinate(idx));
                });
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

    template <class ExecSpace>
    void control_points(
            ExecSpace exec_space,
            Field<Coord<X, Y>, IdxRange<BSplineR, BSplineTheta>, MemorySpace> pts) const
    {
        static_assert(Kokkos::SpaceAccessibility<ExecSpace, MemorySpace>::accessible);
        ddc::parallel_for_each(
                exec_space,
                get_idx_range(pts),
                KOKKOS_CLASS_LAMBDA(Idx<BSplineR, BSplineTheta> idx) {
                    pts(idx) = control_point(idx);
                });
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
struct HasOPoint<DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>> : std::true_type
{
};

template <class X, class Y, class SplineEvaluator, class R, class Theta, class MemorySpace>
struct SingularOPointInvJacobian<DiscreteToCartesian<X, Y, SplineEvaluator, R, Theta, MemorySpace>>
    : std::true_type
{
};

} // namespace mapping_detail
```


