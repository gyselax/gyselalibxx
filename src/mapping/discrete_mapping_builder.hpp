// SPDX-License-Identifier: MIT
#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/splines.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "discrete_to_cartesian.hpp"

/**
 * @brief A class to create a DiscreteToCartesian instance from an analytical mapping.
 * This class creates and stores splines memory spaces describing the analytical mapping.
 * The discrete mapping is then created using the splines without copying data.
 *
 * @tparam X The first Cartesian dimension.
 * @tparam Y The second Cartesian dimension.
 * @tparam SplineBuilder An operator for building spline coefficients.
 * @tparam SplineEvaluator An operator for evaluating a spline.
 */
template <class X, class Y, class SplineBuilder, class SplineEvaluator>
class DiscreteToCartesianBuilder
{
    static_assert(std::is_same_v<
                  typename SplineBuilder::memory_space,
                  typename SplineEvaluator::memory_space>);
    static_assert(std::is_same_v<
                  typename SplineBuilder::exec_space,
                  typename SplineEvaluator::exec_space>);
    static_assert(std::is_same_v<typename SplineBuilder::batch_domain_type, IdxRange<>>);

public:
    /// The type of the mapping that will be created.
    using MappingType = DiscreteToCartesian<X, Y, SplineEvaluator>;

private:
    using ExecSpace = typename SplineBuilder::exec_space;
    using MemorySpace = typename SplineBuilder::memory_space;

    using BSplinesR = typename SplineBuilder::bsplines_type1;
    using BSplinesTheta = typename SplineBuilder::bsplines_type2;

    using GridR = typename SplineBuilder::interpolation_domain_type1;
    using GridTheta = typename SplineBuilder::interpolation_domain_type2;

    using IdxRangeSplines = IdxRange<BSplinesR, BSplinesTheta>;
    using IdxRangeInterpolationPoints = typename SplineBuilder::interpolation_domain_type;
    using IdxInterpolationPoints = typename IdxRangeInterpolationPoints::discrete_element_type;

    using SplineCoeffsMem = DFieldMem<IdxRangeSplines, MemorySpace>;

    using SplineCoeffs = typename SplineCoeffsMem::span_type;

    using InterpolationFieldMem = DFieldMem<IdxRangeInterpolationPoints, MemorySpace>;

    using InterpolationField = typename InterpolationFieldMem::span_type;

private:
    SplineCoeffsMem m_curvilinear_to_x_spline_alloc;
    SplineCoeffsMem m_curvilinear_to_y_spline_alloc;
    SplineEvaluator m_evaluator;
    IdxRange<GridR, GridTheta> m_idx_range_singular_point;

public:
    /**
     * @brief Create an instance of the class capable of providing a DiscreteToCartesian class instance.
     *
     * @param[in] exec_space The execution space where this class runs any for loops.
     * @param[in] analytical_mapping The analytical mapping to be described by this discrete mapping.
     * @param[in] builder A spline builder to be used to create a spline approximating the analytical mapping.
     * @param[in] evaluator A spline evaluator to be used to evaluate a spline approximating the analytical mapping.
     */
    template <class Mapping>
    DiscreteToCartesianBuilder(
            ExecSpace exec_space,
            Mapping const& analytical_mapping,
            SplineBuilder const& builder,
            SplineEvaluator const& evaluator)
        : m_curvilinear_to_x_spline_alloc(get_spline_idx_range(builder))
        , m_curvilinear_to_y_spline_alloc(get_spline_idx_range(builder))
        , m_evaluator(evaluator)
    {
        SplineCoeffs curvilinear_to_x_spline = get_field(m_curvilinear_to_x_spline_alloc);
        SplineCoeffs curvilinear_to_y_spline = get_field(m_curvilinear_to_y_spline_alloc);
        InterpolationFieldMem curvilinear_to_x_vals_alloc(builder.interpolation_domain());
        InterpolationFieldMem curvilinear_to_y_vals_alloc(builder.interpolation_domain());
        InterpolationField curvilinear_to_x_vals = get_field(curvilinear_to_x_vals_alloc);
        InterpolationField curvilinear_to_y_vals = get_field(curvilinear_to_y_vals_alloc);

        set_curvilinear_to_cartesian_values(
                curvilinear_to_x_vals,
                curvilinear_to_y_vals,
                analytical_mapping,
                builder.interpolation_domain());

        builder(get_field(curvilinear_to_x_spline), get_const_field(curvilinear_to_x_vals));
        builder(get_field(curvilinear_to_y_spline), get_const_field(curvilinear_to_y_vals));

        IdxRange<GridR, GridTheta> interp_domain(builder.interpolation_domain());
        IdxStep<GridR, GridTheta>
                n_points_singular_domain(1, interp_domain.template extent<GridTheta>().value());
        m_idx_range_singular_point = interp_domain.take_first(n_points_singular_domain);
    }

    /**
     * @brief Get a DiscreteToCartesian class instance.
     *
     * @return An instance of the mapping.
     */
    DiscreteToCartesian<X, Y, SplineEvaluator> operator()() const
    {
        return DiscreteToCartesian<X, Y, SplineEvaluator>(
                get_const_field(m_curvilinear_to_x_spline_alloc),
                get_const_field(m_curvilinear_to_y_spline_alloc),
                m_evaluator,
                m_idx_range_singular_point);
    }

    /**
     * @brief Fill in the curvilinear fields with interpolation 
     * points mapped with the given analytical mapping. 
     *
     * This function should be private. It is not due to the inclusion of a KOKKOS_LAMBDA
     *
     * @tparam Mapping Type of the analytical mapping. 
     * @param[out] curvilinear_to_x_vals Field of coordinate on X. 
     * @param[out] curvilinear_to_y_vals Field of coordinate on Y. 
     * @param[in] analytical_mapping Analytical mapping. 
     * @param[in] interpolation_idx_range Index range of an interpolation grid. 
     */
    template <class Mapping>
    void set_curvilinear_to_cartesian_values(
            InterpolationField const& curvilinear_to_x_vals,
            InterpolationField const& curvilinear_to_y_vals,
            Mapping const& analytical_mapping,
            IdxRangeInterpolationPoints const& interpolation_idx_range)
    {
        using CurvilinearCoeff
                = Coord<typename Mapping::curvilinear_tag_r,
                        typename Mapping::curvilinear_tag_theta>;
        using CartesianCoeff = Coord<X, Y>;

        ddc::parallel_for_each(
                ExecSpace(),
                interpolation_idx_range,
                KOKKOS_LAMBDA(IdxInterpolationPoints el) {
                    CurvilinearCoeff polar_coord(ddc::coordinate(el));
                    CartesianCoeff cart_coord = analytical_mapping(polar_coord);
                    curvilinear_to_x_vals(el) = ddc::select<X>(cart_coord);
                    curvilinear_to_y_vals(el) = ddc::select<Y>(cart_coord);
                });
    }
};

/**
 * @brief A class to create a DiscreteToCartesian instance from an analytical mapping.
 * This class creates an instance which uses more refined splines than the provided builder and
 * evaluator.
 * This class creates and stores splines memory spaces describing the analytical mapping.
 * The discrete mapping is then created using the splines without copying data.
 *
 * @tparam X The first Cartesian dimension.
 * @tparam Y The second Cartesian dimension.
 * @tparam SplineBuilder An operator for building spline coefficients.
 * @tparam SplineEvaluator An operator for evaluating a spline.
 * @tparam ncells_r The number of cells in the refined spline in the radial direction.
 * @tparam ncells_theta The number of cells in the refined spline in the radial direction.
 */
template <
        class X,
        class Y,
        class SplineBuilder,
        class SplineEvaluator,
        int ncells_r,
        int ncells_theta>
class RefinedDiscreteToCartesianBuilder
{
    static_assert(std::is_same_v<
                  typename SplineBuilder::memory_space,
                  typename SplineEvaluator::memory_space>);
    static_assert(std::is_same_v<
                  typename SplineBuilder::exec_space,
                  typename SplineEvaluator::exec_space>);

private:
    using ExecSpace = typename SplineBuilder::exec_space;
    using MemorySpace = typename SplineBuilder::memory_space;

    using BSplinesROriginal = typename SplineBuilder::bsplines_type1;
    using BSplinesThetaOriginal = typename SplineBuilder::bsplines_type2;

    using R = typename BSplinesROriginal::continuous_dimension_type;
    using Theta = typename BSplinesThetaOriginal::continuous_dimension_type;

public:
    /// @brief The type of the radial B-splines on which the new mapping will be defined.
    struct BSplinesRRefined
        : std::conditional_t<
                  BSplinesROriginal::is_uniform(),
                  ddc::UniformBSplines<R, BSplinesROriginal::degree()>,
                  ddc::NonUniformBSplines<R, BSplinesROriginal::degree()>>
    {
    };

    /// @brief The type of the poloidal B-splines on which the new mapping will be defined.
    struct BSplinesThetaRefined
        : std::conditional_t<
                  BSplinesThetaOriginal::is_uniform(),
                  ddc::UniformBSplines<Theta, BSplinesThetaOriginal::degree()>,
                  ddc::NonUniformBSplines<Theta, BSplinesThetaOriginal::degree()>>
    {
    };

private:
    using GrevillePointsR = ddc::GrevilleInterpolationPoints<
            BSplinesRRefined,
            SplineBuilder::builder_type1::s_bc_xmin,
            SplineBuilder::builder_type1::s_bc_xmax>;

    using GrevillePointsTheta = ddc::GrevilleInterpolationPoints<
            BSplinesThetaRefined,
            SplineBuilder::builder_type2::s_bc_xmin,
            SplineBuilder::builder_type2::s_bc_xmax>;

public:
    /// @brief The type of the grid of radial points on which the new mapping will be defined.
    struct GridRRefined : GrevillePointsR::interpolation_discrete_dimension_type
    {
    };

    /// @brief The type of the grid of poloidal points on which the new mapping will be defined.
    struct GridThetaRefined : GrevillePointsTheta::interpolation_discrete_dimension_type
    {
    };

private:
    using GridROriginal = typename SplineBuilder::interpolation_domain_type1;
    using GridThetaOriginal = typename SplineBuilder::interpolation_domain_type2;

    template <class Builder>
    struct Build_BuilderType;

    template <
            ddc::BoundCond BcLower1,
            ddc::BoundCond BcUpper1,
            ddc::BoundCond BcLower2,
            ddc::BoundCond BcUpper2,
            ddc::SplineSolver Solver>
    struct Build_BuilderType<ddc::SplineBuilder2D<
            ExecSpace,
            MemorySpace,
            BSplinesROriginal,
            BSplinesThetaOriginal,
            GridROriginal,
            GridThetaOriginal,
            BcLower1,
            BcUpper1,
            BcLower2,
            BcUpper2,
            Solver,
            GridROriginal,
            GridThetaOriginal>>
    {
        using type = ddc::SplineBuilder2D<
                ExecSpace,
                MemorySpace,
                BSplinesRRefined,
                BSplinesThetaRefined,
                GridRRefined,
                GridThetaRefined,
                BcLower1,
                BcUpper1,
                BcLower2,
                BcUpper2,
                Solver,
                GridRRefined,
                GridThetaRefined>;
    };

    using RefinedSplineBuilder = typename Build_BuilderType<SplineBuilder>::type;

    using RefinedSplineEvaluator = ddc::SplineEvaluator2D<
            ExecSpace,
            MemorySpace,
            BSplinesRRefined,
            BSplinesThetaRefined,
            GridRRefined,
            GridThetaRefined,
            typename SplineEvaluator::lower_extrapolation_rule_1_type,
            typename SplineEvaluator::upper_extrapolation_rule_1_type,
            typename SplineEvaluator::lower_extrapolation_rule_2_type,
            typename SplineEvaluator::upper_extrapolation_rule_2_type,
            GridRRefined,
            GridThetaRefined>;

    using IdxRangeSplines = IdxRange<BSplinesRRefined, BSplinesThetaRefined>;
    using IdxRangeInterpolationPoints = IdxRange<GridRRefined, GridThetaRefined>;
    using IdxInterpolationPoints = typename IdxRangeInterpolationPoints::discrete_element_type;

    using SplineCoeffsMem = DFieldMem<IdxRangeSplines, MemorySpace>;

    using SplineCoeffs = typename SplineCoeffsMem::span_type;

    using InterpolationFieldMem = DFieldMem<IdxRangeInterpolationPoints, MemorySpace>;

    using InterpolationField = typename InterpolationFieldMem::span_type;

public:
    /// The type of the mapping that will be created.
    using MappingType = DiscreteToCartesian<X, Y, RefinedSplineEvaluator>;

private:
    SplineCoeffsMem m_curvilinear_to_x_spline_alloc;
    SplineCoeffsMem m_curvilinear_to_y_spline_alloc;
    RefinedSplineEvaluator m_evaluator;
    IdxRange<GridRRefined, GridThetaRefined> m_idx_range_singular_point;

public:
    /**
     * @brief Create an instance of the class capable of providing a DiscreteToCartesian class instance.
     *
     * @param[in] exec_space The execution space where this class runs any for loops.
     * @param[in] analytical_mapping The analytical mapping to be described by this discrete mapping.
     * @param[in] builder A spline builder to be used to create a spline approximating the analytical mapping.
     * @param[in] evaluator A spline evaluator to be used to evaluate a spline approximating the analytical mapping.
     */
    template <class Mapping>
    RefinedDiscreteToCartesianBuilder(
            ExecSpace exec_space,
            Mapping const& analytical_mapping,
            SplineBuilder const& builder,
            SplineEvaluator const& evaluator)
        : m_evaluator(
                evaluator.lower_extrapolation_rule_dim_1(),
                evaluator.upper_extrapolation_rule_dim_1(),
                evaluator.lower_extrapolation_rule_dim_2(),
                evaluator.upper_extrapolation_rule_dim_2())
    {
        using CoordR = Coord<R>;
        using CoordTheta = Coord<Theta>;

        const CoordR r_min = ddc::discrete_space<BSplinesROriginal>().rmin();
        const CoordR r_max = ddc::discrete_space<BSplinesROriginal>().rmax();

        const CoordTheta theta_min = ddc::discrete_space<BSplinesThetaOriginal>().rmin();
        const CoordTheta theta_max = ddc::discrete_space<BSplinesThetaOriginal>().rmax();

        IdxStep<BSplinesRRefined> const r_ncells_idx_step(ncells_r);
        IdxStep<BSplinesThetaRefined> const theta_ncells_idx_step(ncells_theta);

        if constexpr (BSplinesRRefined::is_uniform()) {
            ddc::init_discrete_space<BSplinesRRefined>(r_min, r_max, r_ncells_idx_step);
            ddc::init_discrete_space<
                    BSplinesThetaRefined>(theta_min, theta_max, theta_ncells_idx_step);
        } else {
            double const dr((r_max - r_min) / ncells_r);
            double const dp((theta_max - theta_min) / ncells_theta);

            std::vector<CoordR> r_break_points(ncells_r + 1);
            std::vector<CoordTheta> theta_break_points(ncells_theta + 1);

            for (int i(0); i < ncells_r + 1; ++i) {
                r_break_points[i] = r_min + i * dr;
            }
            for (int i(0); i < ncells_theta + 1; ++i) {
                theta_break_points[i] = theta_min + i * dp;
            }

            ddc::init_discrete_space<BSplinesRRefined>(r_break_points);
            ddc::init_discrete_space<BSplinesThetaRefined>(theta_break_points);
        }

        ddc::init_discrete_space<GridRRefined>(
                GrevillePointsR::template get_sampling<GridRRefined>());
        ddc::init_discrete_space<GridThetaRefined>(
                GrevillePointsTheta::template get_sampling<GridThetaRefined>());

        IdxRange<GridRRefined> const interpolation_domain_r(
                GrevillePointsR::template get_domain<GridRRefined>());
        IdxRange<GridThetaRefined> const interpolation_domain_theta(
                GrevillePointsTheta::template get_domain<GridThetaRefined>());
        IdxRange<GridRRefined, GridThetaRefined> const
                refined_domain(interpolation_domain_r, interpolation_domain_theta);

        IdxStep<GridRRefined, GridThetaRefined> n_points_singular_domain(
                1,
                refined_domain.template extent<GridThetaRefined>().value());
        m_idx_range_singular_point = refined_domain.take_first(n_points_singular_domain);

        // Operators on the refined grid
        RefinedSplineBuilder refined_builder(refined_domain);

        IdxRange<BSplinesRRefined, BSplinesThetaRefined> const spline_domain
                = get_spline_idx_range(refined_builder);

        // Compute the B-splines coefficients of the analytical mapping
        m_curvilinear_to_x_spline_alloc = SplineCoeffsMem(spline_domain);
        m_curvilinear_to_y_spline_alloc = SplineCoeffsMem(spline_domain);
        SplineCoeffs curvilinear_to_x_spline = get_field(m_curvilinear_to_x_spline_alloc);
        SplineCoeffs curvilinear_to_y_spline = get_field(m_curvilinear_to_y_spline_alloc);
        InterpolationFieldMem curvilinear_to_x_vals_alloc(refined_domain);
        InterpolationFieldMem curvilinear_to_y_vals_alloc(refined_domain);
        InterpolationField curvilinear_to_x_vals = get_field(curvilinear_to_x_vals_alloc);
        InterpolationField curvilinear_to_y_vals = get_field(curvilinear_to_y_vals_alloc);

        set_curvilinear_to_cartesian_values(
                curvilinear_to_x_vals,
                curvilinear_to_y_vals,
                analytical_mapping,
                refined_builder.interpolation_domain());

        refined_builder(curvilinear_to_x_spline, get_const_field(curvilinear_to_x_vals));
        refined_builder(curvilinear_to_y_spline, get_const_field(curvilinear_to_y_vals));
    }

    /**
     * @brief Get a DiscreteToCartesian class instance.
     *
     * @return An instance of the mapping.
     */
    DiscreteToCartesian<X, Y, RefinedSplineEvaluator> operator()() const
    {
        return DiscreteToCartesian<X, Y, RefinedSplineEvaluator>(
                get_const_field(m_curvilinear_to_x_spline_alloc),
                get_const_field(m_curvilinear_to_y_spline_alloc),
                m_evaluator,
                m_idx_range_singular_point);
    }

    /**
     * @brief Fill in the curvilinear fields with interpolation 
     * points mapped with the given analytical mapping. 
     *
     * This function should be private. It is not due to the inclusion of a KOKKOS_LAMBDA
     *
     * @tparam Mapping Type of the analytical mapping. 
     * @param[out] curvilinear_to_x_vals Field of coordinate on X. 
     * @param[out] curvilinear_to_y_vals Field of coordinate on Y. 
     * @param[in] analytical_mapping Analytical mapping. 
     * @param[in] interpolation_idx_range Index range of an interpolation grid. 
     */
    template <class Mapping>
    void set_curvilinear_to_cartesian_values(
            InterpolationField curvilinear_to_x_vals,
            InterpolationField curvilinear_to_y_vals,
            Mapping const& analytical_mapping,
            IdxRangeInterpolationPoints const& interpolation_idx_range)
    {
        using CurvilinearCoeff
                = Coord<typename Mapping::curvilinear_tag_r,
                        typename Mapping::curvilinear_tag_theta>;
        using CartesianCoeff = Coord<X, Y>;

        ddc::parallel_for_each(
                ExecSpace(),
                interpolation_idx_range,
                KOKKOS_LAMBDA(IdxInterpolationPoints el) {
                    CurvilinearCoeff polar_coord(ddc::coordinate(el));
                    CartesianCoeff cart_coord = analytical_mapping(polar_coord);
                    curvilinear_to_x_vals(el) = ddc::select<X>(cart_coord);
                    curvilinear_to_y_vals(el) = ddc::select<Y>(cart_coord);
                });
    }
};
