// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolator_2d.hpp"
#include "indexed_tensor.hpp"
#include "l_norm_tools.hpp"
#include "metric_tensor_evaluator.hpp"
#include "spline_interpolator_2d.hpp"
#include "spline_polar_foot_finder.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"



/**
 * @brief Define an advection operator on 2D @f$(r, \theta)@f$ domain.
 *
 * The advection operator uses a backward semi-Lagrangian method. The method is based on
 * the property that the solution is constant along the characteristics.
 *
 * For the following equation:
 * @f$\partial_t f(t,x) + A(t, x) \cdot \nabla_x f(t,x) = 0,  @f$
 *
 * we write the characteristics:
 * @f$ \partial_t X(t; s, x) = A(t, X(t; s, x)), \qquad \text{ with } X(s; s, x) = x. @f$
 *
 * Then the property gives us:
 * @f$ f(t, x) = f(0, X(t; 0, x)), \quad \forall t. @f$
 *
 *
 * So the first step of the advection operator is to compute the feet of the characteristics
 * @f$ X(t; t+\Delta t, x_i) @f$ for each mesh point @f$ x_i @f$.
 *
 * For the second step, we interpolate the function at the computed feet of the characteristics, 
 * and obtain the function at the next time step: 
 * @f$ f(t + \Delta t, x) = f(t, X(t; t+\Delta, x))@f$.
 *
 *
 * Different time integration methods are implemented to solve the equation of the characteristics.
 * They are defined in the IPolarFootFinder class.
 *
 * The feet can be advected on different domains (physical domain or pseudo-physical domain)
 * which are determined in the SplinePolarFootFinder operator. 
 *
 * The interpolation of the function is always done in the logical domain,
 * where the B-splines are defined. 
 *
 *
 * @see IPolarFootFinder
 */
template <class FootFinder, class LogicalToPhysicalMapping, class InterpolatorPolar>
class BslAdvectionPolar
{
    using R = typename LogicalToPhysicalMapping::curvilinear_tag_r;
    using Theta = typename LogicalToPhysicalMapping::curvilinear_tag_theta;

    using CoordRTheta = typename LogicalToPhysicalMapping::CoordArg;
    using CoordXY = typename LogicalToPhysicalMapping::CoordResult;

    using CartesianBasis = ddc::to_type_seq_t<CoordXY>;
    using CurvilinearBasis = ddc::to_type_seq_t<CoordRTheta>;

    using IdxRangeBatched = typename FootFinder::IdxRangeOperator;

    using GridR = find_grid_t<R, ddc::to_type_seq_t<IdxRangeBatched>>;
    using GridTheta = find_grid_t<Theta, ddc::to_type_seq_t<IdxRangeBatched>>;

    using IdxRangeRTheta = IdxRange<GridR, GridTheta>;
    using IdxRangeR = IdxRange<GridR>;
    using IdxRangeTheta = IdxRange<GridTheta>;

    using IdxRTheta = typename IdxRangeRTheta::discrete_element_type;
    using IdxBatched = typename IdxRangeBatched::discrete_element_type;
    using IdxStepR = typename IdxRangeR::discrete_vector_type;

    using MemorySpace = typename FootFinder::memory_space;

    using DFieldFDistribu = DField<IdxRangeBatched, MemorySpace>;

    using CFieldMemFeetRTheta = FieldMem<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using CFieldFeetRTheta = Field<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using DVectorFieldMemAdvectionXY
            = DVectorFieldMem<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXY = DVectorField<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorConstFieldAdvectionXY
            = DVectorConstField<IdxRangeBatched, CartesianBasis, MemorySpace>;

    using DVectorFieldAdvectionRTheta
            = DVectorField<IdxRangeBatched, CurvilinearBasis, MemorySpace>;
    using DVectorConstFieldAdvectionRTheta
            = DVectorConstField<IdxRangeBatched, CurvilinearBasis, MemorySpace>;

private:
    InterpolatorPolar const& m_interpolator;

    FootFinder const& m_find_feet;

    LogicalToPhysicalMapping const& m_logical_to_physical_mapping;


public:
    /**
     * @brief Instantiate an advection operator.
     *
     * @param [in] function_interpolator
     *       The polar interpolator to interpolate the function once the
     *      characteristics have been computed.
     * @param[in] foot_finder
     *      An IFootFinder which computes the feet of the characteristics.
     * @param[in] logical_to_physical_mapping
     *      The mapping function from the logical domain to the physical
     *      domain. 
     *
     * @tparam IFootFinder
     *      A child class of IFootFinder.
     */
    BslAdvectionPolar(
            InterpolatorPolar const& function_interpolator,
            FootFinder const& foot_finder,
            LogicalToPhysicalMapping const& logical_to_physical_mapping)
        : m_interpolator(function_interpolator)
        , m_find_feet(foot_finder)
        , m_logical_to_physical_mapping(logical_to_physical_mapping)
    {
    }

    ~BslAdvectionPolar() = default;


    /**
     * @brief Allocate a Field of the advected function.
     *
     * @param [in, out] allfdistribu
     *      A Field containing the values of the function we want to advect.
     * @param [in] advection_field_xy
     *      A field of vectors defined on the Cartesian basis containing the values
     *      of the advection field at each point on the logical grid.
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvectionXY advection_field_xy,
            double dt) const
    {
        // Pre-allocate some memory to prevent allocation later in loop
        std::unique_ptr<IInterpolator2D<IdxRangeRTheta, IdxRangeBatched>> const interpolator_ptr
                = m_interpolator.preallocate();

        // Initialise the feet
        CFieldMemFeetRTheta feet_rtheta_alloc(get_idx_range(advection_field_xy));
        CFieldFeetRTheta feet_rtheta = get_field(feet_rtheta_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(advection_field_xy),
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    IdxRTheta const irtheta(idx);
                    feet_rtheta(idx) = ddc::coordinate(irtheta);
                });

        // Compute the feet of the characteristics at tn -----------------------------------------
        m_find_feet(feet_rtheta, get_const_field(advection_field_xy), dt);

        // Interpolate the function on the feet of the characteristics. --------------------------
        (*interpolator_ptr)(get_field(allfdistribu), get_const_field(feet_rtheta));

        return allfdistribu;
    }


    /**
     * @brief Allocate a Field to the advected function.
     *
     * @param [in, out] allfdistribu
     *      A Field containing the values of the function we want to advect.
     * @param [in] advection_field_rtheta
     *      A field of vectors defined on the Curvilinear basis containing the values
     *      of the advection field at each point on the logical grid.
     *      It is expressed on the contravariant basis.
     * @param [in] advection_field_xy_centre
     *      A vector in the Cartesian basis, containing the value of the advection
     *      field at the O-point.
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvectionRTheta advection_field_rtheta,
            DTensor<CartesianBasis> const& advection_field_xy_centre,
            double dt) const
    {
        using DimX = typename LogicalToPhysicalMapping::cartesian_tag_x;
        using DimY = typename LogicalToPhysicalMapping::cartesian_tag_y;
        using IdxRangeBatchedWithoutR = ddc::remove_dims_of_t<IdxRangeBatched, GridR>;
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeBatched grid(get_idx_range(allfdistribu));
        IdxRangeR radial_grid(grid);
        IdxRangeBatchedWithoutR no_r_grid(grid);

        // Check the first points on R correspond to the O-point.
        CoordXY o_point = m_logical_to_physical_mapping.o_point();
        CoordRTheta point_on_first_row
                = ddc::coordinate(Idx<GridR, GridTheta>(radial_grid.front(), no_r_grid.front()));
        CoordXY diff_points = m_logical_to_physical_mapping(point_on_first_row) - o_point;
        assert(abs(ddc::get<DimX>(diff_points)) < 1e-15);
        assert(abs(ddc::get<DimY>(diff_points)) < 1e-15);

        IdxRangeBatched const
                grid_without_Opoint(radial_grid.remove_first(IdxStep<GridR>(1)), no_r_grid);
        IdxRangeBatched const Opoint_grid(radial_grid.take_first(IdxStepR(1)), no_r_grid);


        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemAdvectionXY advection_field_xy_alloc(grid);
        DVectorFieldAdvectionXY advection_field_xy = get_field(advection_field_xy_alloc);

        LogicalToPhysicalMapping const& logical_to_physical_mapping_proxy
                = m_logical_to_physical_mapping;

        // (Ax, Ay) = J (Ar, Atheta)
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                grid_without_Opoint,
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    IdxRTheta irtheta(idx);
                    CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

                    Tensor J = logical_to_physical_mapping_proxy.jacobian_matrix(coord_rtheta);

                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            idx,
                            tensor_mul(
                                    index<'i', 'j'>(J),
                                    index<'j'>(advection_field_rtheta(idx))));
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                Opoint_grid,
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            idx,
                            advection_field_xy_centre);
                });

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        Kokkos::Profiling::popRegion();

        return allfdistribu;
    }

    /**
     * @brief Allocate a Field to the advected function.
     * 
     * The value at the O-point of the given advection field is not used here. 
     * We compute the advection field on the physical axis at the O-point by 
     * averaging its values at the next interpolation point along r. 
     *
     * @param [in, out] allfdistribu
     *      A Field containing the values of the function we want to advect.
     * @param [in] advection_field_rtheta
     *      A field of vectors defined on the Curvilinear basis containing the values
     *      of the advection field at each point on the logical grid.
     *      It is expressed on the contravariant basis.
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvectionRTheta advection_field_rtheta,
            double dt) const
    {
        using DimX = typename LogicalToPhysicalMapping::cartesian_tag_x;
        using DimY = typename LogicalToPhysicalMapping::cartesian_tag_y;
        using IdxRangeBatchedWithoutR = ddc::remove_dims_of_t<IdxRangeBatched, GridR>;
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeBatched grid(get_idx_range(allfdistribu));
        IdxRangeR radial_grid(grid);
        IdxRangeBatchedWithoutR no_r_grid(grid);

        // Check the first points on R correspond to the O-point.
        CoordXY o_point = m_logical_to_physical_mapping.o_point();
        CoordRTheta point_on_first_row
                = ddc::coordinate(Idx<GridR, GridTheta>(radial_grid.front(), no_r_grid.front()));
        CoordXY diff_points = m_logical_to_physical_mapping(point_on_first_row) - o_point;
        assert(abs(ddc::get<DimX>(diff_points)) < 1e-15);
        assert(abs(ddc::get<DimY>(diff_points)) < 1e-15);

        IdxRangeBatched const
                grid_without_Opoint(radial_grid.remove_first(IdxStep<GridR>(1)), no_r_grid);
        IdxRangeBatched const Opoint_grid(radial_grid.take_first(IdxStepR(1)), no_r_grid);
        IdxRangeBatched const grid_with_points_justafter_Opoint(
                radial_grid.take_first(IdxStepR(2)).remove_first(IdxStep<GridR>(1)),
                no_r_grid);


        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemAdvectionXY advection_field_xy_alloc(grid);
        DVectorFieldAdvectionXY advection_field_xy = get_field(advection_field_xy_alloc);

        LogicalToPhysicalMapping const& logical_to_physical_mapping_proxy
                = m_logical_to_physical_mapping;

        // (Ax, Ay) = J (Ar, Atheta)
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                grid_without_Opoint,
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    IdxRTheta irtheta(idx);
                    CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

                    Tensor J = logical_to_physical_mapping_proxy.jacobian_matrix(coord_rtheta);

                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            idx,
                            tensor_mul(
                                    index<'i', 'j'>(J),
                                    index<'j'>(advection_field_rtheta(idx))));
                });

        // Jacobian ill-defined at the O-point, we average the values around the O-point,
        double advection_field_xy_average_centre_sumx = ddc::parallel_transform_reduce(
                Kokkos::DefaultExecutionSpace(),
                grid_with_points_justafter_Opoint,
                0.,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    return ddcHelper::get<DimX>(advection_field_xy(idx));
                });
        double advection_field_xy_average_centre_sumy = ddc::parallel_transform_reduce(
                Kokkos::DefaultExecutionSpace(),
                grid_with_points_justafter_Opoint,
                0.,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    return ddcHelper::get<DimY>(advection_field_xy(idx));
                });
        DTensor<CartesianBasis> advection_field_xy_average_centre {
                1. / grid_with_points_justafter_Opoint.size()
                * CoordXY(
                        advection_field_xy_average_centre_sumx,
                        advection_field_xy_average_centre_sumy)};

        // and assign the averaged value to all the points at the O-point.
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                Opoint_grid,
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            idx,
                            advection_field_xy_average_centre);
                });

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        Kokkos::Profiling::popRegion();

        return allfdistribu;
    }
};
