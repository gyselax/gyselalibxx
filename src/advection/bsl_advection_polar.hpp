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

    using DimX = typename LogicalToPhysicalMapping::cartesian_tag_x;
    using DimY = typename LogicalToPhysicalMapping::cartesian_tag_y;

    using XYBasis = VectorIndexSet<DimX, DimY>;

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
    using IdxRangeBatch = ddc::remove_dims_of_t<IdxRangeBatched, GridR, GridTheta>;

    using IdxRTheta = typename IdxRangeRTheta::discrete_element_type;
    using IdxBatched = typename IdxRangeBatched::discrete_element_type;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxStepR = typename IdxRangeR::discrete_vector_type;

    using MemorySpace = typename FootFinder::memory_space;
    using ExecSpace = typename FootFinder::ExecSpace;

    using DFieldFDistribu = DField<IdxRangeBatched, MemorySpace>;

    using CFieldMemFeetRTheta = FieldMem<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using CFieldFeetRTheta = Field<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using DVectorFieldMemAdvectionXY
            = DVectorFieldMem<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXY = DVectorField<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorConstFieldAdvectionXY
            = DVectorConstField<IdxRangeBatched, CartesianBasis, MemorySpace>;

    using DVectorFieldMemAdvectionXYOnBatch
            = DVectorFieldMem<IdxRangeBatch, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXYOnBatch = DVectorField<IdxRangeBatch, CartesianBasis, MemorySpace>;

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
                ExecSpace(),
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
     * @warning This operator should be applied if the O-point corresponds to 
     * points of the grid. 
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
        using IdxRangeBatchedWithoutR = ddc::remove_dims_of_t<IdxRangeBatched, GridR>;
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeBatched grid(get_idx_range(allfdistribu));
        IdxRangeR radial_grid(grid);
        IdxRangeBatchedWithoutR no_r_grid(grid);

        // Check the first points on R correspond to the O-point.
        assert(norm_inf(
                       m_logical_to_physical_mapping(ddc::coordinate(
                               Idx<GridR, GridTheta>(radial_grid.front(), no_r_grid.front())))
                       - m_logical_to_physical_mapping.o_point())
               < 1e-13);

        IdxRangeBatched const
                grid_without_Opoint(radial_grid.remove_first(IdxStep<GridR>(1)), no_r_grid);
        IdxRangeBatched const Opoint_grid(radial_grid.take_first(IdxStepR(1)), no_r_grid);


        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemAdvectionXY advection_field_xy_alloc(grid);
        DVectorFieldAdvectionXY advection_field_xy = get_field(advection_field_xy_alloc);

        LogicalToPhysicalMapping const& logical_to_physical_mapping_proxy
                = m_logical_to_physical_mapping;

        // (Ax, Ay) = J (Ar, Atheta)
        copy_to_vector_space<XYBasis>(
                ExecSpace(),
                advection_field_xy[grid_without_Opoint],
                logical_to_physical_mapping_proxy,
                advection_field_rtheta[grid_without_Opoint]);

        ddc::parallel_for_each(
                ExecSpace(),
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
        using IdxRangeBatchedWithoutR = ddc::remove_dims_of_t<IdxRangeBatched, GridR>;
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeBatched grid(get_idx_range(allfdistribu));
        IdxRangeR radial_grid(grid);
        IdxRangeTheta theta_grid(grid);
        IdxRangeBatchedWithoutR no_r_grid(grid);
        IdxRangeBatch no_rtheta_grid(grid);

        LogicalToPhysicalMapping const& logical_to_physical_mapping_proxy
                = m_logical_to_physical_mapping;

        // Test if the first points on R correspond to the O-point.
        CoordXY const Opoint = m_logical_to_physical_mapping.o_point();
        CoordRTheta const point_on_first_row
                = ddc::coordinate(Idx<GridR, GridTheta>(radial_grid.front(), no_r_grid.front()));
        CoordXY const diff_points = m_logical_to_physical_mapping(point_on_first_row) - Opoint;
        bool const first_row_is_o_point
                = ((abs(ddc::get<DimX>(diff_points)) < 1e-15)
                   && (abs(ddc::get<DimY>(diff_points)) < 1e-15));

        /*
            If the O-point is not considered as on the grid, we have to be sure that the grid of 
            the advection field on (R, Theta) matches with the grid of the advected function.  
        */
        assert(first_row_is_o_point
               || (IdxRangeR(get_idx_range(advection_field_rtheta)) == radial_grid));

        IdxRangeBatched const grid_without_Opoint(
                radial_grid.remove_first(IdxStep<GridR>(first_row_is_o_point)),
                no_r_grid);

        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemAdvectionXY advection_field_xy_alloc(grid);
        DVectorFieldAdvectionXY advection_field_xy = get_field(advection_field_xy_alloc);

        // (Ax, Ay) = J (Ar, Atheta)
        copy_to_vector_space<XYBasis>(
                ExecSpace(),
                advection_field_xy[grid_without_Opoint],
                logical_to_physical_mapping_proxy,
                advection_field_rtheta[grid_without_Opoint]);

        // Treatment for the O-point.
        if (first_row_is_o_point) {
            IdxRangeBatched const Opoint_grid(radial_grid.take_first(IdxStepR(1)), no_r_grid);
            IdxRangeRTheta const grid_first_ring(
                    radial_grid.take_first(IdxStepR(2)).remove_first(IdxStep<GridR>(1)),
                    theta_grid);

            DVectorFieldMemAdvectionXYOnBatch advection_field_xy_average_centre_alloc(
                    no_rtheta_grid);
            DVectorFieldAdvectionXYOnBatch advection_field_xy_average_centre(
                    advection_field_xy_average_centre_alloc);

            // Jacobian ill-defined at the O-point, we average the values around the O-point,
            std::size_t ntheta_points = theta_grid.size();
            ddc::for_each(no_rtheta_grid, [&](IdxBatch const idx_batch) {
                CoordXY advection_field_xy_average_on_theta = average_field(
                        get_const_field(ddcHelper::get<DimX>(advection_field_xy)[idx_batch]),
                        get_const_field(ddcHelper::get<DimY>(advection_field_xy)[idx_batch]),
                        grid_first_ring,
                        ntheta_points);

                DTensor<CartesianBasis> advection_field_xy_average_on_theta_tensor(
                        advection_field_xy_average_on_theta);

                ddcHelper::assign_vector_field_element(
                        advection_field_xy_average_centre,
                        idx_batch,
                        advection_field_xy_average_on_theta_tensor);
            });

            // and assign the averaged value to all the points at the O-point.
            ddc::parallel_for_each(
                    ExecSpace(),
                    Opoint_grid,
                    KOKKOS_LAMBDA(IdxBatched const idx) {
                        IdxBatch idx_batch(idx);
                        ddcHelper::get<DimX>(advection_field_xy)(idx) = ddcHelper::get<DimX>(
                                advection_field_xy_average_centre)(idx_batch);
                        ddcHelper::get<DimY>(advection_field_xy)(idx) = ddcHelper::get<DimY>(
                                advection_field_xy_average_centre)(idx_batch);
                    });
        }

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        Kokkos::Profiling::popRegion();

        return allfdistribu;
    }


    /**
     * @brief Average a field on X and a field on Y defined on (R,Theta) 
     * and return a coordinate on (X,Y) of the averaged values. 
     * @param advection_field_x Advection field along X defined on the cross-section (R,Theta).
     * @param advection_field_y Advection field along Y defined on the cross-section (R,Theta).
     * @param grid_first_ring Index range containing the IdxR(1) and all the Theta indices, 
     * @param ntheta_points Number of points on the theta grid. 
     * @return A coordinate on (X,Y) ccontaining (sum_theta A_x(r, theta), sum_theta A_y(r, theta)).
     */
    CoordXY average_field(
            DConstField<IdxRangeRTheta, MemorySpace> advection_field_x,
            DConstField<IdxRangeRTheta, MemorySpace> advection_field_y,
            IdxRangeRTheta grid_first_ring,
            std::size_t ntheta_points) const
    {
        double const sum_x = ddc::parallel_transform_reduce(
                ExecSpace(),
                grid_first_ring,
                0.,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(IdxRTheta const idx_rtheta) {
                    return advection_field_x(idx_rtheta);
                });
        double const sum_y = ddc::parallel_transform_reduce(
                ExecSpace(),
                grid_first_ring,
                0.,
                ddc::reducer::sum<double>(),
                KOKKOS_LAMBDA(IdxRTheta const idx_rtheta) {
                    return advection_field_y(idx_rtheta);
                });
        return 1. / ntheta_points * CoordXY(sum_x, sum_y);
    }
};
