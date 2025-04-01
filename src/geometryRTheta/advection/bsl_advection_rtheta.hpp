// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "i_interpolator_2d.hpp"
#include "indexed_tensor.hpp"
#include "metric_tensor_evaluator.hpp"
#include "spline_interpolator_2d.hpp"
#include "spline_polar_foot_finder.hpp"
#include "vector_field.hpp"
#include "vector_field_mem.hpp"



/**
 * @brief Define an advection operator on 2D @f$(r, \theta)@f$ index range.
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
    using CoordRTheta = typename LogicalToPhysicalMapping::CoordArg;
    using CoordXY = typename LogicalToPhysicalMapping::CoordResult;

    using CartesianBasis = ddc::to_type_seq_t<CoordXY>;
    using CurvilinearBasis = ddc::to_type_seq_t<CoordRTheta>;

    using IdxRangeBatched = typename FootFinder::IdxRangeOperator;

    using IdxRangeRTheta = IdxRange<
            find_grid_t<
                    ddc::type_seq_element_t<0, CurvilinearBasis>,
                    ddc::to_type_seq_t<IdxRangeBatched>>,
            find_grid_t<
                    ddc::type_seq_element_t<1, CurvilinearBasis>,
                    ddc::to_type_seq_t<IdxRangeBatched>>>;

    using IdxRTheta = typename IdxRangeRTheta::discrete_element_type;
    using IdxStepRTheta = typename IdxRangeRTheta::discrete_vector_type;

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
        std::unique_ptr<IInterpolator2D<IdxRangeRTheta, IdxRangeRTheta>> const interpolator_ptr
                = m_interpolator.preallocate();

        // Initialise the feet
        CFieldMemFeetRTheta feet_rtheta_alloc(get_idx_range(advection_field_xy));
        CFieldFeetRTheta feet_rtheta = get_field(feet_rtheta_alloc);
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                get_idx_range(advection_field_xy),
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    feet_rtheta(irtheta) = ddc::coordinate(irtheta);
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
     *      A CoordXY containing the value of the advection field on the 
     *      physical index range axis at the O-point. 
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvectionRTheta advection_field_rtheta,
            CoordXY const& advection_field_xy_centre,
            double dt) const
    {
        Kokkos::Profiling::pushRegion("PolarAdvection");
        IdxRangeRTheta grid(get_idx_range(allfdistribu));

        const int npoints_theta = IdxRangeTheta(grid).size();
        IdxRangeRTheta const grid_without_Opoint(grid.remove_first(IdxStepRTheta(1, 0)));
        IdxRangeRTheta const Opoint_grid(grid.take_first(IdxStepRTheta(1, npoints_theta)));


        // Convert advection field on RTheta to advection field on XY
        DVectorFieldMemAdvectionXY advection_field_xy_alloc(grid);
        DVectorFieldAdvectionXY advection_field_xy = get_field(advection_field_xy_alloc);

        LogicalToPhysicalMapping const& logical_to_physical_mapping_proxy
                = m_logical_to_physical_mapping;

        // (Ax, Ay) = J (Ar, Atheta)
        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                grid_without_Opoint,
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    CoordRTheta const coord_rtheta(ddc::coordinate(irtheta));

                    Tensor J = logical_to_physical_mapping_proxy.jacobian_matrix(coord_rtheta);

                    DTensor<CartesianBasis> advec_field_xy = tensor_mul(
                            index<'i', 'j'>(J),
                            index<'j'>(advection_field_rtheta(irtheta)));
                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            irtheta,
                            advec_field_xy);
                });

        ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                Opoint_grid,
                KOKKOS_LAMBDA(IdxRTheta const irtheta) {
                    ddcHelper::assign_vector_field_element(
                            advection_field_xy,
                            irtheta,
                            advection_field_xy_centre);
                });

        (*this)(get_field(allfdistribu), get_const_field(advection_field_xy), dt);

        Kokkos::Profiling::popRegion();

        return allfdistribu;
    }
};
