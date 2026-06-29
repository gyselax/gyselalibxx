// SPDX-License-Identifier: MIT
#pragma once

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "indexed_tensor.hpp"
#include "l_norm_tools.hpp"
#include "metric_tensor_evaluator.hpp"
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
 * which are determined in the PolarFootFinder operator.
 *
 * The interpolation of the function is always done in the logical domain,
 * where the B-splines are defined. 
 *
 *
 * @see IPolarFootFinder
 */
template <class FootFinder, class LogicalToPhysicalMapping, class Builder2D, class Evaluator2D>
class BslAdvectionPolar
{
    using R = typename LogicalToPhysicalMapping::curvilinear_tag_r;
    using Theta = typename LogicalToPhysicalMapping::curvilinear_tag_theta;

    using DimX = typename LogicalToPhysicalMapping::cartesian_tag_x;
    using DimY = typename LogicalToPhysicalMapping::cartesian_tag_y;

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
    using IdxR = Idx<GridR>;
    using IdxTheta = Idx<GridTheta>;
    using IdxBatched = typename IdxRangeBatched::discrete_element_type;
    using IdxBatch = typename IdxRangeBatch::discrete_element_type;
    using IdxStepR = typename IdxRangeR::discrete_vector_type;

    using MemorySpace = typename FootFinder::memory_space;
    using ExecSpace = typename FootFinder::ExecSpace;

    using IdxRangeBSRTheta =
            typename Builder2D::template batched_spline_domain_type<IdxRangeBatched>;

    using DFieldFDistribu = DField<IdxRangeBatched, MemorySpace>;

    using CFieldMemFeetRTheta = FieldMem<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using CFieldFeetRTheta = Field<CoordRTheta, IdxRangeBatched, MemorySpace>;

    using DVectorFieldMemAdvectionXY
            = DVectorFieldMem<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXY = DVectorField<IdxRangeBatched, CartesianBasis, MemorySpace>;
    using DVectorConstFieldAdvectionXY
            = DVectorConstField<IdxRangeBatched, CartesianBasis, MemorySpace>;

    using DVectorConstFieldAdvection = DVectorConstField<
            IdxRangeBatched,
            VectorIndexSet<typename FootFinder::AdvDim1, typename FootFinder::AdvDim2>,
            MemorySpace>;

    using DVectorFieldMemAdvectionXYOnBatch
            = DVectorFieldMem<IdxRangeBatch, CartesianBasis, MemorySpace>;
    using DVectorFieldAdvectionXYOnBatch = DVectorField<IdxRangeBatch, CartesianBasis, MemorySpace>;

    using DVectorFieldAdvectionRTheta
            = DVectorField<IdxRangeBatched, CurvilinearBasis, MemorySpace>;
    using DVectorConstFieldAdvectionRTheta
            = DVectorConstField<IdxRangeBatched, CurvilinearBasis, MemorySpace>;

private:
    Builder2D const& m_builder_2d;

    Evaluator2D const& m_evaluator_2d;

    FootFinder& m_find_feet_method;

    LogicalToPhysicalMapping const& m_logical_to_physical_mapping;

    std::optional<IdxRangeBatched> m_idx_range_advected_points;

public:
    /**
     * @brief Instantiate an advection operator.
     *
     * @param[in] builder_2d
     *      The 2D builder used to compute interpolations coefficients from
     *      the function values at interpolation points.
     * @param[in] evaluator_2d
     *      The 2D evaluator used to evaluate the interpolating function at the
     *      feet of the characteristics.
     * @param[in] foot_finder
     *      An IFootFinder which computes the feet of the characteristics.
     * @param[in] logical_to_physical_mapping
     *      The mapping function from the logical domain to the physical
     *      domain.
     */
    BslAdvectionPolar(
            Builder2D const& builder_2d,
            Evaluator2D const& evaluator_2d,
            FootFinder& foot_finder,
            LogicalToPhysicalMapping const& logical_to_physical_mapping,
            std::optional<IdxRangeBatched> idx_range_advected_points = std::nullopt)
        : m_builder_2d(builder_2d)
        , m_evaluator_2d(evaluator_2d)
        , m_find_feet_method(foot_finder)
        , m_logical_to_physical_mapping(logical_to_physical_mapping)
        , m_idx_range_advected_points(idx_range_advected_points)
    {
    }

    ~BslAdvectionPolar() = default;


    /**
     * @brief Advect a function over a time step dt with the given advection field
     * along the physical directions. 
     *
     * @param [in, out] allfdistribu
     *      A Field containing the values of the function we want to advect.
     * @param [in] advection_field
     *      A field of vectors defined on the Cartesian basis containing the values
     *      of the advection field at each point on the logical grid.
     * @param [in] dt
     *      A time step used.
     *
     * @return A Field to allfdistribu advected on the time step given.
     */
    DFieldFDistribu operator()(
            DFieldFDistribu allfdistribu,
            DVectorConstFieldAdvection advection_field,
            double dt) const
    {
        // Pre-allocate spline coefficient storage
        DFieldMem<IdxRangeBSRTheta, MemorySpace> coefs_alloc(
                m_builder_2d.batched_spline_domain(get_idx_range(allfdistribu)));

        // Compute the feet of the characteristics at tn -----------------------------------------
        typename FootFinder::ElementwiseOperator find_foot_alloc
                = m_find_feet_method(get_const_field(advection_field));
        typename FootFinder::ElementwiseOperator::GPUCompat find_foot = find_foot_alloc(dt);

        // Interpolate the function on the feet of the characteristics. --------------------------
        Kokkos::Profiling::pushRegion("(GSLX) BslAdvectionPolar/Interpolation");
        m_builder_2d(get_field(coefs_alloc), get_const_field(allfdistribu));

        Evaluator2D const& evaluator_2d_proxy = m_evaluator_2d;

        DConstField<IdxRangeBSRTheta, MemorySpace> coefs = get_const_field(coefs_alloc);

        IdxRangeBatched idx_range_advected_points = m_idx_range_advected_points ? *m_idx_range_advected_points : get_idx_range(allfdistribu);

        const std::source_location location = std::source_location::current();
        ddc::parallel_for_each(
                location.function_name(),
                ExecSpace(),
                idx_range_advected_points,
                KOKKOS_LAMBDA(IdxBatched const idx) {
                    IdxBatch batch_idx(idx);
                    allfdistribu(idx) = evaluator_2d_proxy(find_foot(idx), coefs[batch_idx]);
                });
        Kokkos::Profiling::popRegion();

        unify_value_at_centre_pt(allfdistribu);

        return allfdistribu;
    }


    /**
     * @brief Replace the value at @f$  (r=0, \theta)@f$  point
     *  by the value at @f$ (r=0,0) @f$ for all @f$ \theta @f$.
     *
     *  For polar geometry, to ensure continuity at the centre point, we
     *  have to be sure that all the points for @f$ r = 0 @f$ have the same value.
     *  As the computation of the values of a table can induces machine errors,
     *  this function is useful to reset the values at the central point at
     *  the same value.
     *
     *  @param[in, out] values
     *      The table of values we want to unify at the central point.
     */
    template <class T>
    static void unify_value_at_centre_pt(Field<T, IdxRangeBatched, MemorySpace> values)
    {
        IdxRangeBatched full_idx_range = get_idx_range(values);
        IdxRangeBatch const batched_idx_range(full_idx_range);
        IdxRangeR const r_idx_range(full_idx_range);
        IdxRangeTheta const theta_idx_range(full_idx_range);
        IdxR r0_idx = r_idx_range.front();
        IdxTheta theta0_idx = theta_idx_range.front();
        if (std::fabs(ddc::coordinate(r0_idx)) < 1e-15) {
            const std::source_location location = std::source_location::current();
            ddc::parallel_for_each(
                    location.function_name(),
                    ExecSpace(),
                    batched_idx_range,
                    KOKKOS_LAMBDA(const IdxBatch ib) {
                        for (IdxTheta itheta : theta_idx_range) {
                            values(ib, r0_idx, itheta) = values(ib, r0_idx, theta0_idx);
                        }
                    });
        }
    }
};
