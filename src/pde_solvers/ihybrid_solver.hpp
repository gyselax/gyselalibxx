// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "vector_field.hpp"

template <class IdxRangeHybrid, class IdxRangeFull, class MemorySpace, class LayoutSpace>
class IHybridSolver;

/**
 * An abstract class from which a Poisson solver can inherit.
 * Classes inheriting from this must implement a way to solve the following equation:
 * @f$ -\Delta \phi = \rho @f$
 *
 * @tparam IdxRangeLaplacian The index range on which the equation is defined.
 * @tparam IdxRangeFull The index range on which the operator() acts. This is equal to the
 *                      IdxRangeLaplacian plus any batched dimensions.
 * @tparam LayoutSpace The layout space of the Fields passed to operator().
 * @tparam MemorySpace The space (CPU/GPU) where the Fields passed to operator()
 *                      are saved.
 */
template <class... ODims, class IdxRangeFull, class MemorySpace, class LayoutSpace>
class IHybridSolver<IdxRange<ODims...>, IdxRangeFull, MemorySpace, LayoutSpace>
{
protected:
    /// @brief The tags describing the real dimensions in the equation.
    using real_hybrid_tags = ddc::detail::TypeSeq<typename ODims::continuous_dimension_type...>;
    /// @brief The tags describing the discrete dimensions in the equation.
    using hybrid_tags = ddc::detail::TypeSeq<ODims...>;
    /// @brief The tags describing the dimensions of the index range on which the operator acts.
    using space_tags = ddc::to_type_seq_t<IdxRangeFull>;
    /// @brief The tags describing the batched dimensions.
    using batch_tags = ddc::type_seq_remove_t<space_tags, hybrid_tags>;

protected:
    /// @brief Indicates whether the gradient is represented by a VectorField or a Field.
    static constexpr bool using_vector_field = ddc::type_seq_size_v<hybrid_tags> == 1;

public:
    /// @brief The Field type of the arguments to operator().
    using field_type = DField<IdxRangeFull, MemorySpace, LayoutSpace>;
    /// @brief The const Field type of the arguments to operator().
    using const_field_type = DConstField<IdxRangeFull, MemorySpace, LayoutSpace>;

    /// @brief The type of the derivative of @f$ \phi @f$.
    using vector_field_type = std::conditional_t<
            ddc::type_seq_size_v<hybrid_tags> == 1,
            field_type,
            VectorField<double, IdxRangeFull, real_hybrid_tags, MemorySpace, LayoutSpace>>;

    /// @brief The index range type describing the batch dimensions.
    using batch_idx_range_type =
            typename ddc::detail::convert_type_seq_to_discrete_domain_t<batch_tags>;
    /// @brief The index for indexing a batch dimension.
    using batch_index_type = typename batch_idx_range_type::discrete_element_type;

    /// @brief The type of the index range on which the equation is defined.
    using hybrid_idx_range_type = IdxRange<ODims...>;

    /// @brief The layout space of the Fields passed to operator().
    using layout_space = LayoutSpace;
    /// @brief The space (CPU/GPU) where the Fields passed to operator() are saved.
    using memory_space = MemorySpace;

public:
    /**
     * @brief An operator which calculates the solution @f$\phi@f$ to Poisson's equation:
     * @f$ -\Delta \phi = \rho @f$
     *
     * @param[out] phi The solution to Poisson's equation.
     * @param[in] rho The right-hand side of Poisson's equation.
     *
     * @return A reference to the solution to Poisson's equation.
     */
    virtual field_type operator()(field_type phi, field_type rho) const = 0;

    /**
     * @brief An operator which calculates the solution @f$\phi@f$ to Poisson's equation and
     * its derivative:
     * @f$ - \Delta \phi = \rho @f$
     * @f$ E = - \nabla \phi @f$
     *
     * @param[out] phi The solution to Poisson's equation.
     * @param[out] E The derivative of the solution to Poisson's equation.
     * @param[in] rho The right-hand side of Poisson's equation.
     *
     * @return A reference to the solution to Poisson's equation.
     */
    virtual field_type operator()(field_type phi, vector_field_type E, field_type rho) const = 0;

    /**
     * @brief An operator which calculates the solution @f$ magnetic_field_z, pressure_field_z@f$ to sub-step pvb 
     * in the hybrid model
     *
     * @param[out] pressure_field The pressure, the result of the nonlinear solver.
     * @param[in] pressure_field_old The pressure in last time step, the input of the nonlinear solver.
     * @param[in] pressure_field_mid The pressure in middle time step, the intermidiate value of the nonlinear solver.
     * @param[in] pressure_field_previous The pressure in previous iteration, the intermidiate value of the nonlinear solver.
     * @param[out] magnetic_field_z The magnetic field, the result of the nonlinear solver.
     * @param[in] magnetic_field_z_old The magnetic field in last time step, the input of the nonlinear solver.
     * @param[in] magnetic_field_z_mid The magnetic field in middle time step, the intermidiate value of the nonlinear solver.
     * @param[in] mean_velocity_x_mid The intermidiate value of the nonlinear solver, the average of ux during a time step.
     * @param[in] mean_velocity_y_mid The intermidiate value of the nonlinear solver, the average of uy during a time step.
     * @param[out] rhs_1 the intermidiate value of the nonlinear solver, and the velocity frame shift in vx.
     * @param[out] rhs_2 the intermidiate value of the nonlinear solver, and the velocity frame shift in vy.
     * @param[in] rhs_3 the intermidiate value of the nonlinear solver.
     * @param[in] rhs_4 the intermidiate value of the nonlinear solver.
     * @param[in] mean_velocity_x The mean velocity in vx in the last time step.
     * @param[in] mean_velocity_y The mean velocity in vy in the last time step.
     * @param[in] mean_velocity_x_each The mean velocity in vx in the last time step for each species ions.
     * @param[in] mean_velocity_x_each The mean velocity in vx in the last time step for each species ions.
     * @param[in] rho_each The charge density in the last time step for each species ions.
     * @param[in] gradx_magnetic The derivative along x of the magnetic field, the intermidiate value of the nonlinear solver.
     * @param[in] grady_magnetic The derivative along y of the magnetic field, the intermidiate value of the nonlinear solver.
     * @param[in] gradx_pressure The derivative along x of the pressure field, the intermidiate value of the nonlinear solver.
     * @param[in] grady_pressure The derivative along y of the pressure field, the intermidiate value of the nonlinear solver.
     * @param[in] rho The charge density for all species ions.
     * @param[in] dt The time step.
     *
     * @return A reference to the solution of sub-step pvb.
     */
    virtual field_type operator()(field_type pressure_field,
        field_type pressure_field_old,
        field_type pressure_field_mid,
        field_type pressure_field_previous,
        field_type magnetic_field_z,
        field_type magnetic_field_z_old,
        field_type magnetic_field_z_mid,
        field_type magnetic_field_z_previous,
        field_type mean_velocity_x_mid,
        field_type mean_velocity_y_mid,
        field_type rhs_1,
        field_type rhs_2,
        field_type rhs_3,
        field_type rhs_4,
        field_type mean_velocity_x,
        field_type mean_velocity_y,
        DFieldSpXY mean_velocity_x_each,
        DFieldSpXY mean_velocity_y_each,
        DFieldSpXY rho_each,
        field_type gradx_magnetic,
        field_type grady_magnetic,
        field_type gradx_pressure,
        field_type grady_pressure,
        field_type rho,
        double dt) const = 0;
};
