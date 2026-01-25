// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "vector_field.hpp"
#include "geometry.hpp"

template <class IdxRangeHybrid, class IdxRangeFull, class MemorySpace, class LayoutSpace>
class IHybridSolver1d;

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
class IHybridSolver1d<IdxRange<ODims...>, IdxRangeFull, MemorySpace, LayoutSpace>
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
     * @param[in] dt The time step.
     *
     * @return A reference to the solution of sub-step pvb.
     */

    
    virtual field_type operator()(
        field_type magnetic_field_y,
        field_type magnetic_field_y_old,
        field_type magnetic_field_y_mid, 
        field_type magnetic_field_y_previous,
        field_type magnetic_field_z, 
        field_type magnetic_field_z_old, 
        field_type magnetic_field_z_mid, 
        field_type magnetic_field_z_previous,
        field_type rho,
        field_type magnetic_field_x,
        field_type gradx_magnetic_field_y_mid,
        field_type gradx_magnetic_field_z_mid,
        double dt) const = 0;


    virtual field_type operator()(
        field_type magnetic_field_y, field_type magnetic_field_y_old,
        field_type magnetic_field_y_mid, field_type magnetic_field_y_previous,
        field_type magnetic_field_z, field_type magnetic_field_z_old, 
        field_type magnetic_field_z_mid, field_type magnetic_field_z_previous,
        DFieldSpX u_old_x, DFieldSpX u_old_y, DFieldSpX u_old_z, 
        field_type u_bar_x, field_type u_bar_y, field_type u_bar_z, 
        field_type rho, DFieldSpX rho_each,
        field_type magnetic_field_x,
        field_type gradx_rho,
        field_type gradx_magnetic_field_y_mid,
        field_type gradx_magnetic_field_z_mid,
        field_type rhs_1,
        field_type rhs_2,
        field_type rhs_3,
        field_type rhs_5,
        field_type rhs_6,
        field_type Mxx, field_type Mxy, field_type Mxz,
        field_type Myx, field_type Myy, field_type Myz,
        field_type Mzx, field_type Mzy, field_type Mzz,
        field_type weighted_u_x, field_type weighted_u_y, field_type weighted_u_z,
        field_type weighted_p_para_x, field_type weighted_p_para_y, field_type weighted_p_para_z,
        field_type qx, field_type qy, field_type qz,
        field_type p_parallel_x, field_type p_parallel_y, field_type p_parallel_z,
        double const electron_temperature,
        double dt) const = 0;
};
