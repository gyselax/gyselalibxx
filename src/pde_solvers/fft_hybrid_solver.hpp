// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "ihybrid_solver.hpp"
#include "vector_index_tools.hpp"
#include "geometry.hpp"
#include "l_norm_tools.hpp"
#include <iostream>
#include <cmath>
#include "species_info.hpp"

/**
 * See @ref FFTPoissonSolverImplementation.
 */
template <
        class IdxRangeHybrid,
        class IdxRangeFull = IdxRangeHybrid,
        class ExecSpace = Kokkos::DefaultExecutionSpace,
        class LayoutSpace = Kokkos::layout_right>
class FFTHybridSolver;

/**
 * @brief A class to solve the following equation:
 * @f$ -\Delta \phi = \rho @f$
 * using a Fourier transform.
 *
 * The implementation of this class can be found at FFTPoissonSolver< IdxRange<GridPDEDim1D...>, IdxRangeFull, ExecSpace, LayoutSpace >.
 * @anchor FFTPoissonSolverImplementation
 *
 * @tparam IdxRangeLaplacian The index range on which the equation is defined.
 * @tparam IdxRangeFull The index range on which the operator() acts. This is equal to the
 *                      IdxRangeLaplacian plus any batched dimensions.
 * @tparam ExecSpace The space (CPU/GPU) where the calculations will take place.
 * @tparam LayoutSpace The layout space of the Fields passed to operator().
 */
template <class... GridPDEDim1D, class IdxRangeFull, class ExecSpace, class LayoutSpace>
class FFTHybridSolver<IdxRange<GridPDEDim1D...>, IdxRangeFull, ExecSpace, LayoutSpace>
    : public IHybridSolver<
              IdxRange<GridPDEDim1D...>,
              IdxRangeFull,
              typename ExecSpace::memory_space,
              LayoutSpace>
{
private:
    using base_type = IHybridSolver<
            IdxRange<GridPDEDim1D...>,
            IdxRangeFull,
            typename ExecSpace::memory_space,
            LayoutSpace>;

public:
    template <class Dim>
    struct GridFourier : ddc::PeriodicSampling<ddc::Fourier<Dim>>
    {
    };

public:
    /// @brief The Field type of the arguments to operator().
    using field_type = typename base_type::field_type;
    /// @brief The const Field type of the arguments to operator().
    using const_field_type = typename base_type::const_field_type;

    /// @brief The type of the derivative of @f$ \phi @f$.
    using vector_field_type = typename base_type::vector_field_type;

    /// @brief The index range type describing the batch dimensions.
    using batch_idx_range_type = typename base_type::batch_idx_range_type;
    /// @brief The index type for indexing a batch dimension.
    using batch_index_type = typename base_type::batch_index_type;

    /// @brief The type of the index range on which the equation is defined.
    using hybrid_idx_range_type = typename base_type::hybrid_idx_range_type;

    /// @brief The layout space of the Fields passed to operator().
    using layout_space = typename base_type::layout_space;
    /// @brief The space (CPU/GPU) where the Fields passed to operator() are saved.
    using memory_space = typename base_type::memory_space;

    /// @brief The type of the Fourier space index range.
    using fourier_idx_range_type
            = IdxRange<GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>;
    /// @brief The type of an index of the Fourier space index range.
    using fourier_index_type = typename fourier_idx_range_type::discrete_element_type;

    /// @brief The type of a Field storing the Fourier transform of a function.
    using fourier_field_mem_type
            = FieldMem<Kokkos::complex<double>, fourier_idx_range_type, memory_space>;
    /// @brief The type of a Field storing the Fourier transform of a function.
    using fourier_field_type = typename fourier_field_mem_type::span_type;

private:
    /// @brief The normalisation used for the Fourier transform
    static constexpr ddc::FFT_Normalization m_norm = ddc::FFT_Normalization::BACKWARD;

private:
    /**
     * @brief The multiplicative factor corresponding to the Laplace operator @f$ \Delta @f$
     * expressed in Fourier space at a given mode.
     *
     * @param index The index of the Fourier mode.
     */
    template <class... FDim>
    KOKKOS_FUNCTION static double get_laplace_operator(Idx<FDim...> index)
    {
        return (((double)ddc::coordinate(ddc::select<FDim>(index))
                 * (double)ddc::coordinate(ddc::select<FDim>(index)))
                + ...);
    }

    /**
     * @brief Differentiate an expression from its representation in Fourier space by multiplying
     * by -i * k and then converting back to real space.
     *
     * @param[out] derivative The Field where the derivative will be saved.
     * @param[out] fourier_derivative The derivative of the function in Fourier space will be saved.
     * @param[in] values The Field containing the values of the function in Fourier space.
     *
     * @tparam Dim The dimension along which the expression is differentiated.
     */
    template <class Dim>
    void differentiate_and_invert_fourier_values(
            DField<hybrid_idx_range_type, memory_space, LayoutSpace> derivative,
            fourier_field_type fourier_derivative,
            fourier_field_type values) const
    {
        negative_differentiate_equation<Dim>(fourier_derivative, values);
        // Perform the inverse 1D FFT of the solution to deduce the electric field
        ddc::
                ifft(ExecSpace(),
                     get_field(derivative),
                     get_field(fourier_derivative),
                     ddc::kwArgs_fft {m_norm});
    }

    /**
     * @brief Get the gradient of a 1D expression from its representation in Fourier space.
     *
     * @param[out] gradient The Field where the derivative will be saved.
     * @param[out] fourier_derivative The derivative of the function in Fourier space will be saved.
     * @param[in] values The Field containing the values of the function in Fourier space.
     */
    void get_gradient(
            DField<hybrid_idx_range_type, memory_space, LayoutSpace> gradient,
            fourier_field_type fourier_derivative,
            fourier_field_type values) const
    {
        using Dim =
                typename ddc::type_seq_element_t<0, ddc::to_type_seq_t<hybrid_idx_range_type>>::
                        continuous_dimension_type;
        using FourierDim = GridFourier<Dim>;
        differentiate_and_invert_fourier_values<FourierDim>(gradient, fourier_derivative, values);
    }

    /**
     * @brief Get the derivative in a specific direction of a space defined function.
     *
     * @param[out] gradient The Field where the derivative will be saved.
     * @param[out] fourier_derivative The derivative of the function in Fourier space will be saved.
     * @param[in] values The Field containing the values of the function in real space.
     */
    template <class Dim>
    void get_gradient(
            fourier_field_type intermediate_chunk,
            field_type values) const
    {   
        hybrid_idx_range_type idx_range(get_idx_range(values));
        batch_idx_range_type batch_idx_range(get_idx_range(values));

        ddc::for_each(batch_idx_range, [&](batch_index_type ib) {
                // Perform the inverse FFT
                ddc::fft(ExecSpace(), get_field(intermediate_chunk), values[ib], ddc::kwArgs_fft {m_norm});
                /// Derivative at the direction Dim
                differentiate_equation<GridFourier<Dim>>(intermediate_chunk);
                // Perform the inverse FFT
                ddc::ifft(ExecSpace(), values[ib], get_field(intermediate_chunk), ddc::kwArgs_fft {m_norm});
            });
    }


    template <class Dim>
    void get_gradient(
            fourier_field_type intermediate_chunk,
            field_type values_in, field_type values_out) const
    {   
        hybrid_idx_range_type idx_range(get_idx_range(values_in));
        batch_idx_range_type batch_idx_range(get_idx_range(values_in));

        ddc::for_each(batch_idx_range, [&](batch_index_type ib) {
                // Perform the inverse FFT
                ddc::fft(ExecSpace(), get_field(intermediate_chunk), values_in[ib], ddc::kwArgs_fft {m_norm});
                /// Derivative at the direction Dim
                differentiate_equation<GridFourier<Dim>>(intermediate_chunk);
                // Perform the inverse FFT
                ddc::ifft(ExecSpace(), values_out[ib], get_field(intermediate_chunk), ddc::kwArgs_fft {m_norm});
            });
    }

    /**
     * @brief Get the gradient of a multi-dimensional expression from its representation in Fourier space.
     *
     * @param[out] gradient The VectorField where the gradient will be saved.
     * @param[out] fourier_derivative The derivative of the function in Fourier space will be saved.
     * @param[in] values The Field containing the values of the function in Fourier space.
     */
    template <class... Dims>
    void get_gradient(
            VectorField<
                    double,
                    hybrid_idx_range_type,
                    VectorIndexSet<Dims...>,
                    memory_space,
                    layout_space> gradient,
            fourier_field_type fourier_derivative,
            fourier_field_type values) const
    {
        ((differentiate_and_invert_fourier_values<
                 GridFourier<Dims>>(ddcHelper::get<Dims>(gradient), fourier_derivative, values)),
         ...);
    }

    template <class Grid1D>
    void init_fourier_space(IdxRange<Grid1D> idx_range)
    {
        using GridFFT = GridFourier<typename Grid1D::continuous_dimension_type>;
        ddc::init_discrete_space<GridFFT>(ddc::init_fourier_space<GridFFT>(idx_range));
    }

public:
    /**
     * @brief A function to solve the Poisson equation in Fourier space
     * This function should be private. It is not due to the inclusion of a KOKKOS_LAMBDA
     *
     * @param[out] intermediate_chunk The solution to the Poisson equation in Fourier space.
     * @param[in] rho The right-hand side of the Poisson equation in Fourier space.
     */
    template <class Layout>
    void solve_poisson_equation(
            fourier_field_type intermediate_chunk,
            DField<hybrid_idx_range_type, memory_space, Layout> rho) const
    {
        // Compute FFT(rho)
        ddc::fft(ExecSpace(), intermediate_chunk, rho, ddc::kwArgs_fft {m_norm});

        fourier_idx_range_type const k_mesh = get_idx_range(intermediate_chunk);

        // Solve Poisson's equation -\Delta phi = -(\sum_j \partial_j^2) \phi = rho
        //   in Fourier space as -(\sum_j i*k_i * i*k_i) FFT(Phi) = FFT(rho))
        ddc::parallel_for_each(
                ExecSpace(),
                k_mesh,
                KOKKOS_LAMBDA(fourier_index_type const ik) {
                    if (ik != k_mesh.front()) {
                        intermediate_chunk(ik) = intermediate_chunk(ik) / get_laplace_operator(ik);
                    } else {
                        intermediate_chunk(ik) = 0.;
                    }
                });
    }

    /**
     * @brief Differentiate and multiply by -1 an expression in Fourier space by multiplying by -i * k
     * This function should be private. It is not due to the inclusion of a KOKKOS_LAMBDA
     *
     * @param derivative The Field where the derivative will be saved.
     * @param values The Field containing the values of the function in Fourier space.
     *
     * @tparam Dim The dimension along which the expression is differentiated.
     */
    template <class Dim>
    void negative_differentiate_equation(fourier_field_type derivative, fourier_field_type values)
            const
    {
        Kokkos::complex<double> imaginary_unit(0.0, 1.0);
        ddc::parallel_for_each(
                ExecSpace(),
                get_idx_range(values),
                KOKKOS_LAMBDA(fourier_index_type const ik) {
                    Idx<Dim> const ikx = ddc::select<Dim>(ik);
                    derivative(ik) = -imaginary_unit * ddc::coordinate(ikx) * values(ik);
                });
    }

    template <class Dim>
    void differentiate_equation(fourier_field_type values)
            const
    {
        Kokkos::complex<double> imaginary_unit(0.0, 1.0);
        ddc::parallel_for_each(
                ExecSpace(),
                get_idx_range(values),
                KOKKOS_LAMBDA(fourier_index_type const ik) {
                    Idx<Dim> const ikx = ddc::select<Dim>(ik);
                    values(ik) = imaginary_unit * ddc::coordinate(ikx) * values(ik);
                });
    }

public:
    /**
     * @brief A constructor for the FFT Poisson solver.
     * This constructor calls ddc::init_discrete_space so it should only be called once per
     * simulation.
     *
     * @param hybrid_idx_range The index range on which the equation should be solved.
     */
    explicit FFTHybridSolver(hybrid_idx_range_type hybrid_idx_range)
    {
        ((init_fourier_space<GridPDEDim1D>(ddc::select<GridPDEDim1D>(hybrid_idx_range))), ...);
    }

    /**
     * @brief An operator which calculates the solution @f$\phi@f$ to Poisson's equation:
     * @f$ - \Delta \phi = \rho @f$
     *
     * @param[out] phi The solution to Poisson's equation.
     * @param[in] rho The right-hand side of Poisson's equation.
     *
     * @return A reference to the solution to Poisson's equation.
     */
    virtual field_type operator()(field_type phi, field_type rho) const final
    {
        Kokkos::Profiling::pushRegion("FFTPoissonSolver");

        hybrid_idx_range_type idx_range(get_idx_range(phi));
        batch_idx_range_type batch_idx_range(get_idx_range(phi));

        // Build a mesh in the fourier space, for N points
        fourier_idx_range_type const k_mesh = ddc::fourier_mesh<
                GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>(idx_range, false);

        fourier_field_mem_type intermediate_chunk_alloc(k_mesh);
        fourier_field_type intermediate_chunk = get_field(intermediate_chunk_alloc);

        ddc::for_each(batch_idx_range, [&](batch_index_type ib) {
            solve_poisson_equation(intermediate_chunk, rho[ib]);

            // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
            ddc::
                    ifft(ExecSpace(),
                         phi[ib],
                         get_field(intermediate_chunk),
                         ddc::kwArgs_fft {m_norm});
        });

        Kokkos::Profiling::popRegion();
        return phi;
    }

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
    virtual field_type operator()(field_type phi, vector_field_type E, field_type rho) const final
    {
        Kokkos::Profiling::pushRegion("FFTHybridSolver");

        hybrid_idx_range_type idx_range(get_idx_range(phi));
        batch_idx_range_type batch_idx_range(get_idx_range(phi));

        // Build a mesh in the fourier space, for N points
        fourier_idx_range_type const k_mesh = ddc::fourier_mesh<
                GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>(idx_range, false);

        fourier_field_mem_type intermediate_chunk_alloc(k_mesh);
        fourier_field_mem_type fourier_efield_alloc(k_mesh);

        fourier_field_type intermediate_chunk = get_field(intermediate_chunk_alloc);
        fourier_field_type fourier_efield = get_field(fourier_efield_alloc);

        ddc::for_each(batch_idx_range, [&](batch_index_type ib) {
            solve_poisson_equation(intermediate_chunk, rho[ib]);
            get_gradient(E[ib], fourier_efield, intermediate_chunk);

            // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
            ddc::
                    ifft(ExecSpace(),
                         phi[ib],
                         get_field(intermediate_chunk),
                         ddc::kwArgs_fft {m_norm});
        });
        Kokkos::Profiling::popRegion();
        return phi;
    }

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

    // used for the sub-step pvb
    virtual field_type operator()(field_type pressure_field, field_type pressure_field_old,
                                  field_type pressure_field_mid, field_type pressure_field_previous,
                                  field_type magnetic_field_z, field_type magnetic_field_z_old,
                                  field_type magnetic_field_z_mid, field_type magnetic_field_z_previous,
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
                                  double const dt) const final
    {
        Kokkos::Profiling::pushRegion("FFTHybridSolver_SUBSTEP_PVB");

        hybrid_idx_range_type idx_range(get_idx_range(pressure_field));
        batch_idx_range_type batch_idx_range(get_idx_range(pressure_field));

        // Build a mesh in the fourier space, for N points
        fourier_idx_range_type const k_mesh = ddc::fourier_mesh<
                GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>(idx_range, false);

        fourier_field_mem_type intermediate_chunk_alloc(k_mesh);

        fourier_field_type intermediate_chunk = get_field(intermediate_chunk_alloc);
        
        // the charge and mass of each kinds on ions
        
        IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(mean_velocity_x_each);
        
        host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
        host_t<DConstFieldSp> const kinetic_charges_host
                = charges_host[get_idx_range<Species>(mean_velocity_x_each)];

        auto const kinetic_charges_alloc
                = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
        DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);

        host_t<DConstFieldSp> const masses_host = ddc::host_discrete_space<Species>().masses();
        host_t<DConstFieldSp> const kinetic_masses_host
                = masses_host[get_idx_range<Species>(mean_velocity_x_each)];

        auto const kinetic_masses_alloc
                = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_masses_host);
        DConstFieldSp kinetic_masses = get_const_field(kinetic_masses_alloc);
        
        // Before the picard iteration, set field_old = field
        ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                pressure_field_old(ixy) = pressure_field(ixy);
                magnetic_field_z_old(ixy) = magnetic_field_z(ixy);
            });
        // The picard iteration starts here 
        ///////////////////////////////////////////////////////////////////////////////////////////
        int iter = 0;
        int iter_number = 20;
        for (; iter < iter_number; ++iter) {
            // 0th step, compute the p at the mid-time, and store the p value at pressure_field_previous 
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        pressure_field_previous(ixy) = pressure_field(ixy);
                        pressure_field_mid(ixy) = 0.5 * (pressure_field_old(ixy) + pressure_field(ixy));
                        magnetic_field_z_previous(ixy) = magnetic_field_z(ixy);
                        magnetic_field_z_mid(ixy) = 0.5 * (magnetic_field_z_old(ixy) + magnetic_field_z(ixy));
                        //std::cout << "rho is: " << std::setprecision(16) << mean_velocity_y(ixy) << std::endl;
                    });

            // compute mean_velocity_x_mid, mean_velocity_y_mid, and magnetic_field_z_mid
                /// mean_velocity_x_mid
            get_gradient<X>(intermediate_chunk, pressure_field_mid, gradx_pressure);
            get_gradient<Y>(intermediate_chunk, pressure_field_mid, grady_pressure);
            get_gradient<X>(intermediate_chunk, magnetic_field_z_mid, gradx_magnetic);
            get_gradient<Y>(intermediate_chunk, magnetic_field_z_mid, grady_magnetic);

            
            // set rhs 1, 2, 3, 4 as zero
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const ixy) {
                    rhs_1(ixy) = 0.0;
                    rhs_2(ixy) = 1.0;
                    rhs_3(ixy) = mean_velocity_x(ixy);
                    rhs_4(ixy) = mean_velocity_y(ixy);
                });
            
            ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
                ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        
                        double theta = dt * kinetic_charges(isp) / kinetic_masses(isp) * magnetic_field_z_mid(ixy);
                        double first = 2.0 * theta * (Kokkos::sin(0.5*theta) / theta) * (Kokkos::sin(0.5*theta) / theta) * rho_each(isp,ixy)/rho(ixy);
                        double second = -(1.0 - std::sin(theta)/theta) * rho_each(isp,ixy)/rho(ixy);
                        double ui_plus_curlb_minus_qx = mean_velocity_x_each(isp,ixy) + grady_magnetic(ixy)/rho(ixy) + grady_pressure(ixy)/rho(ixy)/magnetic_field_z_mid(ixy) ;
                        double ui_plus_curlb_minus_qy = mean_velocity_y_each(isp,ixy) - gradx_magnetic(ixy)/rho(ixy) - gradx_pressure(ixy)/rho(ixy)/magnetic_field_z_mid(ixy) ;
                        
                        rhs_1(ixy) += first;
                        rhs_2(ixy) += second;
                        rhs_3(ixy) += first * ui_plus_curlb_minus_qy + second * ui_plus_curlb_minus_qx;
                        rhs_4(ixy) += - first * ui_plus_curlb_minus_qx + second * ui_plus_curlb_minus_qy;

                    });
            }); 
            
            
            
            
            
            
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const ixy) {
                    double determinant = rhs_2(ixy) * rhs_2(ixy) + rhs_1(ixy)*rhs_1(ixy);
                    double upper_left = rhs_2(ixy) / determinant;
                    double upper_right = -rhs_1(ixy) / determinant;
                    double low_left = rhs_1(ixy) / determinant;
                    double low_right = rhs_2(ixy) / determinant;
                    mean_velocity_x_mid(ixy) = upper_left * rhs_3(ixy) + upper_right * rhs_4(ixy);
                    mean_velocity_y_mid(ixy) = low_left * rhs_3(ixy) + low_right * rhs_4(ixy);

                });
            
            
            
            
            
            
            
            /*
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const ixy) {
                    double B3_mid = magnetic_field_z_mid(ixy);
                    double theta = dt * B3_mid;
                    double second = 2.0 * theta * (Kokkos::sin(0.5*theta) / theta) * (Kokkos::sin(0.5*theta) / theta);
                    double third = 1.0 - Kokkos::sin(theta)/theta;
                    double curlb_minus_q_x = grady_magnetic(ixy) / rho(ixy) + grady_pressure(ixy) / rho(ixy) / B3_mid;
                    double curlb_minus_q_y = -gradx_magnetic(ixy) / rho(ixy) - gradx_pressure(ixy) / rho(ixy) / B3_mid;

                    mean_velocity_x_mid(ixy) = mean_velocity_x(ixy) + second * curlb_minus_q_y + third * curlb_minus_q_x;
                    mean_velocity_y_mid(ixy) = mean_velocity_y(ixy) - second * curlb_minus_q_x + third * curlb_minus_q_y;

                    //std::cout << "mean_velocity_x_mid 1 is: " << std::setprecision(16) << mean_velocity_x_mid(ixy) - mean1 << std::endl;
                    //std::cout << "mean_velocity_y_mid 2 is: " << std::setprecision(16) << mean_velocity_y_mid(ixy) - mean2 << std::endl;
                });
            */
            
            
            
            
            
                
            
            
            // next we assemble the two terms on the right side of scheme of Bz
            /// the first term rhs_1
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const ixy) {
                    rhs_1(ixy) = -gradx_pressure(ixy)/rho(ixy) - mean_velocity_y_mid(ixy) * magnetic_field_z_mid(ixy) - magnetic_field_z_mid(ixy) * gradx_magnetic(ixy) / rho(ixy);
                });
            get_gradient<Y>(intermediate_chunk, rhs_1);

            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const ixy) {
                    rhs_2(ixy) = -grady_pressure(ixy)/rho(ixy) + mean_velocity_x_mid(ixy) * magnetic_field_z_mid(ixy) - magnetic_field_z_mid(ixy) * grady_magnetic(ixy) / rho(ixy);
                });
            get_gradient<X>(intermediate_chunk, rhs_2);

            /// update the magnetic field in the picard iteration 
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const ixy) {
                    magnetic_field_z(ixy) = magnetic_field_z_old(ixy) - dt * (-rhs_1(ixy) + rhs_2(ixy) );
                });
            
            
            // next we assemble the four terms on the right side of scheme of p
             /// the first term rhs_1
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        rhs_1(ixy) = (mean_velocity_x_mid(ixy) - grady_magnetic(ixy) / rho(ixy)) * pressure_field_mid(ixy);
                    });

            get_gradient<X>(intermediate_chunk, rhs_1);

             /// the second term rhs_2
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        rhs_2(ixy) = (mean_velocity_y_mid(ixy) + gradx_magnetic(ixy) / rho(ixy)) * pressure_field_mid(ixy);
                    });
            
            get_gradient<Y>(intermediate_chunk, rhs_2);

             /// the thrid term rhs_3
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        rhs_3(ixy) = mean_velocity_x_mid(ixy) - grady_magnetic(ixy) / rho(ixy);
                    });

            get_gradient<X>(intermediate_chunk, rhs_3);
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        rhs_3(ixy) = rhs_3(ixy) * pressure_field_mid(ixy) * (5.0/3.0 - 1);
                    });

             /// the fourth term rhs_4
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        rhs_4(ixy) = mean_velocity_y_mid(ixy) + gradx_magnetic(ixy) / rho(ixy);
                    });
            
            get_gradient<Y>(intermediate_chunk, rhs_4);
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        rhs_4(ixy) = rhs_4(ixy) * pressure_field_mid(ixy) * (5.0/3.0 - 1);
                    });

            /// update the pressure field in the picard iteration 
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        pressure_field(ixy) = pressure_field_old(ixy) - dt * (rhs_1(ixy) + rhs_2(ixy) + rhs_3(ixy) + rhs_4(ixy) );
                    });
            

            // compute the difference bewtween the new p and the p from last iteration
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxXY const ixy) {
                    magnetic_field_z_previous(ixy)
                            = magnetic_field_z_previous(ixy) - magnetic_field_z(ixy);
                });

            // compute the difference bewtween the new p and the p from last iteration
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxXY const ixy) {
                        pressure_field_previous(ixy)
                                = pressure_field_previous(ixy) - pressure_field(ixy);
                    });

            double max_val = norm_inf(Kokkos::DefaultExecutionSpace(), get_const_field(pressure_field_previous))
                           + norm_inf(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_z_previous));
            std::cout << "The iteration error of sub step pvb is: " << max_val << std::endl;

            double tol = 1e-15;
            if (max_val < tol) {
                std::cout << "The picard iteration of sub step pvb is finished." << std::endl;
                break;
            }
        }
        Kokkos::Profiling::popRegion();
        // prepare the advection in vx
        ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                rhs_1(ixy) = grady_pressure(ixy)/magnetic_field_z_mid(ixy)/rho(ixy) - mean_velocity_x_mid(ixy) + grady_magnetic(ixy) / rho(ixy) ;
            });
        
        // prepare the advection in vy
        ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxXY const ixy) {
                rhs_2(ixy) = - gradx_pressure(ixy)/magnetic_field_z_mid(ixy)/rho(ixy) - mean_velocity_y_mid(ixy) - gradx_magnetic(ixy) / rho(ixy) ;
            });
        
            //a = (grid_v2[j] - mean_u2[k] - pressure_mid_gradient[k]/density_array[k]/B3_mid[k]  -  B3_mid_gradient[k]/density_array[k]) * np.tan(-0.5*ddt*0.5*(B3_old[k] + B3[k]))
            //a = (grid_v1[j] - mean_u1[k]) * np.sin(ddt*0.5*(B3_old[k] + B3[k]))

        return pressure_field;
    }
};
