// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include "ddc_alias_inline_functions.hpp"
#include "ddc_aliases.hpp"
#include "ddc_helper.hpp"
#include "ihybrid_solver_1d.hpp"
#include "vector_index_tools.hpp"
#include "geometry.hpp"
#include "l_norm_tools.hpp"
#include <iostream>
#include <cmath>
#include "species_info.hpp"

/**
 * See @ref 1DFFTHybridSolverImplementation.
 */
template <
        class IdxRangeHybrid,
        class IdxRangeFull = IdxRangeHybrid,
        class ExecSpace = Kokkos::DefaultExecutionSpace,
        class LayoutSpace = Kokkos::layout_right>
class FFTHybridSolver1D;

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
class FFTHybridSolver1D<IdxRange<GridPDEDim1D...>, IdxRangeFull, ExecSpace, LayoutSpace>
    : public IHybridSolver1d<
              IdxRange<GridPDEDim1D...>,
              IdxRangeFull,
              typename ExecSpace::memory_space,
              LayoutSpace>
{
private:
    using base_type = IHybridSolver1d<
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
    explicit FFTHybridSolver1D(hybrid_idx_range_type hybrid_idx_range)
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
        Kokkos::Profiling::pushRegion("FFTHybridSolver1d");

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
     * @brief An operator which calculates the solution @f$ magnetic_field_y, magnetic_field_z@f$ to sub-step bb 
     * in the hybrid model
     * @param[out] magnetic_field_y The By, the result of the nonlinear solver.
     * @param[in] magnetic_field_y_old The By in last time step, the input of the nonlinear solver.
     * @param[in] magnetic_field_y_mid The By in middle time step, the intermidiate value of the nonlinear solver.
     * @param[in] magnetic_field_y_previous The By in previous iteration, the intermidiate value of the nonlinear solver.
     * @param[out] magnetic_field_z The Bz field, the result of the nonlinear solver.
     * @param[in] magnetic_field_z_old The Bz field in last time step, the input of the nonlinear solver.
     * @param[in] magnetic_field_z_mid The Bz field in middle time step, the intermidiate value of the nonlinear solver.
     * @param[in] magnetic_field_z_previous The Bz value of the nonlinear solver, the intermidiate value of the nonlinear solver.
     * @param[in] rho The charge density of all ions.
     * @param[in] magnetic_field_x The Bx value in the last time step.
     * @param[in] gradx_magnetic_field_y_mid The partial_x By, intermidate value.
     * @param[in] gradx_magnetic_field_z_mid The partial_x Bz, intermidate value.
     * @param[in] dt The time step.
     * 
     * @return A reference to the solution of sub-step bb.
    */

    // used for the sub-step bb
    virtual field_type operator()(field_type magnetic_field_y, field_type magnetic_field_y_old,
                                  field_type magnetic_field_y_mid, field_type magnetic_field_y_previous,
                                  field_type magnetic_field_z, field_type magnetic_field_z_old, 
                                  field_type magnetic_field_z_mid, field_type magnetic_field_z_previous,
                                  field_type rho,
                                  field_type magnetic_field_x,
                                  field_type gradx_magnetic_field_y_mid,
                                  field_type gradx_magnetic_field_z_mid,
                                  double const dt) const final
    {
        Kokkos::Profiling::pushRegion("FFTHybridSolver_SUBSTEP_bb1d3v");

        hybrid_idx_range_type idx_range(get_idx_range(magnetic_field_y));
        batch_idx_range_type batch_idx_range(get_idx_range(magnetic_field_y));

        // Build a mesh in the fourier space, for N points
        fourier_idx_range_type const k_mesh = ddc::fourier_mesh<
                GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>(idx_range, false);

        fourier_field_mem_type intermediate_chunk_alloc(k_mesh);

        fourier_field_type intermediate_chunk = get_field(intermediate_chunk_alloc);
        
        // Before the picard iteration, set field_old = field
        ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxX const ixy) {
                magnetic_field_y_old(ixy) = magnetic_field_y(ixy);
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
                    KOKKOS_LAMBDA(IdxX const ixy) {
                        magnetic_field_y_previous(ixy) = magnetic_field_y(ixy);
                        magnetic_field_z_previous(ixy) = magnetic_field_z(ixy);
                        magnetic_field_y_mid(ixy) = 0.5 * (magnetic_field_y_old(ixy) + magnetic_field_y(ixy));
                        magnetic_field_z_mid(ixy) = 0.5 * (magnetic_field_z_old(ixy) + magnetic_field_z(ixy));
                        //std::cout << "rho is: " << std::setprecision(16) << mean_velocity_y(ixy) << std::endl;
                    });

            // compute fft for the magnetic_field_y_mid and magnetic_field_z_mid
            get_gradient<X>(intermediate_chunk, magnetic_field_y_mid, gradx_magnetic_field_y_mid);
            get_gradient<X>(intermediate_chunk, magnetic_field_z_mid, gradx_magnetic_field_z_mid);

            
            // compute the Bx / rho * gradx Bz or gradx By on the right hand side
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxX const ixy) {
                    magnetic_field_y_mid(ixy) = magnetic_field_x(ixy) / rho(ixy) * gradx_magnetic_field_z_mid(ixy);
                    magnetic_field_z_mid(ixy) = magnetic_field_x(ixy) / rho(ixy) * gradx_magnetic_field_y_mid(ixy);;
                });

            // compute fft for the Bx / rho * gradx Bz or gradx By on the right hand side 
            get_gradient<X>(intermediate_chunk, magnetic_field_y_mid, gradx_magnetic_field_y_mid);
            get_gradient<X>(intermediate_chunk, magnetic_field_z_mid, gradx_magnetic_field_z_mid);

            // update By and Bz 
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxX const ixy) {
                    magnetic_field_y(ixy) = magnetic_field_y_old(ixy) + dt * gradx_magnetic_field_y_mid(ixy);
                    magnetic_field_z(ixy) = magnetic_field_z_old(ixy) - dt * gradx_magnetic_field_z_mid(ixy);
                });

            
            // compute the difference bewtween the new By, Bz and the By, Bz from last iteration
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxX const ixy) {
                        magnetic_field_y_previous(ixy)
                                = magnetic_field_y(ixy) - magnetic_field_y_previous(ixy);
                        magnetic_field_z_previous(ixy)
                                = magnetic_field_z(ixy) - magnetic_field_z_previous(ixy);
                    });

            double max_val = norm_inf(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_y_previous))
                           + norm_inf(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_z_previous));
            std::cout << "The iteration error of sub step bb is: " << max_val << std::endl;

            double tol = 1e-15;
            if (max_val < tol) {
                std::cout << "The picard iteration of sub step bb is finished." << std::endl;
                break;
            }

        }
        Kokkos::Profiling::popRegion();
        
        return magnetic_field_y;
    }



    /**  Isothermal electron case.
     * @brief An operator which calculates the solution @f$ magnetic_field_y, magnetic_field_z@f$ to sub-step fb.
     * in the hybrid model
     * @param[out] magnetic_field_y The By, the result of the nonlinear solver.
     * @param[in] magnetic_field_y_old The By in last time step, the input of the nonlinear solver.
     * @param[in] magnetic_field_y_mid The By in middle time step, the intermidiate value of the nonlinear solver.
     * @param[in] magnetic_field_y_previous The By in previous iteration, the intermidiate value of the nonlinear solver.
     * @param[out] magnetic_field_z The Bz field, the result of the nonlinear solver.
     * @param[in] magnetic_field_z_old The Bz field in last time step, the input of the nonlinear solver.
     * @param[in] magnetic_field_z_mid The Bz field in middle time step, the intermidiate value of the nonlinear solver.
     * @param[in] magnetic_field_z_previous The Bz value of the nonlinear solver, the intermidiate value of the nonlinear solver.
     * @param[in] rho The charge density of all ions.
     * @param[in] magnetic_field_x The Bx value in the last time step.
     * @param[in] gradx_magnetic_field_y_mid The partial_x By, intermidate value.
     * @param[in] gradx_magnetic_field_z_mid The partial_x Bz, intermidate value.
     * @param[in] dt The time step.
     * 
     * @return A reference to the solution of sub-step fb.
    */

    // used for the sub-step fb
    virtual field_type operator()(field_type magnetic_field_y, field_type magnetic_field_y_old,
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
                                  double const dt) const final
    {
        Kokkos::Profiling::pushRegion("FFTHybridSolver_SUBSTEP_fb1d3v");

        hybrid_idx_range_type idx_range(get_idx_range(magnetic_field_y));
        batch_idx_range_type batch_idx_range(get_idx_range(magnetic_field_y));

        // Build a mesh in the fourier space, for N points
        fourier_idx_range_type const k_mesh = ddc::fourier_mesh<
                GridFourier<typename GridPDEDim1D::continuous_dimension_type>...>(idx_range, false);

        fourier_field_mem_type intermediate_chunk_alloc(k_mesh);

        fourier_field_type intermediate_chunk = get_field(intermediate_chunk_alloc);

        // the charge and mass of each kinds on ions
        
        IdxRangeSp const kin_species_idx_range = get_idx_range<Species>(u_old_x);
        
        host_t<DConstFieldSp> const charges_host = ddc::host_discrete_space<Species>().charges();
        host_t<DConstFieldSp> const kinetic_charges_host
                = charges_host[get_idx_range<Species>(u_old_x)];

        auto const kinetic_charges_alloc
                = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_charges_host);
        DConstFieldSp kinetic_charges = get_const_field(kinetic_charges_alloc);

        host_t<DConstFieldSp> const masses_host = ddc::host_discrete_space<Species>().masses();
        host_t<DConstFieldSp> const kinetic_masses_host
                = masses_host[get_idx_range<Species>(u_old_x)];

        auto const kinetic_masses_alloc
                = create_mirror_view_and_copy(Kokkos::DefaultExecutionSpace(), kinetic_masses_host);
        DConstFieldSp kinetic_masses = get_const_field(kinetic_masses_alloc);
        
        // Before the picard iteration, set field_old = field
        ddc::parallel_for_each(
            Kokkos::DefaultExecutionSpace(),
            idx_range,
            KOKKOS_LAMBDA(IdxX const ixy) {
                magnetic_field_y_old(ixy) = magnetic_field_y(ixy);
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
                    KOKKOS_LAMBDA(IdxX const ixy) {
                        magnetic_field_y_previous(ixy) = magnetic_field_y(ixy);
                        magnetic_field_z_previous(ixy) = magnetic_field_z(ixy);
                        magnetic_field_y_mid(ixy) = 0.5 * (magnetic_field_y_old(ixy) + magnetic_field_y(ixy));
                        magnetic_field_z_mid(ixy) = 0.5 * (magnetic_field_z_old(ixy) + magnetic_field_z(ixy));

                        Mxx(ixy) = 0.0;
                        Mxy(ixy) = 0.0;
                        Mxz(ixy) = 0.0;

                        Myx(ixy) = 0.0;
                        Myy(ixy) = 0.0;
                        Myz(ixy) = 0.0;

                        Mzx(ixy) = 0.0;
                        Mzy(ixy) = 0.0;
                        Mzz(ixy) = 0.0;

                        weighted_u_x(ixy) = 0.0;
                        weighted_u_y(ixy) = 0.0;
                        weighted_u_z(ixy) = 0.0;

                        weighted_p_para_x(ixy) = 0.0;
                        weighted_p_para_y(ixy) = 0.0;
                        weighted_p_para_z(ixy) = 0.0;


                        //std::cout << "rho is: " << std::setprecision(16) << mean_velocity_y(ixy) << std::endl;
                    });


            // compute fft for the magnetic_field_y_mid, magnetic_field_z_mid, and rho.
            get_gradient<X>(intermediate_chunk, magnetic_field_y_mid, gradx_magnetic_field_y_mid);
            get_gradient<X>(intermediate_chunk, magnetic_field_z_mid, gradx_magnetic_field_z_mid);
            get_gradient<X>(intermediate_chunk, rho, gradx_rho);


            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxX const ixy) {
                        double B_square = magnetic_field_x(ixy) * magnetic_field_x(ixy) + magnetic_field_y_mid(ixy) * magnetic_field_y_mid(ixy) + magnetic_field_z_mid(ixy) * magnetic_field_z_mid(ixy);
                        qx(ixy) = 0.0;
                        qy(ixy) = magnetic_field_z_mid(ixy) * gradx_rho(ixy) * electron_temperature / B_square / rho(ixy);
                        qz(ixy) = - magnetic_field_y_mid(ixy) * gradx_rho(ixy) * electron_temperature / B_square / rho(ixy);
                        p_parallel_x(ixy) = gradx_rho(ixy) * magnetic_field_x(ixy) * electron_temperature / rho(ixy) * magnetic_field_x(ixy) / B_square;
                        p_parallel_y(ixy) = gradx_rho(ixy) * magnetic_field_x(ixy) * electron_temperature / rho(ixy) * magnetic_field_y_mid(ixy) / B_square;
                        p_parallel_z(ixy) = gradx_rho(ixy) * magnetic_field_x(ixy) * electron_temperature / rho(ixy) * magnetic_field_z_mid(ixy) / B_square;
 
                    });

            // decomposition of p/rho, and compute u_bar.
            ddc::for_each(kin_species_idx_range, [&](IdxSp isp) {
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxX const ixy) {
                    double B_square = magnetic_field_x(ixy) * magnetic_field_x(ixy) + magnetic_field_y_mid(ixy) * magnetic_field_y_mid(ixy) + magnetic_field_z_mid(ixy) * magnetic_field_z_mid(ixy);
                    double q_over_m = - kinetic_charges(isp) / kinetic_masses(isp);
                    double rho_i_over_rho = rho_each(isp,ixy) / rho(ixy);
                    // matrix M 
                    double coe1 = 2.0 * std::sin(0.5*dt*q_over_m*std::sqrt(B_square)) * std::sin(0.5*dt*q_over_m*std::sqrt(B_square)) / dt / B_square/q_over_m; 
                    double coe2 = (1.0 - std::sin( dt * q_over_m * std::sqrt(B_square) ) / ( dt*q_over_m*std::sqrt(B_square) ) ) / B_square;

                    double Mxx_temp = 1.0 - coe2 * (magnetic_field_y_mid(ixy) * magnetic_field_y_mid(ixy) + magnetic_field_z_mid(ixy) * magnetic_field_z_mid(ixy));
                    double Mxy_temp = coe1 * magnetic_field_z_mid(ixy) + coe2 * magnetic_field_x(ixy) * magnetic_field_y_mid(ixy);
                    double Mxz_temp = - coe1 * magnetic_field_y_mid(ixy) + coe2 * magnetic_field_x(ixy) * magnetic_field_z_mid(ixy);

                    double Myx_temp = - coe1 *  magnetic_field_z_mid(ixy) + coe2 * magnetic_field_x(ixy) * magnetic_field_y_mid(ixy);
                    double Myy_temp = 1.0 - coe2 * (magnetic_field_x(ixy) * magnetic_field_x(ixy) + magnetic_field_z_mid(ixy) * magnetic_field_z_mid(ixy));
                    double Myz_temp = coe1 * magnetic_field_x(ixy) + coe2 * magnetic_field_y_mid(ixy) * magnetic_field_z_mid(ixy);

                    double Mzx_temp = coe1 * magnetic_field_y_mid(ixy) + coe2 * magnetic_field_x(ixy) * magnetic_field_z_mid(ixy);
                    double Mzy_temp = -coe1 * magnetic_field_x(ixy) + coe2 * magnetic_field_y_mid(ixy) *  magnetic_field_z_mid(ixy);
                    double Mzz_temp = 1.0 - coe2 * (magnetic_field_x(ixy)*magnetic_field_x(ixy) + magnetic_field_y_mid(ixy) * magnetic_field_y_mid(ixy));

                    Mxx(ixy) += rho_i_over_rho * Mxx_temp;
                    Mxy(ixy) += rho_i_over_rho * Mxy_temp;
                    Mxz(ixy) += rho_i_over_rho * Mxz_temp;

                    Myx(ixy) += rho_i_over_rho * Myx_temp;
                    Myy(ixy) += rho_i_over_rho * Myy_temp;
                    Myz(ixy) += rho_i_over_rho * Myz_temp;

                    Mzx(ixy) += rho_i_over_rho * Mzx_temp;
                    Mzy(ixy) += rho_i_over_rho * Mzy_temp;
                    Mzz(ixy) += rho_i_over_rho * Mzz_temp;

                    weighted_u_x(ixy) += rho_i_over_rho * (Mxx_temp * u_old_x(isp,ixy) + Mxy_temp * u_old_y(isp,ixy) + Mxz_temp * u_old_z(isp,ixy));
                    weighted_u_y(ixy) += rho_i_over_rho * (Myx_temp * u_old_x(isp,ixy) + Myy_temp * u_old_y(isp,ixy) + Myz_temp * u_old_z(isp,ixy));
                    weighted_u_z(ixy) += rho_i_over_rho * (Mzx_temp * u_old_x(isp,ixy) + Mzy_temp * u_old_y(isp,ixy) + Mzz_temp * u_old_z(isp,ixy));

                    weighted_p_para_x(ixy) += 0.5 * dt * q_over_m * rho_i_over_rho * p_parallel_x(ixy);
                    weighted_p_para_y(ixy) += 0.5 * dt * q_over_m * rho_i_over_rho * p_parallel_y(ixy);
                    weighted_p_para_z(ixy) += 0.5 * dt * q_over_m * rho_i_over_rho * p_parallel_z(ixy);
                });
            });


            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxX const ixy) {
                    // the inverse of M, denote by N
                    double det_M = Mxx(ixy) * (Myy(ixy) * Mzz(ixy) - Myz(ixy) * Mzy(ixy)) - Mxy(ixy) * (Myx(ixy) * Mzz(ixy) - Myz(ixy) * Mzx(ixy)) + Mxz(ixy) * (Myx(ixy) * Mzy(ixy) - Myy(ixy) * Mzx(ixy));
                    
                    double Nxx = (Myy(ixy) * Mzz(ixy) - Myz(ixy) * Mzy(ixy)) / det_M;
                    double Nxy = (Mxz(ixy) * Mzy(ixy) - Mxy(ixy) * Mzz(ixy)) / det_M;
                    double Nxz = (Mxy(ixy) * Myz(ixy) - Mxz(ixy) * Myy(ixy)) / det_M;

                    double Nyx = (Myz(ixy) * Mzx(ixy) - Myx(ixy) * Mzz(ixy)) / det_M;
                    double Nyy = (Mxx(ixy) * Mzz(ixy) - Mxz(ixy) * Mzx(ixy)) / det_M;
                    double Nyz = (Mxz(ixy) * Myx(ixy) - Mxx(ixy) * Myz(ixy)) / det_M;

                    double Nzx = (Myx(ixy) * Mzy(ixy) - Myy(ixy) * Mzx(ixy)) / det_M;
                    double Nzy = (Mxy(ixy) * Mzx(ixy) - Mxx(ixy) * Mzy(ixy)) / det_M;
                    double Nzz = (Mxx(ixy) * Myy(ixy) - Mxy(ixy) * Myx(ixy)) / det_M;

                    // hall over rho minus q
                    double horq_x = -qx(ixy);
                    double horq_y = -gradx_magnetic_field_z_mid(ixy)/rho(ixy) - qy(ixy);
                    double horq_z = gradx_magnetic_field_y_mid(ixy)/rho(ixy) - qz(ixy);

                    u_bar_x(ixy) = horq_x + Nxx * (weighted_u_x(ixy) - horq_x - weighted_p_para_x(ixy)) + Nxy * (weighted_u_y(ixy) - horq_y - weighted_p_para_y(ixy)) + Nxz * (weighted_u_z(ixy) - horq_z - weighted_p_para_z(ixy));
                    u_bar_y(ixy) = horq_y + Nyx * (weighted_u_x(ixy) - horq_x - weighted_p_para_x(ixy)) + Nyy * (weighted_u_y(ixy) - horq_y - weighted_p_para_y(ixy)) + Nyz * (weighted_u_z(ixy) - horq_z - weighted_p_para_z(ixy));
                    u_bar_z(ixy) = horq_z + Nzx * (weighted_u_x(ixy) - horq_x - weighted_p_para_x(ixy)) + Nzy * (weighted_u_y(ixy) - horq_y - weighted_p_para_y(ixy)) + Nzz * (weighted_u_z(ixy) - horq_z - weighted_p_para_z(ixy));

                });


            
            // compute the Bx / rho * gradx Bz or gradx By on the right hand side
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxX const ixy) {
                    rhs_2(ixy) = u_bar_z(ixy) * magnetic_field_x(ixy) - u_bar_x(ixy) * magnetic_field_z_mid(ixy) - magnetic_field_x(ixy) * gradx_magnetic_field_y_mid(ixy) / rho(ixy);
                    rhs_3(ixy) = u_bar_x(ixy) * magnetic_field_y_mid(ixy) - u_bar_y(ixy) * magnetic_field_x(ixy) - magnetic_field_x(ixy) * gradx_magnetic_field_z_mid(ixy) / rho(ixy);
                });

            // compute fft for the Bx / rho * gradx Bz or gradx By on the right hand side 
            get_gradient<X>(intermediate_chunk, rhs_2, rhs_5);
            get_gradient<X>(intermediate_chunk, rhs_3, rhs_6);

            // update By and Bz 
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxX const ixy) {
                    magnetic_field_y(ixy) = magnetic_field_y_old(ixy) - dt * rhs_6(ixy);
                    magnetic_field_z(ixy) = magnetic_field_z_old(ixy) + dt * rhs_5(ixy);
                });

            
            // compute the difference bewtween the new By, Bz and the By, Bz from last iteration
            ddc::parallel_for_each(
                    Kokkos::DefaultExecutionSpace(),
                    idx_range,
                    KOKKOS_LAMBDA(IdxX const ixy) {
                        magnetic_field_y_previous(ixy)
                                = magnetic_field_y(ixy) - magnetic_field_y_previous(ixy);
                        magnetic_field_z_previous(ixy)
                                = magnetic_field_z(ixy) - magnetic_field_z_previous(ixy);
                    });

            double max_val = norm_inf(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_y_previous))
                           + norm_inf(Kokkos::DefaultExecutionSpace(), get_const_field(magnetic_field_z_previous));
            std::cout << "The iteration error of sub step bb is: " << max_val << std::endl;

            double tol = 1e-15;
            if (max_val < tol) {
                std::cout << "The picard iteration of sub step bb is finished." << std::endl;
                break;
            }

        }
        Kokkos::Profiling::popRegion();

        // prepare for the shifted rotation
            ddc::parallel_for_each(
                Kokkos::DefaultExecutionSpace(),
                idx_range,
                KOKKOS_LAMBDA(IdxX const ixy) {
                    rhs_1(ixy) = u_bar_x(ixy) + qx(ixy);
                    rhs_2(ixy) = u_bar_y(ixy) + qy(ixy) + gradx_magnetic_field_z_mid(ixy) / rho(ixy);
                    rhs_3(ixy) = u_bar_z(ixy) + qz(ixy) - gradx_magnetic_field_y_mid(ixy) / rho(ixy);

                    //std::cout << "yl_x: " << std::setprecision(16) << rhs_1(ixy) << std::endl;
                    //std::cout << "yl_y: " << std::setprecision(16) << rhs_2(ixy) << std::endl;
                    //std::cout << "yl_z: " << std::setprecision(16) << rhs_3(ixy) << std::endl;
                });
        
        return magnetic_field_y;


    }




};
