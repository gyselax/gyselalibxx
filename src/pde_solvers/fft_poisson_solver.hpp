#pragma once

#include <ddc/ddc.hpp>
#include <ddc/kernels/fft.hpp>

#include <ddc_helper.hpp>
#include <directional_tag.hpp>

#include "ipoisson_solver.hpp"

/**
 * See @ref FFTPoissonSolverImplementation.
 */
template <
        class LaplacianDomain,
        class FullDomain,
        class ExecSpace,
        class LayoutSpace = std::experimental::layout_right>
class FFTPoissonSolver;

/**
 * @brief A class to solve the following equation:
 * @f$ -\Delta \phi = \rho @f$
 * using a Fourier transform.
 *
 * The implementation of this class can be found at FFTPoissonSolver< ddc::DiscreteDomain<DDims...>, FullDomain, ExecSpace, LayoutSpace >.
 * @anchor FFTPoissonSolverImplementation
 *
 * @tparam LaplacianDomain The domain on which the equation is defined.
 * @tparam FullDomain The domain on which the operator() acts. This is equal to the
 *                      LaplacianDomain plus any batched dimensions.
 * @tparam ExecSpace The space (CPU/GPU) where the calculations will take place.
 * @tparam LayoutSpace The layout space of the ChunkSpans passed to operator().
 */
template <class... DDims, class FullDomain, class ExecSpace, class LayoutSpace>
class FFTPoissonSolver<ddc::DiscreteDomain<DDims...>, FullDomain, ExecSpace, LayoutSpace>
    : public IPoissonSolver<
              ddc::DiscreteDomain<DDims...>,
              FullDomain,
              LayoutSpace,
              typename ExecSpace::memory_space>
{
private:
    using base_type = IPoissonSolver<
            ddc::DiscreteDomain<DDims...>,
            FullDomain,
            LayoutSpace,
            typename ExecSpace::memory_space>;

public:
    template <class RDim>
    struct IDimFourier : ddc::PeriodicSampling<ddc::Fourier<RDim>>
    {
    };

public:
    /// @brief The ChunkSpan type of the arguments to operator().
    using chunk_span_type = typename base_type::chunk_span_type;
    /// @brief The const ChunkSpan type of the arguments to operator().
    using chunk_view_type = typename base_type::chunk_view_type;

    /// @brief The type of the derivative of @f$ \phi @f$.
    using vector_span_type = typename base_type::vector_span_type;

    /// @brief The DiscreteDomain describing the batch dimensions.
    using batch_domain_type = typename base_type::batch_domain_type;
    /// @brief The DiscreteElement for indexing a batch dimension.
    using batch_element_type = typename base_type::batch_element_type;

    /// @brief The type of the domain on which the equation is defined.
    using laplacian_domain_type = typename base_type::laplacian_domain_type;

    /// @brief The layout space of the ChunkSpans passed to operator().
    using layout_space = typename base_type::layout_space;
    /// @brief The space (CPU/GPU) where the ChunkSpans passed to operator() are saved.
    using memory_space = typename base_type::memory_space;

    /// @brief The type of the Fourier space domain.
    using fourier_domain_type
            = ddc::DiscreteDomain<IDimFourier<typename DDims::continuous_dimension_type>...>;
    /// @brief The type of an index of the Fourier space domain.
    using fourier_element_type = typename fourier_domain_type::discrete_element_type;

    /// @brief The type of a Chunk storing the Fourier transform of a function.
    using fourier_chunk_type = ddc::Chunk<
            Kokkos::complex<double>,
            fourier_domain_type,
            ddc::KokkosAllocator<Kokkos::complex<double>, memory_space>>;
    /// @brief The type of a ChunkSpan storing the Fourier transform of a function.
    using fourier_span_type = typename fourier_chunk_type::span_type;

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
    KOKKOS_FUNCTION static double get_laplace_operator(ddc::DiscreteElement<FDim...> index)
    {
        return (((double)ddc::coordinate(ddc::select<FDim>(index))
                 * (double)ddc::coordinate(ddc::select<FDim>(index)))
                + ...);
    }

    /**
     * @brief Differentiate an expression from its representation in Fourier space by multiplying
     * by -i * k and then converting back to real space.
     *
     * @param[out] derivative The ChunkSpan where the derivative will be saved.
     * @param[out] fourier_derivative The derivative of the function in Fourier space will be saved.
     * @param[in] values The ChunkSpan containing the values of the function in Fourier space.
     *
     * @tparam Dim The dimension along which the expression is differentiated.
     */
    template <class Dim>
    void differentiate_and_invert_fourier_values(
            ddc::ChunkSpan<double, laplacian_domain_type, LayoutSpace, memory_space> derivative,
            fourier_span_type fourier_derivative,
            fourier_span_type values) const
    {
        negative_differentiate_equation<Dim>(fourier_derivative, values);
        // Perform the inverse 1D FFT of the solution to deduce the electric field
        ddc::
                ifft(ExecSpace(),
                     derivative.span_view(),
                     fourier_derivative.span_view(),
                     ddc::kwArgs_fft {m_norm});
    }

    /**
     * @brief Get the gradient of a 1D expression from its representation in Fourier space.
     *
     * @param[out] gradient The ChunkSpan where the derivative will be saved.
     * @param[out] fourier_derivative The derivative of the function in Fourier space will be saved.
     * @param[in] values The ChunkSpan containing the values of the function in Fourier space.
     */
    void get_gradient(
            ddc::ChunkSpan<double, laplacian_domain_type, LayoutSpace, memory_space> gradient,
            fourier_span_type fourier_derivative,
            fourier_span_type values) const
    {
        using RDim =
                typename ddc::type_seq_element_t<0, ddc::to_type_seq_t<laplacian_domain_type>>::
                        continuous_dimension_type;
        using FourierDim = IDimFourier<RDim>;
        differentiate_and_invert_fourier_values<FourierDim>(gradient, fourier_derivative, values);
    }

    /**
     * @brief Get the gradient of a multi-dimensional expression from its representation in Fourier space.
     *
     * @param[out] gradient The VectorFieldSpan where the gradient will be saved.
     * @param[out] fourier_derivative The derivative of the function in Fourier space will be saved.
     * @param[in] values The ChunkSpan containing the values of the function in Fourier space.
     */
    template <class... Dims>
    void get_gradient(
            VectorFieldSpan<
                    double,
                    laplacian_domain_type,
                    NDTag<Dims...>,
                    layout_space,
                    memory_space> gradient,
            fourier_span_type fourier_derivative,
            fourier_span_type values) const
    {
        ((differentiate_and_invert_fourier_values<
                 IDimFourier<Dims>>(ddcHelper::get<Dims>(gradient), fourier_derivative, values)),
         ...);
    }

    template <class IDim>
    void init_fourier_space(ddc::DiscreteDomain<IDim> dom)
    {
        using IDimF = IDimFourier<typename IDim::continuous_dimension_type>;
        ddc::init_discrete_space<IDimF>(ddc::init_fourier_space<IDimF>(dom));
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
            fourier_span_type intermediate_chunk,
            ddc::ChunkSpan<double, laplacian_domain_type, Layout, memory_space> rho) const
    {
        // Compute FFT(rho)
        ddc::fft(ExecSpace(), intermediate_chunk, rho, ddc::kwArgs_fft {m_norm});

        fourier_domain_type const k_mesh = intermediate_chunk.domain();

        // Solve Poisson's equation -\Delta phi = -(\sum_j \partial_j^2) \phi = rho
        //   in Fourier space as -(\sum_j i*k_i * i*k_i) FFT(Phi) = FFT(rho))
        ddc::parallel_for_each(
                ExecSpace(),
                k_mesh,
                KOKKOS_LAMBDA(fourier_element_type const ik) {
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
     * @param derivative The ChunkSpan where the derivative will be saved.
     * @param values The ChunkSpan containing the values of the function in Fourier space.
     *
     * @tparam Dim The dimension along which the expression is differentiated.
     */
    template <class Dim>
    void negative_differentiate_equation(fourier_span_type derivative, fourier_span_type values)
            const
    {
        Kokkos::complex<double> imaginary_unit(0.0, 1.0);
        ddc::parallel_for_each(
                ExecSpace(),
                values.domain(),
                KOKKOS_LAMBDA(fourier_element_type const ik) {
                    ddc::DiscreteElement<Dim> const ikx = ddc::select<Dim>(ik);
                    derivative(ik) = -imaginary_unit * ddc::coordinate(ikx) * values(ik);
                });
    }

public:
    /**
     * @brief A constructor for the FFT Poisson solver.
     * This constructor calls ddc::init_discrete_space so it should only be called once per
     * simulation.
     *
     * @param laplacian_domain The domain on which the equation should be solved.
     */
    FFTPoissonSolver(laplacian_domain_type laplacian_domain)
    {
        ((init_fourier_space<DDims>(ddc::DiscreteDomain<DDims>(laplacian_domain))), ...);
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
    virtual chunk_span_type operator()(chunk_span_type phi, chunk_span_type rho) const final
    {
        Kokkos::Profiling::pushRegion("FFTPoissonSolver");

        laplacian_domain_type domain(phi.domain());
        batch_domain_type batch_domain(phi.domain());

        // Build a mesh in the fourier space, for N points
        fourier_domain_type const k_mesh = ddc::FourierMesh<
                IDimFourier<typename DDims::continuous_dimension_type>...>(domain, false);

        fourier_chunk_type intermediate_chunk_alloc(k_mesh);
        fourier_span_type intermediate_chunk = intermediate_chunk_alloc.span_view();

        ddc::for_each(batch_domain, [&](batch_element_type ib) {
            solve_poisson_equation(intermediate_chunk, rho[ib]);

            // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
            ddc::
                    ifft(ExecSpace(),
                         phi[ib],
                         intermediate_chunk.span_view(),
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
    virtual chunk_span_type operator()(chunk_span_type phi, vector_span_type E, chunk_span_type rho)
            const final
    {
        Kokkos::Profiling::pushRegion("FFTPoissonSolver");

        laplacian_domain_type domain(phi.domain());
        batch_domain_type batch_domain(phi.domain());

        // Build a mesh in the fourier space, for N points
        fourier_domain_type const k_mesh = ddc::FourierMesh<
                IDimFourier<typename DDims::continuous_dimension_type>...>(domain, false);

        fourier_chunk_type intermediate_chunk_alloc(k_mesh);
        fourier_chunk_type fourier_efield_alloc(k_mesh);

        fourier_span_type intermediate_chunk = intermediate_chunk_alloc.span_view();
        fourier_span_type fourier_efield = fourier_efield_alloc.span_view();

        ddc::for_each(batch_domain, [&](batch_element_type ib) {
            solve_poisson_equation(intermediate_chunk, rho[ib]);
            get_gradient(E[ib], fourier_efield, intermediate_chunk);

            // Perform the inverse 1D FFT of the solution to deduce the electrostatic potential
            ddc::
                    ifft(ExecSpace(),
                         phi[ib],
                         intermediate_chunk.span_view(),
                         ddc::kwArgs_fft {m_norm});
        });
        Kokkos::Profiling::popRegion();
        return phi;
    }
};
