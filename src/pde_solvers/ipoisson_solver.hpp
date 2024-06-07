#pragma once

#include <ddc/ddc.hpp>

#include <vector_field_span.hpp>

template <class LaplacianDomain, class FullDomain, class LayoutSpace, class MemorySpace>
class IPoissonSolver;

/**
 * An abstract class from which a Poisson solver can inherit.
 * Classes inheriting from this must implement a way to solve the following equation:
 * @f$ -\Delta \phi = \rho @f$
 *
 * @tparam LaplacianDomain The domain on which the equation is defined.
 * @tparam FullDomain The domain on which the operator() acts. This is equal to the
 *                      LaplacianDomain plus any batched dimensions.
 * @tparam LayoutSpace The layout space of the ChunkSpans passed to operator().
 * @tparam MemorySpace The space (CPU/GPU) where the ChunkSpans passed to operator()
 *                      are saved.
 */
template <class... ODims, class FullDomain, class LayoutSpace, class MemorySpace>
class IPoissonSolver<ddc::DiscreteDomain<ODims...>, FullDomain, LayoutSpace, MemorySpace>
{
protected:
    /// @brief The tags describing the real dimensions in the equation.
    using real_laplacian_tags = ddc::detail::TypeSeq<typename ODims::continuous_dimension_type...>;
    /// @brief The tags describing the discrete dimensions in the equation.
    using laplacian_tags = ddc::detail::TypeSeq<ODims...>;
    /// @brief The tags describing the dimensions of the domain on which the operator acts.
    using space_tags = ddc::to_type_seq_t<FullDomain>;
    /// @brief The tags describing the batched dimensions.
    using batch_tags = ddc::type_seq_remove_t<space_tags, laplacian_tags>;

protected:
    /// @brief Indicates whether the gradient is represented by a VectorFieldSpan or a ChunkSpan.
    static constexpr bool using_vector_span = ddcHelper::type_seq_length_v<laplacian_tags> == 1;

public:
    /// @brief The ChunkSpan type of the arguments to operator().
    using chunk_span_type = ddc::ChunkSpan<double, FullDomain, LayoutSpace, MemorySpace>;
    /// @brief The const ChunkSpan type of the arguments to operator().
    using chunk_view_type = ddc::ChunkView<double, FullDomain, LayoutSpace, MemorySpace>;

    /// @brief The type of the derivative of @f$ \phi @f$.
    using vector_span_type = std::conditional_t<
            ddcHelper::type_seq_length_v<laplacian_tags> == 1,
            chunk_span_type,
            VectorFieldSpan<double, FullDomain, real_laplacian_tags, LayoutSpace, MemorySpace>>;

    /// @brief The DiscreteDomain describing the batch dimensions.
    using batch_domain_type = typename ddc::detail::convert_type_seq_to_discrete_domain<batch_tags>;
    /// @brief The DiscreteElement for indexing a batch dimension.
    using batch_element_type = typename batch_domain_type::discrete_element_type;

    /// @brief The type of the domain on which the equation is defined.
    using laplacian_domain_type = ddc::DiscreteDomain<ODims...>;

    /// @brief The layout space of the ChunkSpans passed to operator().
    using layout_space = LayoutSpace;
    /// @brief The space (CPU/GPU) where the ChunkSpans passed to operator() are saved.
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
    virtual chunk_span_type operator()(chunk_span_type phi, chunk_span_type rho) const = 0;

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
            const = 0;
};
