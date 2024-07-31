// SPDX-License-Identifier: MIT
#pragma once
#include <ddc/ddc.hpp>

#include <impilayout.hpp>

/**
 * @brief A class describing a way in which data may be laid out across MPI processes.
 *
 * This class describes the simplest way of laying out data. In this layout the data is
 * distributed along dimensions in order until the data is distributed. It is possible that
 * the data may not be distributed along some of the requested dimensions if such distribution
 * is not necessary.
 * For example if we distribute dimensions R and Theta of a (R, Theta, Phi, Vpar, Mu) grid
 * of size 256 x 1024 x 256 x 64 x 8 over 64 processes, the theta dimension will not be
 * distributed.
 *
 * The data is distributed such that it is maximally distributed over each dimension
 * (in order) such that the size of each local domain along that dimension is the same for
 * all processes.
 * For example if we distribute dimensions X and Y of a (X, Y, Z) grid of size (10, 15, 4)
 * over 6 processes, the X dimension will be distributed over 2 processes and the Y
 * dimension will be distributed over 3 processes.
 *
 * @tparam Domain The DiscreteDomain on which the data is defined.
 * @tparam DistributedDim The tags of the discrete dimensions which are distributed
 *              across MPI processes.
 */
template <class Domain, class... DistributedDim>
class MPILayout : public IMPILayout<Domain, DistributedDim...>
{
    static_assert(ddc::is_discrete_domain_v<Domain>);

    using base_type = IMPILayout<Domain, DistributedDim...>;

public:
    /// The domain of the data
    using discrete_domain_type = Domain;
    /// The domain of the distributed section of the data
    using distributed_subdomain = typename base_type::distributed_subdomain;
    /// A type sequence describing the dimensions which are distributed across MPI processes.
    using distributed_type_seq = typename base_type::distributed_type_seq;

public:
    /**
     * @brief Get the distributed domain which follows the chosen layout.
     *
     * @param[in] global_domain The global (non-distributed) domain.
     * @param[in] comm_size The number of MPI processes over which the data is distributed.
     * @param[in] rank The rank of the current MPI process.
     *
     * @returns The distributed domain.
     */
    discrete_domain_type distribute_domain(
            discrete_domain_type global_domain,
            int comm_size,
            int rank)
    {
        return internal_distribute_domain(global_domain, comm_size, rank);
    }



protected:
    /**
     * @brief Distribute a 1D domain over the MPI processes.
     *
     * @param[in] global_domain The domain to be distributed.
     * @param[in] comm_size The number of processes over which the data should be disctributed.
     * @param[in] rank The rank of the process within the current group of processes
     *
     * @returns The distributed domain.
     */
    template <class HeadTag>
    ddc::DiscreteDomain<HeadTag> internal_distribute_domain(
            ddc::DiscreteDomain<HeadTag> global_domain,
            int comm_size,
            int rank)
    {
        if constexpr (ddc::in_tags_v<HeadTag, distributed_type_seq>) {
            assert(global_domain.size() % comm_size == 0);
            ddc::DiscreteVector<HeadTag> elems_on_dim(global_domain.size() / comm_size);
            ddc::DiscreteDomain<HeadTag>
                    local_domain(global_domain.front() + rank * elems_on_dim, elems_on_dim);
            return local_domain;
        } else {
            // Data is not actually distributed as it handles the case of a domain which is not defined on a distributed dimension.
            assert(comm_size == 1);
            return global_domain;
        }
    }

    /**
     * @brief Distribute the domain over the MPI processes.
     *
     * This function is called recursively. At each pass it distributes the domain over the
     * first dimension in the discrete domain. The remaining dimensions and processes are
     * then handled in the recursive call.
     *
     * @param[in] domain The domain to be distributed.
     * @param[in] comm_size The number of processes over which the data should be disctributed.
     * @param[in] rank The rank of the process within the current group of processes
     *
     * @returns The distributed domain.
     */
    template <class HeadTag, class... Tags, std::enable_if_t<(sizeof...(Tags) > 0), bool> = true>
    ddc::DiscreteDomain<HeadTag, Tags...> internal_distribute_domain(
            ddc::DiscreteDomain<HeadTag, Tags...> domain,
            int comm_size,
            int rank)
    {
        ddc::DiscreteDomain<HeadTag> global_domain_along_dim = ddc::select<HeadTag>(domain);
        ddc::DiscreteDomain<HeadTag> local_domain_along_dim;
        ddc::DiscreteDomain<Tags...> remaining_domain;

        if constexpr (ddc::in_tags_v<HeadTag, distributed_type_seq>) {
            // The number of MPI processes along this dimension
            int n_ranks_along_dim = std::gcd(comm_size, global_domain_along_dim.size());
            // The number of MPI processes along all subsequent dimensions
            int n_elems_lower_dims = comm_size / n_ranks_along_dim;
            // The rank index for the MPI process along this dimension
            int rank_along_dim = rank / n_elems_lower_dims;
            // The rank index for the MPI process along all subsequent dimensions
            int remaining_rank = rank % n_elems_lower_dims;
            // Calculate the local domain
            ddc::DiscreteVector<HeadTag> elems_on_dim(
                    global_domain_along_dim.size() / n_ranks_along_dim);
            ddc::DiscreteElement<HeadTag> distrib_start(
                    global_domain_along_dim.front() + rank_along_dim * elems_on_dim);
            local_domain_along_dim = ddc::DiscreteDomain<HeadTag>(distrib_start, elems_on_dim);
            // Calculate the domain for the subsequent dimensions
            ddc::DiscreteDomain<Tags...> remaining_dims = ddc::select<Tags...>(domain);
            remaining_domain = internal_distribute_domain(
                    remaining_dims,
                    n_elems_lower_dims,
                    remaining_rank);
        } else {
            // Calculate the local domain
            local_domain_along_dim = global_domain_along_dim;
            // Calculate the domain for the subsequent dimensions
            remaining_domain
                    = internal_distribute_domain(ddc::select<Tags...>(domain), comm_size, rank);
        }
        return ddc::DiscreteDomain<HeadTag, Tags...>(local_domain_along_dim, remaining_domain);
    }
};
