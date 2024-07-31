// SPDX-License-Identifier: MIT
#pragma once
#include <mpi.h>

#include "ddc_aliases.hpp"

/**
 * @brief A superclass describing an operator for converting from/to different MPI layouts.
 *
 * @tparam Layout1 One of the MPI layouts.
 * @tparam Layout2 The other MPI layouts.
 */
template <class Layout1, class Layout2>
class IMPITranspose
{
    static_assert(ddc::type_seq_same_v<
                  ddc::to_type_seq_t<typename Layout1::discrete_domain_type>,
                  ddc::to_type_seq_t<typename Layout2::discrete_domain_type>>);

public:
    /// The type of the index range of the first MPI layout.
    using idx_range_type1 = typename Layout1::discrete_domain_type;
    /// The type of the index range of the second MPI layout.
    using idx_range_type2 = typename Layout2::discrete_domain_type;

protected:
    /// The MPI communicator
    MPI_Comm m_comm;

public:
    /**
     * @brief A constructor for the class.
     *
     * @param comm The MPI communicator
     */
    explicit IMPITranspose(MPI_Comm comm) : m_comm(comm) {}
};
