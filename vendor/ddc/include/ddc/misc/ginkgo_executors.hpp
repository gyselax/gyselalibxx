#pragma once

#include <memory>
#include <type_traits>

#include <ginkgo/ginkgo.hpp>

inline std::shared_ptr<gko::Executor> create_default_host_executor()
{
#ifdef KOKKOS_ENABLE_SERIAL
    if constexpr (std::is_same_v<Kokkos::DefaultHostExecutionSpace, Kokkos::Serial>) {
        return gko::ReferenceExecutor::create();
    }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    if constexpr (std::is_same_v<Kokkos::DefaultHostExecutionSpace, Kokkos::OpenMP>) {
        return gko::OmpExecutor::create();
    }
#endif
} // Comes from "Basic Kokkos Extension" Ginkgo MR

template <typename ExecSpace>
inline std::shared_ptr<gko::Executor> create_gko_exec()
{
#ifdef KOKKOS_ENABLE_SERIAL
    if (std::is_same_v<ExecSpace, Kokkos::Serial>) {
        return gko::ReferenceExecutor::create();
    }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    if (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
        return gko::OmpExecutor::create();
    }
#endif
#ifdef KOKKOS_ENABLE_CUDA
    if (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
        return gko::CudaExecutor::create(0, create_default_host_executor());
    }
#endif
#ifdef KOKKOS_ENABLE_HIP
    if (std::is_same_v<ExecSpace, Kokkos::HIP>) {
        return gko::HipExecutor::create(0, create_default_host_executor());
    }
#endif
}
