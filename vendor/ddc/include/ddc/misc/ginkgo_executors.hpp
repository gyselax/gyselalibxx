#pragma once

#include <memory>
#include <type_traits>

#include <ginkgo/ginkgo.hpp>

#include <Kokkos_Core.hpp>

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
    if constexpr (std::is_same_v<ExecSpace, Kokkos::Serial>) {
        return gko::ReferenceExecutor::create();
    }
#endif
#ifdef KOKKOS_ENABLE_OPENMP
    if constexpr (std::is_same_v<ExecSpace, Kokkos::OpenMP>) {
        return gko::OmpExecutor::create();
    }
#endif
#ifdef KOKKOS_ENABLE_CUDA
    if constexpr (std::is_same_v<ExecSpace, Kokkos::Cuda>) {
        ExecSpace exec_space;
        return gko::CudaExecutor::
                create(exec_space.cuda_device(),
                       create_default_host_executor(),
                       std::make_shared<gko::CudaAllocator>(),
                       exec_space.cuda_stream());
    }
#endif
#ifdef KOKKOS_ENABLE_HIP
    if constexpr (std::is_same_v<ExecSpace, Kokkos::HIP>) {
        ExecSpace exec_space;
        return gko::HipExecutor::
                create(exec_space.hip_device(),
                       create_default_host_executor(),
                       std::make_shared<gko::HipAllocator>(),
                       exec_space.hip_stream());
    }
#endif
}
