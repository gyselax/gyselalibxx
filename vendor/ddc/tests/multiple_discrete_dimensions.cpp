// Copyright (C) The DDC development team, see COPYRIGHT.md file
//
// SPDX-License-Identifier: MIT

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include <Kokkos_Core.hpp>

namespace DDC_HIP_5_7_ANONYMOUS_NAMESPACE_WORKAROUND(MULTIPLE_DISCRETE_DIMENSIONS_CPP)
{
    class SingleValueDiscreteDimension
    {
    public:
        using discrete_dimension_type = SingleValueDiscreteDimension;

    public:
        template <class DDim, class MemorySpace>
        class Impl
        {
            template <class ODDim, class OMemorySpace>
            friend class Impl;

        private:
            int m_value;

        public:
            using discrete_dimension_type = SingleValueDiscreteDimension;

            Impl() = default;

            Impl(Impl const&) = delete;

            template <class OriginMemorySpace>
            explicit Impl(Impl<DDim, OriginMemorySpace> const& impl) : m_value(impl.m_value)
            {
            }

            Impl(Impl&&) = default;

            explicit Impl(int value) : m_value(value) {}

            ~Impl() = default;

            KOKKOS_FUNCTION int value() const noexcept
            {
                return m_value;
            }
        };
    };

    struct SVDD1 : SingleValueDiscreteDimension
    {
    };

    struct SVDD2 : SingleValueDiscreteDimension
    {
    };

} // namespace )

TEST(MultipleDiscreteDimensions, Value)
{
    ddc::init_discrete_space<SVDD1>(239);
    ddc::init_discrete_space<SVDD2>(928);
    EXPECT_EQ(ddc::discrete_space<SVDD1>().value(), 239);
    EXPECT_EQ(ddc::discrete_space<SVDD2>().value(), 928);
}
