// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

using T = double;
using A = ddc::AlignedAllocator<T, 64>;
using U = char;
using B = std::allocator_traits<A>::rebind_alloc<U>;

TEST(AlignedAllocatorTest, Traits)
{
    using traits = std::allocator_traits<A>;
    EXPECT_TRUE((std::is_same_v<traits::allocator_type, A>));
    EXPECT_TRUE((std::is_same_v<traits::value_type, T>));
    EXPECT_TRUE((std::is_same_v<traits::pointer, T*>));
    EXPECT_TRUE((std::is_same_v<traits::const_pointer, T const*>));
    EXPECT_TRUE((std::is_same_v<traits::void_pointer, void*>));
    EXPECT_TRUE((std::is_same_v<traits::const_void_pointer, void const*>));
    EXPECT_TRUE((std::is_same_v<traits::difference_type, std::ptrdiff_t>));
    EXPECT_TRUE((std::is_same_v<traits::size_type, std::size_t>));
    EXPECT_TRUE((std::is_same_v<traits::propagate_on_container_copy_assignment, std::false_type>));
    EXPECT_TRUE((std::is_same_v<traits::propagate_on_container_move_assignment, std::false_type>));
    EXPECT_TRUE((std::is_same_v<traits::propagate_on_container_swap, std::false_type>));
    EXPECT_TRUE((std::is_same_v<traits::rebind_alloc<U>, ddc::AlignedAllocator<U, 64>>));
    EXPECT_TRUE((std::is_same_v<traits::is_always_equal, std::true_type>));
}

TEST(AlignedAllocatorTest, DefaultConstructor)
{
    EXPECT_TRUE(std::is_default_constructible_v<A>);
}

TEST(AlignedAllocatorTest, CopyConstructor)
{
    EXPECT_TRUE(std::is_copy_constructible_v<A>);
}

TEST(AlignedAllocatorTest, RebindCopyConstructor)
{
    EXPECT_TRUE((std::is_constructible_v<A, B const&>));
}

TEST(AlignedAllocatorTest, MoveConstructor)
{
    EXPECT_TRUE(std::is_move_constructible_v<A>);
}

TEST(AlignedAllocatorTest, RebindMoveConstructor)
{
    EXPECT_TRUE((std::is_constructible_v<A, B&&>));
}

TEST(AlignedAllocatorTest, CopyAssignment)
{
    EXPECT_TRUE(std::is_copy_assignable_v<A>);
}

TEST(AlignedAllocatorTest, RebindCopyAssignment)
{
    EXPECT_TRUE((std::is_assignable_v<A, B>));
}

TEST(AlignedAllocatorTest, MoveAssignment)
{
    EXPECT_TRUE(std::is_move_assignable_v<A>);
}

TEST(AlignedAllocatorTest, RebindMoveAssignment)
{
    EXPECT_TRUE((std::is_assignable_v<A, B&&>));
}

TEST(AlignedAllocatorTest, Comparison)
{
    constexpr A a1;
    constexpr A a2;
    EXPECT_TRUE(a1 == a2);
    EXPECT_FALSE((a1 != a2));
}

TEST(AlignedAllocatorTest, RebindComparison)
{
    constexpr B b;
    constexpr A a(b);
    EXPECT_EQ(b, B(a));
    EXPECT_EQ(a, A(b));
}

TEST(AlignedAllocatorTest, AllocationType)
{
    EXPECT_TRUE((std::is_same_v<
                 decltype(std::declval<A>().allocate(1)),
                 std::allocator_traits<A>::pointer>));
}
