// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

using T = double;
using A = AlignedAllocator<T, 64>;
using U = char;
using B = std::allocator_traits<A>::rebind_alloc<U>;

TEST(AlignedAllocatorTest, Traits)
{
    using traits = std::allocator_traits<A>;
    ASSERT_TRUE((std::is_same_v<traits::allocator_type, A>));
    ASSERT_TRUE((std::is_same_v<traits::value_type, T>));
    ASSERT_TRUE((std::is_same_v<traits::pointer, T*>));
    ASSERT_TRUE((std::is_same_v<traits::const_pointer, T const*>));
    ASSERT_TRUE((std::is_same_v<traits::void_pointer, void*>));
    ASSERT_TRUE((std::is_same_v<traits::const_void_pointer, void const*>));
    ASSERT_TRUE((std::is_same_v<traits::difference_type, std::ptrdiff_t>));
    ASSERT_TRUE((std::is_same_v<traits::size_type, std::size_t>));
    ASSERT_TRUE((std::is_same_v<traits::propagate_on_container_copy_assignment, std::false_type>));
    ASSERT_TRUE((std::is_same_v<traits::propagate_on_container_move_assignment, std::false_type>));
    ASSERT_TRUE((std::is_same_v<traits::propagate_on_container_swap, std::false_type>));
    ASSERT_TRUE((std::is_same_v<traits::rebind_alloc<U>, AlignedAllocator<U, 64>>));
    ASSERT_TRUE((std::is_same_v<traits::is_always_equal, std::true_type>));
}

TEST(AlignedAllocatorTest, DefaultConstructor)
{
    ASSERT_TRUE(std::is_default_constructible_v<A>);
}

TEST(AlignedAllocatorTest, CopyConstructor)
{
    ASSERT_TRUE(std::is_copy_constructible_v<A>);
}

TEST(AlignedAllocatorTest, RebindCopyConstructor)
{
    ASSERT_TRUE((std::is_constructible_v<A, B const&>));
}

TEST(AlignedAllocatorTest, MoveConstructor)
{
    ASSERT_TRUE(std::is_move_constructible_v<A>);
}

TEST(AlignedAllocatorTest, RebindMoveConstructor)
{
    ASSERT_TRUE((std::is_constructible_v<A, B&&>));
}

TEST(AlignedAllocatorTest, CopyAssignment)
{
    ASSERT_TRUE(std::is_copy_assignable_v<A>);
}

TEST(AlignedAllocatorTest, RebindCopyAssignment)
{
    ASSERT_TRUE((std::is_assignable_v<A, B>));
}

TEST(AlignedAllocatorTest, MoveAssignment)
{
    ASSERT_TRUE(std::is_move_assignable_v<A>);
}

TEST(AlignedAllocatorTest, RebindMoveAssignment)
{
    ASSERT_TRUE((std::is_assignable_v<A, B&&>));
}

TEST(AlignedAllocatorTest, Comparison)
{
    constexpr A a1;
    constexpr A a2;
    ASSERT_TRUE(a1 == a2);
    ASSERT_FALSE((a1 != a2));
}

TEST(AlignedAllocatorTest, RebindComparison)
{
    constexpr B b;
    constexpr A a(b);
    ASSERT_EQ(b, B(a));
    ASSERT_EQ(a, A(b));
}

TEST(AlignedAllocatorTest, AllocationType)
{
    ASSERT_TRUE((std::is_same_v<
                 decltype(std::declval<A>().allocate(1)),
                 std::allocator_traits<A>::pointer>));
}
