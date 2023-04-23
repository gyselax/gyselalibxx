// SPDX-License-Identifier: MIT
#include <memory>
#include <type_traits>

#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

namespace {

struct a;
struct b;
struct c;
struct y;
struct z;

} // namespace

TEST(TypeSeqTest, Rank)
{
    using A = ddc::detail::TypeSeq<a, b, c>;
    EXPECT_EQ((ddc::type_seq_rank_v<a, A>), 0);
    EXPECT_EQ((ddc::type_seq_rank_v<b, A>), 1);
    EXPECT_EQ((ddc::type_seq_rank_v<c, A>), 2);
}

TEST(TypeSeqTest, Element)
{
    using A = ddc::detail::TypeSeq<a, b, c>;
    EXPECT_TRUE((std::is_same_v<a, ddc::type_seq_element_t<0, A>>));
    EXPECT_TRUE((std::is_same_v<b, ddc::type_seq_element_t<1, A>>));
    EXPECT_TRUE((std::is_same_v<c, ddc::type_seq_element_t<2, A>>));
}

TEST(TypeSeqTest, SameTags)
{
    using A = ddc::detail::TypeSeq<a, b, c>;
    using B = ddc::detail::TypeSeq<z, c, y>;
    using C = ddc::detail::TypeSeq<c, b, a>;
    EXPECT_FALSE((ddc::type_seq_same_v<A, B>));
    EXPECT_FALSE((ddc::type_seq_same_v<B, A>));
    EXPECT_TRUE((ddc::type_seq_same_v<A, C>));
    EXPECT_TRUE((ddc::type_seq_same_v<C, A>));
    EXPECT_FALSE((ddc::type_seq_same_v<B, C>));
    EXPECT_FALSE((ddc::type_seq_same_v<C, B>));
}

TEST(TypeSeqTest, Remove)
{
    using A = ddc::detail::TypeSeq<a, b, c>;
    using B = ddc::detail::TypeSeq<z, c, y>;
    using R = ddc::type_seq_remove_t<A, B>;
    using ExpectedR = ddc::detail::TypeSeq<a, b>;
    EXPECT_TRUE((ddc::type_seq_same_v<R, ExpectedR>));
}

TEST(TypeSeqTest, Merge)
{
    using A = ddc::detail::TypeSeq<a, b, c>;
    using B = ddc::detail::TypeSeq<z, c, y>;
    using R = ddc::type_seq_merge_t<A, B>;
    using ExpectedR = ddc::detail::TypeSeq<a, b, c, z, y>;
    EXPECT_TRUE((ddc::type_seq_same_v<R, ExpectedR>));
}
