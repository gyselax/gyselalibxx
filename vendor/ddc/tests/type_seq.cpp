#include <memory>
#include <type_traits>

#include <gtest/gtest.h>

#include "type_seq.h"

using namespace detail;

struct a;
struct b;
struct c;
struct y;
struct z;

TEST(TypeSeqTest, rank)
{
    using A = TypeSeq<a, b, c>;
    EXPECT_EQ((type_seq_rank_v<a, A>), 0);
    EXPECT_EQ((type_seq_rank_v<b, A>), 1);
    EXPECT_EQ((type_seq_rank_v<c, A>), 2);
}

TEST(TypeSeqTest, element)
{
    using A = TypeSeq<a, b, c>;
    EXPECT_TRUE((std::is_same_v<a, type_seq_element_t<0, A>>));
    EXPECT_TRUE((std::is_same_v<b, type_seq_element_t<1, A>>));
    EXPECT_TRUE((std::is_same_v<c, type_seq_element_t<2, A>>));
}

TEST(TypeSeqTest, same_tags)
{
    using A = TypeSeq<a, b, c>;
    using B = TypeSeq<z, c, y>;
    using C = TypeSeq<c, b, a>;
    EXPECT_FALSE((type_seq_same_v<A, B>));
    EXPECT_FALSE((type_seq_same_v<B, A>));
    EXPECT_TRUE((type_seq_same_v<A, C>));
    EXPECT_TRUE((type_seq_same_v<C, A>));
    EXPECT_FALSE((type_seq_same_v<B, C>));
    EXPECT_FALSE((type_seq_same_v<C, B>));
}

TEST(TypeSeqTest, remove)
{
    using A = TypeSeq<a, b, c>;
    using B = TypeSeq<z, c, y>;
    using R = type_seq_remove_t<A, B>;
    using ExpectedR = TypeSeq<a, b>;
    EXPECT_TRUE((type_seq_same_v<R, ExpectedR>));
}

TEST(TypeSeqTest, merge)
{
    using A = TypeSeq<a, b, c>;
    using B = TypeSeq<z, c, y>;
    using R = type_seq_merge_t<A, B>;
    using ExpectedR = TypeSeq<a, b, c, z, y>;
    EXPECT_TRUE((type_seq_same_v<R, ExpectedR>));
}
