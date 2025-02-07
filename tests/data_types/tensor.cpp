// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "tensor.hpp"
#include "tensor_index_tools.hpp"
#include "vector_index_tools.hpp"

using namespace tensor_tools;

namespace {
struct R_cov;
struct Theta_cov;

struct R
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = R_cov;
};
struct Theta
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = false;
    static bool constexpr IS_CONTRAVARIANT = true;
    using Dual = Theta_cov;
};

struct R_cov
{
    static bool constexpr PERIODIC = false;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    using Dual = R;
};
struct Theta_cov
{
    static bool constexpr PERIODIC = true;
    static bool constexpr IS_COVARIANT = true;
    static bool constexpr IS_CONTRAVARIANT = false;
    using Dual = Theta;
};

int dot_product(
        tensor_tools::Vector<int, R_cov, Theta_cov> a,
        tensor_tools::Vector<int, R, Theta> b)
{
    return ddcHelper::get<R_cov>(a) * ddcHelper::get<R>(b)
           + ddcHelper::get<Theta_cov>(a) * ddcHelper::get<Theta>(b);
}

TEST(TensorTest, ExplicitDotProduct)
{
    tensor_tools::Vector<int, R_cov, Theta_cov> a;
    tensor_tools::Vector<int, R, Theta> b;
    ddcHelper::get<R_cov>(a) = -6;
    ddcHelper::get<Theta_cov>(a) = 8;
    ddcHelper::get<R>(b) = 5;
    ddcHelper::get<Theta>(b) = 12;
    int val = dot_product(a, b);
    EXPECT_EQ(val, 66);
}

TEST(TensorTest, TensorScalarMul)
{
    tensor_tools::Vector<int, R_cov, Theta_cov> a;
    ddcHelper::get<R_cov>(a) = -6;
    ddcHelper::get<Theta_cov>(a) = 8;
    a *= 4;
    int val = ddcHelper::get<R_cov>(a);
    EXPECT_EQ(val, -24);
    val = ddcHelper::get<Theta_cov>(a);
    EXPECT_EQ(val, 32);
}

TEST(TensorTools, TensorIndexSet)
{
    using IdxSet = VectorIndexSet<R, Theta>;
    using IdxSet_cov = VectorIndexSet<R_cov, Theta_cov>;
    using TestIndexSet = ddc::detail::TypeSeq<IdxSet, IdxSet, IdxSet_cov>;
    using TestIndexElement0 = get_nth_tensor_index_element_t<0, TestIndexSet>;
    static_assert(TestIndexElement0::rank() == 3);
    static_assert(std::is_same_v<TestIndexElement0, TensorIndexElement<TestIndexSet, R, R, R_cov>>);
    static_assert(TestIndexElement0::index() == 0);
    using TestIndexElement1 = get_nth_tensor_index_element_t<1, TestIndexSet>;
    static_assert(TestIndexElement1::rank() == 3);
    static_assert(
            std::is_same_v<TestIndexElement1, TensorIndexElement<TestIndexSet, R, R, Theta_cov>>);
    static_assert(TestIndexElement1::index() == 1);
    using TestIndexElement2 = get_nth_tensor_index_element_t<2, TestIndexSet>;
    static_assert(
            std::is_same_v<TestIndexElement2, TensorIndexElement<TestIndexSet, R, Theta, R_cov>>);
    static_assert(TestIndexElement2::index() == 2);
    using TestIndexElement3 = get_nth_tensor_index_element_t<3, TestIndexSet>;
    static_assert(std::is_same_v<
                  TestIndexElement3,
                  TensorIndexElement<TestIndexSet, R, Theta, Theta_cov>>);
    static_assert(TestIndexElement3::index() == 3);
    using TestIndexElement4 = get_nth_tensor_index_element_t<4, TestIndexSet>;
    static_assert(
            std::is_same_v<TestIndexElement4, TensorIndexElement<TestIndexSet, Theta, R, R_cov>>);
    static_assert(TestIndexElement4::index() == 4);
    using TestIndexElement5 = get_nth_tensor_index_element_t<5, TestIndexSet>;
    static_assert(std::is_same_v<
                  TestIndexElement5,
                  TensorIndexElement<TestIndexSet, Theta, R, Theta_cov>>);
    static_assert(TestIndexElement5::index() == 5);
    using TestIndexElement6 = get_nth_tensor_index_element_t<6, TestIndexSet>;
    static_assert(std::is_same_v<
                  TestIndexElement6,
                  TensorIndexElement<TestIndexSet, Theta, Theta, R_cov>>);
    static_assert(TestIndexElement6::index() == 6);
    using TestIndexElement7 = get_nth_tensor_index_element_t<7, TestIndexSet>;
    static_assert(std::is_same_v<
                  TestIndexElement7,
                  TensorIndexElement<TestIndexSet, Theta, Theta, Theta_cov>>);
    static_assert(TestIndexElement7::index() == 7);
}


} // namespace
