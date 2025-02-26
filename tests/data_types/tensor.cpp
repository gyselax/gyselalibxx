// SPDX-License-Identifier: MIT
#include <ddc/ddc.hpp>

#include <gtest/gtest.h>

#include "indexed_tensor.hpp"
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

int dot_product(Vector<int, R_cov, Theta_cov> a, Vector<int, R, Theta> b)
{
    return ddcHelper::get<R_cov>(a) * ddcHelper::get<R>(b)
           + ddcHelper::get<Theta_cov>(a) * ddcHelper::get<Theta>(b);
}

} // namespace

TEST(TensorTest, ExplicitDotProduct)
{
    Vector<int, R_cov, Theta_cov> a;
    Vector<int, R, Theta> b;
    ddcHelper::get<R_cov>(a) = -6;
    ddcHelper::get<Theta_cov>(a) = 8;
    ddcHelper::get<R>(b) = 5;
    ddcHelper::get<Theta>(b) = 12;
    double val = dot_product(a, b);
    EXPECT_EQ(val, 66);
}

TEST(TensorTest, TensorScalarDiv)
{
    Vector<int, R_cov, Theta_cov> a;
    ddcHelper::get<R_cov>(a) = -6;
    ddcHelper::get<Theta_cov>(a) = 8;
    a /= 2;
    int val = ddcHelper::get<R_cov>(a);
    EXPECT_EQ(val, -3);
    val = ddcHelper::get<Theta_cov>(a);
    EXPECT_EQ(val, 4);
}

TEST(TensorTest, TensorAdd)
{
    Vector<int, R_cov, Theta_cov> a;
    Vector<int, R_cov, Theta_cov> b;
    ddcHelper::get<R_cov>(a) = -6;
    ddcHelper::get<Theta_cov>(a) = 8;
    ddcHelper::get<R_cov>(b) = 6;
    ddcHelper::get<Theta_cov>(b) = -16;
    a += b;
    int val = ddcHelper::get<R_cov>(a);
    EXPECT_EQ(val, 0);
    val = ddcHelper::get<Theta_cov>(a);
    EXPECT_EQ(val, -8);
}

TEST(TensorTest, TensorMinus)
{
    Vector<int, R_cov, Theta_cov> a;
    Vector<int, R_cov, Theta_cov> b;
    ddcHelper::get<R_cov>(a) = -6;
    ddcHelper::get<Theta_cov>(a) = 8;
    ddcHelper::get<R_cov>(b) = 6;
    ddcHelper::get<Theta_cov>(b) = -16;
    a -= b;
    int val = ddcHelper::get<R_cov>(a);
    EXPECT_EQ(val, -12);
    val = ddcHelper::get<Theta_cov>(a);
    EXPECT_EQ(val, 24);
}

TEST(TensorTools, CharOccurences)
{
    using CharTypeSeq = ddc::detail::TypeSeq<
            std::integral_constant<char, 'i'>,
            std::integral_constant<char, 'j'>,
            std::integral_constant<char, 'i'>,
            std::integral_constant<char, 'k'>,
            std::integral_constant<char, 'k'>,
            std::integral_constant<char, 'k'>>;
    std::size_t n_i = char_occurences_v<'i', CharTypeSeq>;
    std::size_t n_j = char_occurences_v<'j', CharTypeSeq>;
    std::size_t n_k = char_occurences_v<'k', CharTypeSeq>;
    EXPECT_EQ(n_i, 2);
    EXPECT_EQ(n_j, 1);
    EXPECT_EQ(n_k, 3);
}

TEST(TensorTools, TypeSeqUnique)
{
    using CharTypeSeq = ddc::detail::TypeSeq<
            std::integral_constant<char, 'i'>,
            std::integral_constant<char, 'j'>,
            std::integral_constant<char, 'i'>,
            std::integral_constant<char, 'k'>,
            std::integral_constant<char, 'k'>,
            std::integral_constant<char, 'k'>>;
    using ExpectedTypeSeq = ddc::detail::TypeSeq<
            std::integral_constant<char, 'i'>,
            std::integral_constant<char, 'j'>,
            std::integral_constant<char, 'k'>>;
    using UniqueCharTypeSeq = type_seq_unique_t<CharTypeSeq>;
    static_assert(std::is_same_v<UniqueCharTypeSeq, ExpectedTypeSeq>);
}

TEST(TensorTools, RepeatedIndices)
{
    using IdxSet = VectorIndexSet<R, Theta>;
    using IdxSet_cov = VectorIndexSet<R_cov, Theta_cov>;
    using TypeSeqVectorIdxMap = ddc::detail::TypeSeq<
            VectorIndexIdMap<'i', IdxSet>,
            VectorIndexIdMap<'j', IdxSet>,
            VectorIndexIdMap<'j', IdxSet_cov>,
            VectorIndexIdMap<'k', IdxSet_cov>,
            VectorIndexIdMap<'k', IdxSet>,
            VectorIndexIdMap<'l', IdxSet_cov>>;
    using NonRepeatedIdxs = non_repeated_indices_t<TypeSeqVectorIdxMap>;
    using UniqueIdxs = unique_indices_t<TypeSeqVectorIdxMap>;
    static_assert(std::is_same_v<
                  NonRepeatedIdxs,
                  ddc::detail::TypeSeq<
                          VectorIndexIdMap<'i', IdxSet>,
                          VectorIndexIdMap<'l', IdxSet_cov>>>);
    static_assert(std::is_same_v<
                  UniqueIdxs,
                  ddc::detail::TypeSeq<
                          VectorIndexIdMap<'i', IdxSet>,
                          VectorIndexIdMap<'j', IdxSet>,
                          VectorIndexIdMap<'k', IdxSet>,
                          VectorIndexIdMap<'l', IdxSet>>>);
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

TEST(TensorTools, TensorIndexMap)
{
    using IdxSet = VectorIndexSet<R, Theta>;
    using IdxSet_cov = VectorIndexSet<R_cov, Theta_cov>;
    using TestIndexSet = ddc::detail::TypeSeq<IdxSet, IdxSet_cov, IdxSet, IdxSet>;
    using GlobalTensorIndexIdMap = ddc::detail::TypeSeq<
            VectorIndexIdMap<'i', IdxSet>,
            VectorIndexIdMap<'j', IdxSet_cov>,
            VectorIndexIdMap<'j', IdxSet>,
            VectorIndexIdMap<'k', IdxSet>>;
    using TestIndexElement0 = get_nth_tensor_index_element_from_map_t<0, GlobalTensorIndexIdMap>;
    static_assert(TestIndexElement0::rank() == 4);
    static_assert(
            std::is_same_v<TestIndexElement0, TensorIndexElement<TestIndexSet, R, R_cov, R, R>>);
    using TestIndexElement1 = get_nth_tensor_index_element_from_map_t<1, GlobalTensorIndexIdMap>;
    static_assert(TestIndexElement1::rank() == 4);
    static_assert(std::is_same_v<
                  TestIndexElement1,
                  TensorIndexElement<TestIndexSet, R, R_cov, R, Theta>>);
    using TestIndexElement2 = get_nth_tensor_index_element_from_map_t<2, GlobalTensorIndexIdMap>;
    static_assert(TestIndexElement2::rank() == 4);
    static_assert(std::is_same_v<
                  TestIndexElement2,
                  TensorIndexElement<TestIndexSet, R, Theta_cov, Theta, R>>);
    using TestIndexElement3 = get_nth_tensor_index_element_from_map_t<3, GlobalTensorIndexIdMap>;
    static_assert(TestIndexElement3::rank() == 4);
    static_assert(std::is_same_v<
                  TestIndexElement3,
                  TensorIndexElement<TestIndexSet, R, Theta_cov, Theta, Theta>>);
    using TestIndexElement4 = get_nth_tensor_index_element_from_map_t<4, GlobalTensorIndexIdMap>;
    static_assert(TestIndexElement4::rank() == 4);
    static_assert(std::is_same_v<
                  TestIndexElement4,
                  TensorIndexElement<TestIndexSet, Theta, R_cov, R, R>>);
    using TestIndexElement5 = get_nth_tensor_index_element_from_map_t<5, GlobalTensorIndexIdMap>;
    static_assert(TestIndexElement5::rank() == 4);
    static_assert(std::is_same_v<
                  TestIndexElement5,
                  TensorIndexElement<TestIndexSet, Theta, R_cov, R, Theta>>);
    using TestIndexElement6 = get_nth_tensor_index_element_from_map_t<6, GlobalTensorIndexIdMap>;
    static_assert(TestIndexElement6::rank() == 4);
    static_assert(std::is_same_v<
                  TestIndexElement6,
                  TensorIndexElement<TestIndexSet, Theta, Theta_cov, Theta, R>>);
    using TestIndexElement7 = get_nth_tensor_index_element_from_map_t<7, GlobalTensorIndexIdMap>;
    static_assert(TestIndexElement7::rank() == 4);
    static_assert(std::is_same_v<
                  TestIndexElement7,
                  TensorIndexElement<TestIndexSet, Theta, Theta_cov, Theta, Theta>>);
}

TEST(TensorTools, ExtractElement)
{
    using IdxSet = VectorIndexSet<R, Theta>;
    using IdxSet_cov = VectorIndexSet<R_cov, Theta_cov>;
    using GlobalTensorIndexSet = ddc::detail::TypeSeq<IdxSet, IdxSet_cov, IdxSet, IdxSet>;
    using LocalTensorIndexSet = ddc::detail::TypeSeq<IdxSet, IdxSet_cov>;
    using GlobalTensorIndexIdMap = ddc::detail::TypeSeq<
            VectorIndexIdMap<'i', IdxSet>,
            VectorIndexIdMap<'j', IdxSet_cov>,
            VectorIndexIdMap<'j', IdxSet>,
            VectorIndexIdMap<'k', IdxSet>>;
    using LocalTensorIndexIdMap = ddc::detail::
            TypeSeq<VectorIndexIdMap<'i', IdxSet>, VectorIndexIdMap<'j', IdxSet_cov>>;
    using GlobalTensorIndexElement = TensorIndexElement<GlobalTensorIndexSet, R, R_cov, R, R>;
    using LocalTensorIndexElement = extract_sub_tensor_element_t<
            GlobalTensorIndexIdMap,
            LocalTensorIndexIdMap,
            GlobalTensorIndexElement>;
    static_assert(LocalTensorIndexElement::rank() == 2);
    static_assert(std::is_same_v<
                  LocalTensorIndexElement,
                  TensorIndexElement<LocalTensorIndexSet, R, R_cov>>);
    using GlobalTensorIndexElementTest1
            = TensorIndexElement<GlobalTensorIndexSet, R, Theta_cov, Theta, R>;
    using LocalTensorIndexElement1 = extract_sub_tensor_element_t<
            GlobalTensorIndexIdMap,
            LocalTensorIndexIdMap,
            GlobalTensorIndexElementTest1>;
    static_assert(LocalTensorIndexElement1::rank() == 2);
    static_assert(std::is_same_v<
                  LocalTensorIndexElement1,
                  TensorIndexElement<LocalTensorIndexSet, R, Theta_cov>>);
    using GlobalTensorIndexElementTest2
            = TensorIndexElement<GlobalTensorIndexSet, R, Theta_cov, Theta, R>;
    using LocalTensorIndexIdMap2
            = ddc::detail::TypeSeq<VectorIndexIdMap<'i', IdxSet>, VectorIndexIdMap<'j', IdxSet>>;
    using LocalTensorIndexSet2 = ddc::detail::TypeSeq<IdxSet, IdxSet>;
    using LocalTensorIndexElement2 = extract_sub_tensor_element_t<
            GlobalTensorIndexIdMap,
            LocalTensorIndexIdMap2,
            GlobalTensorIndexElementTest2>;
    static_assert(LocalTensorIndexElement2::rank() == 2);
    static_assert(std::is_same_v<
                  LocalTensorIndexElement2,
                  TensorIndexElement<LocalTensorIndexSet2, R, Theta>>);
    static_assert(std::is_same_v<LocalTensorIndexElement2::tensor_index_set, LocalTensorIndexSet2>);
}

TEST(TensorTools, IndexedTensorBuild)
{
    using IdxSet = VectorIndexSet<R, Theta>;
    using IdxSet_cov = VectorIndexSet<R_cov, Theta_cov>;
    using Tensor3D = Tensor<int, IdxSet, IdxSet_cov, IdxSet>;
    Tensor3D G;
    IndexedTensor G_idx = index<'i', 'j', 'k'>(G);
    static_assert(std::is_same_v<
                  typename decltype(G_idx)::index_pattern,
                  ddc::detail::TypeSeq<
                          VectorIndexIdMap<'i', IdxSet>,
                          VectorIndexIdMap<'j', IdxSet_cov>,
                          VectorIndexIdMap<'k', IdxSet>>>);
}

TEST(TensorTest, Mul)
{
    using Tensor2D_A = Tensor<int, VectorIndexSet<R, Theta>, VectorIndexSet<R, Theta>>;
    using Tensor2D_B
            = Tensor<int, VectorIndexSet<R_cov, Theta_cov>, VectorIndexSet<R_cov, Theta_cov>>;
    using Tensor2D_C = Tensor<int, VectorIndexSet<R, Theta>, VectorIndexSet<R_cov, Theta_cov>>;
    Tensor2D_A A;
    Tensor2D_B B;
    ddcHelper::get<R, R>(A) = 1;
    ddcHelper::get<R, Theta>(A) = 0;
    ddcHelper::get<Theta, R>(A) = 2;
    ddcHelper::get<Theta, Theta>(A) = 4;
    ddcHelper::get<R_cov, R_cov>(B) = 6;
    ddcHelper::get<R_cov, Theta_cov>(B) = 8;
    ddcHelper::get<Theta_cov, R_cov>(B) = 4;
    ddcHelper::get<Theta_cov, Theta_cov>(B) = 5;
    Tensor2D_C C = tensor_mul(index<'i', 'j'>(A), index<'j', 'k'>(B));
    int val = ddcHelper::get<R, R_cov>(C);
    EXPECT_EQ(val, 6);
    val = ddcHelper::get<R, Theta_cov>(C);
    EXPECT_EQ(val, 8);
    val = ddcHelper::get<Theta, R_cov>(C);
    EXPECT_EQ(val, 28);
    val = ddcHelper::get<Theta, Theta_cov>(C);
    EXPECT_EQ(val, 36);
    Tensor2D_C D = tensor_mul(index<'i', 'j'>(A), index<'k', 'j'>(B));
    val = ddcHelper::get<R, R_cov>(D);
    EXPECT_EQ(val, 6);
    val = ddcHelper::get<R, Theta_cov>(D);
    EXPECT_EQ(val, 4);
    val = ddcHelper::get<Theta, R_cov>(D);
    EXPECT_EQ(val, 44);
    val = ddcHelper::get<Theta, Theta_cov>(D);
    EXPECT_EQ(val, 28);
    int E = tensor_mul(index<'i', 'j'>(A), index<'j', 'i'>(B));
    EXPECT_EQ(E, 42);
}
