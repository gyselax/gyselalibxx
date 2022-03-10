/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 2.0
//              Copyright (2019) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Christian R. Trott (crtrott@sandia.gov)
//
// ************************************************************************
//@HEADER
*/

#ifndef MDARRAY_INCLUDE_EXPERIMENTAL_BITS_BASIC_MDARRAY_HPP_
#define MDARRAY_INCLUDE_EXPERIMENTAL_BITS_BASIC_MDARRAY_HPP_

#include "container_policy_basic.hpp"
#include "const_wrapped_accessor_policy.hpp"

#include <experimental/__p0009_bits/macros.hpp>
#include <experimental/__p0009_bits/layout_right.hpp>
#include <experimental/__p0009_bits/extents.hpp>
#include <experimental/__p0009_bits/basic_mdspan.hpp>

namespace std {
namespace experimental {
inline namespace __mdarray_version_0 {

namespace __detail {

template <class Derived, class ExtentsIdxs>
struct _basic_mdarray_crtp_helper;

// Workaround for not being able to give explicit template parameters to lambdas in older
// versions of C++, thus making expanding parameter packs with indices more difficult
template <class Derived, size_t... ExtIdxs>
struct _basic_mdarray_crtp_helper<
  Derived, std::index_sequence<ExtIdxs...>
>
{
 protected:
  MDSPAN_FORCE_INLINE_FUNCTION Derived& __self() noexcept { return *static_cast<Derived*>(this); }
  MDSPAN_FORCE_INLINE_FUNCTION Derived const& __self() const noexcept { return *static_cast<Derived const*>(this); }
  MDSPAN_FORCE_INLINE_FUNCTION constexpr ptrdiff_t __size() const noexcept {
    return _MDSPAN_FOLD_TIMES_RIGHT((__self().map_.extents().template __extent<ExtIdxs>()), /* * ... * */ 1);
  }
  template <class ReferenceType, class IndexType, size_t N>
  MDSPAN_FORCE_INLINE_FUNCTION constexpr ReferenceType __callop(const array<IndexType, N>& indices) noexcept {
    return __self().cp_.access(__self().c_, __self().map_(indices[ExtIdxs]...));
  }
  template <class ReferenceType, class IndexType, size_t N>
  MDSPAN_FORCE_INLINE_FUNCTION constexpr ReferenceType __callop(const array<IndexType, N>& indices) const noexcept {
    return __self().cp_.access(__self().c_, __self().map_(indices[ExtIdxs]...));
  }
};


} // end namespace __detail

template <
  class ElementType,
  class Extents,
  class LayoutPolicy=layout_right,
  class ContainerPolicy=typename __detail::__container_policy_select<ElementType, LayoutPolicy, Extents>::type
>
class basic_mdarray;

template <
  class ElementType,
  ptrdiff_t... Exts,
  class LayoutPolicy,
  class ContainerPolicy
>
class basic_mdarray<
  ElementType, std::experimental::extents<Exts...>,
  LayoutPolicy, ContainerPolicy
> : __detail::_basic_mdarray_crtp_helper<
      basic_mdarray<ElementType, std::experimental::extents<Exts...>, LayoutPolicy, ContainerPolicy>,
      make_index_sequence<sizeof...(Exts)>
    >
{
private:

  using __crtp_base_t = __detail::_basic_mdarray_crtp_helper<
    basic_mdarray<ElementType, std::experimental::extents<Exts...>, LayoutPolicy, ContainerPolicy>,
    make_index_sequence<sizeof...(Exts)>
  >;

public:
  using element_type = ElementType;
  using extents_type = std::experimental::extents<Exts...>;
  using layout_type = LayoutPolicy;
  using mapping_type = typename layout_type::template mapping<extents_type>; // TODO @proposal-bug typo in synopsis
  using value_type = remove_cv_t<element_type>;
  using index_type = ptrdiff_t;
  using difference_type = ptrdiff_t;
  using container_policy_type = ContainerPolicy;
  using container_type = typename container_policy_type::container_type;
  using pointer = typename container_policy_type::pointer; // TODO @proposal-bug this is misspelled in the synopsis
  using const_pointer = typename container_policy_type::const_pointer;
  using reference = typename container_policy_type::reference;
  using const_reference = typename container_policy_type::const_reference;
  using view_type =
    basic_mdspan<element_type, extents_type, layout_type, container_policy_type>;
  using const_view_type =
    basic_mdspan<const element_type, extents_type, layout_type,
      __detail::__const_wrapped_accessor_policy<container_policy_type>
    >;

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr basic_mdarray() noexcept(std::is_nothrow_default_constructible<container_type>::value) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr basic_mdarray(basic_mdarray const&) noexcept(std::is_nothrow_copy_constructible<container_type>::value) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr basic_mdarray(basic_mdarray&&) noexcept(std::is_nothrow_move_constructible<container_type>::value) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED basic_mdarray& operator=(basic_mdarray&&) noexcept(std::is_nothrow_move_assignable<container_type>::value) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED basic_mdarray& operator=(basic_mdarray const&) noexcept(std::is_nothrow_copy_assignable<container_type>::value) = default;

  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~basic_mdarray() noexcept(std::is_nothrow_destructible<container_type>::value) = default;

  // TODO noexcept clause
  MDSPAN_TEMPLATE_REQUIRES(
    class... IndexType,
    /* requires */ (
      _MDSPAN_FOLD_AND(_MDSPAN_TRAIT(is_convertible, IndexType, index_type) /* && ... */) &&
      (sizeof...(IndexType) == extents_type::rank_dynamic()) &&
      _MDSPAN_TRAIT(is_constructible, mapping_type, extents_type) &&
      _MDSPAN_TRAIT(is_default_constructible, container_policy_type)
      // TODO constraint on create without allocator being available, if we don't change to CP owning the allocator
    )
  )
  MDSPAN_INLINE_FUNCTION
  constexpr explicit
  basic_mdarray(IndexType... dynamic_extents)
    : cp_(),
      map_(extents_type(dynamic_extents...)),
      c_(cp_.create(map_.required_span_size()))
  { }

  // TODO noexcept specification
  MDSPAN_FUNCTION_REQUIRES(
    (MDSPAN_INLINE_FUNCTION constexpr explicit),
    basic_mdarray, (const mapping_type& m), noexcept,
    /* requires */ (_MDSPAN_TRAIT(is_default_constructible, container_policy_type))
  ) : cp_(),
      map_(m),
      c_(cp_.create(map_.required_span_size()))
  { }

  // TODO noexcept specification
  // TODO @proposal-bug add this to the proposal?
  MDSPAN_FUNCTION_REQUIRES(
    (MDSPAN_INLINE_FUNCTION constexpr explicit),
    basic_mdarray, (mapping_type&& m), noexcept,
    /* requires */ (_MDSPAN_TRAIT(is_default_constructible, container_policy_type))
  ) : cp_(),
      map_(std::move(m)),
      c_(cp_.create(map_.required_span_size()))
  { }

  MDSPAN_INLINE_FUNCTION
  constexpr
  basic_mdarray(mapping_type const& m, container_policy_type const& cp)
    : cp_(cp),
      map_(m),
      c_(cp_.create(map_.required_span_size()))
  { }

  // TODO noexcept specification
  MDSPAN_TEMPLATE_REQUIRES(
    class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherCP,
    /* requires */ (
      _MDSPAN_TRAIT(is_convertible, typename OtherLayoutPolicy::template mapping<OtherExtents>, mapping_type) &&
      _MDSPAN_TRAIT(is_constructible, container_policy_type, OtherCP const&) &&
      _MDSPAN_TRAIT(is_constructible, container_type, typename OtherCP::container_type const&) &&
      _MDSPAN_TRAIT(is_convertible, OtherExtents, extents_type)
    )
  )
  constexpr basic_mdarray(
    basic_mdarray<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherCP> const& other
  ) noexcept
    : cp_(other.cp_),
      map_(other.map_),
      c_(other.c_)
  { }

  // TODO noexcept specification
  MDSPAN_TEMPLATE_REQUIRES(
    class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherCP,
    /* requires */ (
      _MDSPAN_TRAIT(is_convertible, typename OtherLayoutPolicy::template mapping<OtherExtents>, mapping_type) &&
      _MDSPAN_TRAIT(is_constructible, container_policy_type, OtherCP&&) &&
      _MDSPAN_TRAIT(is_constructible, container_type, typename OtherCP::container_type&&) &&
      _MDSPAN_TRAIT(is_convertible, OtherExtents, extents_type)
    )
  )
  constexpr basic_mdarray(
    basic_mdarray<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherCP>&& other
  ) noexcept
    : cp_(std::move(other.cp_)),
      map_(std::move(other.map_)),
      c_(std::move(other.c_))
  { }

  // TODO noexcept specification
  MDSPAN_TEMPLATE_REQUIRES(
    class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherCP,
    /* requires */ (
      _MDSPAN_TRAIT(is_convertible, typename OtherLayoutPolicy::template mapping<OtherExtents>, mapping_type) &&
      _MDSPAN_TRAIT(is_assignable, container_policy_type, OtherCP const&) &&
      _MDSPAN_TRAIT(is_assignable, container_type, typename OtherCP::container_type const&) &&
      _MDSPAN_TRAIT(is_convertible, OtherExtents, extents_type)
    )
  )
  _MDSPAN_CONSTEXPR_14 basic_mdarray& operator=(
    basic_mdarray<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherCP> const& other
  ) noexcept
  {
    cp_ = other.cp_;
    map_ = other.map_;
    c_ = other.c_;
    return *this;
  }

  // TODO noexcept specification
  MDSPAN_TEMPLATE_REQUIRES(
    class OtherElementType, class OtherExtents, class OtherLayoutPolicy, class OtherCP,
    /* requires */ (
      _MDSPAN_TRAIT(is_convertible, typename OtherLayoutPolicy::template mapping<OtherExtents>, mapping_type) &&
      _MDSPAN_TRAIT(is_assignable, container_policy_type, OtherCP&&) &&
      _MDSPAN_TRAIT(is_assignable, container_type, typename OtherCP::container_type&&) &&
      _MDSPAN_TRAIT(is_convertible, OtherExtents, extents_type)
    )
  )
  _MDSPAN_CONSTEXPR_14 basic_mdarray& operator=(
    basic_mdarray<OtherElementType, OtherExtents, OtherLayoutPolicy, OtherCP>&& other
  ) noexcept
  {
    cp_ = std::move(other.cp_);
    map_ = std::move(other.map_);
    c_ = std::move(other.c_);
    return *this;
  }

  //==========================================================================

  MDSPAN_TEMPLATE_REQUIRES(
    class Index,
    /* requires */ (
      _MDSPAN_TRAIT(is_convertible, Index, index_type) &&
        sizeof...(Exts) == 1
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  _MDSPAN_CONSTEXPR_14 reference operator[](Index idx) noexcept
  {
    return cp_.access(c_, map_(index_type(idx)));
  }

  MDSPAN_TEMPLATE_REQUIRES(
    class Index,
    /* requires */ (
      _MDSPAN_TRAIT(is_convertible, Index, index_type) &&
      sizeof...(Exts) == 1
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr reference operator[](Index idx) const noexcept
  {
    return cp_.access(c_, map_(index_type(idx)));
  }

  MDSPAN_TEMPLATE_REQUIRES(
    class... IndexType,
    /* requires */ (
      _MDSPAN_FOLD_AND(_MDSPAN_TRAIT(is_convertible, IndexType, index_type) /* && ... */) &&
      sizeof...(Exts) == extents_type::rank()
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  _MDSPAN_CONSTEXPR_14 reference operator()(IndexType... indices) noexcept
  {
    return cp_.access(c_, map_(index_type(indices)...));
  }

  MDSPAN_TEMPLATE_REQUIRES(
    class... IndexType,
    /* requires */ (
      _MDSPAN_FOLD_AND(_MDSPAN_TRAIT(is_convertible, IndexType, index_type) /* && ... */) &&
        sizeof...(Exts) == extents_type::rank()
    )
  )
  MDSPAN_FORCE_INLINE_FUNCTION
  constexpr const_reference operator()(IndexType... indices) const noexcept
  {
    return cp_.access(c_, map_(index_type(indices)...));
  }


  // TODO array version of dereference op

  //==========================================================================

  MDSPAN_INLINE_FUNCTION static constexpr int rank() noexcept { return extents_type::rank(); }
  MDSPAN_INLINE_FUNCTION static constexpr int rank_dynamic() noexcept { return extents_type::rank_dynamic(); }
  MDSPAN_INLINE_FUNCTION static constexpr index_type static_extent(size_t r) noexcept { return extents_type::static_extent(r); }

  MDSPAN_INLINE_FUNCTION constexpr extents_type extents() const noexcept { return map_.extents(); };
  MDSPAN_INLINE_FUNCTION constexpr index_type extent(size_t r) const noexcept { return map_.extents().extent(r); };
  // TODO basic_mdarray.size()
  MDSPAN_INLINE_FUNCTION constexpr index_type size() const noexcept {
    return this->__crtp_base_t::template __size();
  };

  MDSPAN_INLINE_FUNCTION _MDSPAN_CONSTEXPR_14 index_type unique_size() const noexcept {
    if(map_.is_unique()) {
      return size();
    }
    else if(map_.is_contiguous()) {
      return map_.required_span_size();
    }
    else {
      // ??? guess, for now, until this gets fixed in the proposal ???
      return map_.required_span_size();
    }
  }

  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_unique() noexcept { return mapping_type::is_always_unique(); };
  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_contiguous() noexcept { return mapping_type::is_always_contiguous(); };
  MDSPAN_INLINE_FUNCTION static constexpr bool is_always_strided() noexcept { return mapping_type::is_always_strided(); };

  MDSPAN_INLINE_FUNCTION constexpr mapping_type mapping() const noexcept { return map_; };
  MDSPAN_INLINE_FUNCTION constexpr bool is_unique() const noexcept { return map_.is_unique(); };
  MDSPAN_INLINE_FUNCTION constexpr bool is_contiguous() const noexcept { return map_.is_contiguous(); };
  MDSPAN_INLINE_FUNCTION constexpr bool is_strided() const noexcept { return map_.is_strided(); };
  MDSPAN_INLINE_FUNCTION constexpr index_type stride(size_t r) const { return map_.stride(r); };

  //==========================================================================

  // TODO noexcept specification
  MDSPAN_INLINE_FUNCTION
  _MDSPAN_CONSTEXPR_14 view_type view() noexcept {
    return view_type(c_.data(), map_, cp_);
  }
  // TODO noexcept specification
  MDSPAN_INLINE_FUNCTION
  constexpr const const_view_type view() const noexcept {
    return const_view_type(c_.data(), map_, cp_);
  }

  // TODO noexcept specification
  MDSPAN_INLINE_FUNCTION
  _MDSPAN_CONSTEXPR_14 pointer data() noexcept {
    return cp_.data();
  }

  // TODO noexcept specification
  MDSPAN_INLINE_FUNCTION
  constexpr const_pointer data() const noexcept {
    return cp_.data();
  }

  MDSPAN_INLINE_FUNCTION
  constexpr container_policy_type container_policy() const noexcept { return cp_; }

private:

  template <class, class, class, class>
  friend class basic_mdarray;

  // TODO @proposal-bug these should be in the reverse order in the proposal, even as exposition only
  _MDSPAN_NO_UNIQUE_ADDRESS container_policy_type cp_; // TODO @proposal-bug this is misspelled in the synopsis
  _MDSPAN_NO_UNIQUE_ADDRESS mapping_type map_;
  container_type c_;

};

template <class T, ptrdiff_t... Exts>
using mdarray = basic_mdarray<T, std::experimental::extents<Exts...>>;

} // end namespace __mdarray_version_0
} // end namespace experimental
} // end namespace std

#endif //MDARRAY_INCLUDE_EXPERIMENTAL_BITS_BASIC_MDARRAY_HPP_
