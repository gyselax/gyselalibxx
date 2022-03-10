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

#ifndef MDARRAY_INCLUDE_EXPERIMENTAL_BITS_CONTAINER_POLICY_BASIC_HPP_
#define MDARRAY_INCLUDE_EXPERIMENTAL_BITS_CONTAINER_POLICY_BASIC_HPP_

#include "const_wrapped_accessor_policy.hpp"

#include <experimental/__p0009_bits/macros.hpp>
#include <experimental/__p0009_bits/dynamic_extent.hpp>
#include <experimental/__p0009_bits/extents.hpp>

#include <vector>
#include <array>
#include <memory> // std::allocator
#include <stdexcept>


// TODO @proposal-bug should the container policy own the allocator?

namespace std {
namespace experimental {
inline namespace __mdarray_version_0 {

namespace __detail {

// TODO this is inefficient with a stateful allocator. Make different class for the offset policy portion
template <class Container, class Derived>
class __container_policy_common {
public:
  using container_type = Container;
  // TODO @proposal-bug This should be value_type (both here and in P0009)
  using element_type = typename container_type::value_type;
  using pointer = typename container_type::pointer;
  using const_pointer = typename container_type::const_pointer;
  using reference = typename container_type::reference;
  using const_reference = typename container_type::const_reference;

  MDSPAN_FORCE_INLINE_FUNCTION
  static constexpr reference access(container_type& c, ptrdiff_t i)
    noexcept(noexcept(c[i]))
  {
    return c[size_t(i)];
  }
  MDSPAN_FORCE_INLINE_FUNCTION
  static constexpr const_reference access(container_type const& c, ptrdiff_t i)
    noexcept(noexcept(c[i]))
  {
    return c[size_t(i)];
  }
  MDSPAN_FORCE_INLINE_FUNCTION
  static constexpr reference access(pointer c, ptrdiff_t i) noexcept {
    return c[size_t(i)];
  }
  MDSPAN_FORCE_INLINE_FUNCTION
  static constexpr const_reference access(const_pointer c, ptrdiff_t i) noexcept {
    return c[size_t(i)];
  }

  MDSPAN_INLINE_FUNCTION
  static constexpr pointer offset(pointer p, ptrdiff_t i) noexcept {
    return p[size_t(i)];
  }
  MDSPAN_INLINE_FUNCTION
  static constexpr const_pointer offset(const_pointer p, ptrdiff_t i) noexcept {
    return p[size_t(i)];
  }

  // TODO converting constructors (and make sure they work for conversion to const offset policy
  //MDSPAN_INLINE_FUNCTION_DEFAULTED _MDSPAN_CONSTEXPR_14_DEFAULTED
  //__container_policy_common() noexcept = default;
  //MDSPAN_INLINE_FUNCTION_DEFAULTED _MDSPAN_CONSTEXPR_14_DEFAULTED
  //__container_policy_common(__container_policy_common const&) noexcept = default;
  //MDSPAN_INLINE_FUNCTION_DEFAULTED _MDSPAN_CONSTEXPR_14_DEFAULTED
  //__container_policy_common(__container_policy_common&&) noexcept = default;
  //MDSPAN_INLINE_FUNCTION_DEFAULTED _MDSPAN_CONSTEXPR_14_DEFAULTED
  //__container_policy_common& operator=(__container_policy_common const&) noexcept = default;
  //MDSPAN_INLINE_FUNCTION_DEFAULTED _MDSPAN_CONSTEXPR_14_DEFAULTED
  //__container_policy_common& operator=(__container_policy_common&&) noexcept = default;
  //MDSPAN_INLINE_FUNCTION_DEFAULTED ~__container_policy_common() noexcept = default;
  //
  //MDSPAN_INLINE_FUNCTION
  //__container_policy_common(Derived const&) noexcept { }

};

} // end namespace __detail

template <class T, class Allocator=std::allocator<T>>
class vector_container_policy
  : public __detail::__container_policy_common<
      std::vector<T, Allocator>, vector_container_policy<T, Allocator>
    >
{
private:
  using __base_t_ = __detail::__container_policy_common<
    std::vector<T, Allocator>, vector_container_policy<T, Allocator>
  >;
public:

  using offset_policy = vector_container_policy;
  using const_offset_policy = __detail::__const_wrapped_accessor_policy<__base_t_>;

  // TODO noexcept clause
  MDSPAN_FUNCTION_REQUIRES(
    (MDSPAN_INLINE_FUNCTION static constexpr typename __base_t_::container_type),
    create, (size_t n), noexcept,
    /* requires */ (
      _MDSPAN_TRAIT(is_constructible, typename __base_t_::container_type, size_t, typename __base_t_::element_type)
    )
  )
  {
    return typename __base_t_::container_type(n, typename __base_t_::element_type{});
  }

  // TODO noexcept clause
  MDSPAN_FUNCTION_REQUIRES(
    (MDSPAN_INLINE_FUNCTION static constexpr typename __base_t_::container_type),
    create, (size_t n, Allocator& alloc), noexcept,
    /* requires */ (
      _MDSPAN_TRAIT(
        is_constructible,
        typename __base_t_::container_type, size_t, typename __base_t_::element_type, Allocator&
      )
    )
  )
  {
    return typename __base_t_::container_type(n, typename __base_t_::element_type{}, alloc);
  }
};

template <class T, size_t N>
class array_container_policy
 : public __detail::__container_policy_common<
    std::array<T, N>, array_container_policy<T, N>
   >
{
private:
  using __base_t_ = __detail::__container_policy_common<
    std::array<T, N>, array_container_policy<T, N>
  >;
public:

  using offset_policy = array_container_policy;
  using const_offset_policy = __detail::__const_wrapped_accessor_policy<__base_t_>;

  MDSPAN_INLINE_FUNCTION
  static constexpr typename __base_t_::container_type
  create(size_t n) {
    if(n > N) throw length_error("array_container_policy asked for an array larger than static size");
    return typename __base_t_::container_type(n);
  }
};

namespace __detail {

//==============================================================================
// <editor-fold desc="container policy select implementation"> {{{1

template <ptrdiff_t> using __void_ptrdiff_t = void;

// TODO make an alternative to this that works with MSVC?
template <class Map, class=void>
struct __has_constexpr_required_span_size
  : std::false_type
{ static constexpr auto size = 0; };

template <class Map>
struct __has_constexpr_required_span_size<
  Map, __void_ptrdiff_t<(Map{}.required_span_size(), 0)>
> : std::true_type
{ static constexpr auto size = Map{}.required_span_size(); };

template <class, class, class> struct __container_policy_select;

template <class T, class LP, ptrdiff_t... Extents>
struct __container_policy_select<
  T, LP, std::experimental::extents<Extents...>
>
{

  // Avoid even instantiating the constexpr check if there are dynamic
  // extents (it's extra work, and it doesn't mean anything, because
  // the default constructed map doesn't have all of the information
  // to get the right size)
  struct __has_constexpr_req_size_delay {
    template <class Map>
    using __apply = __has_constexpr_required_span_size<Map>;
  };

  struct __false_zero_delay {
    template <class Map>
    struct __apply {
      static constexpr auto value = false;
      static constexpr auto size = 0;
    };
  };

  using __use_array = typename std::conditional<
    _MDSPAN_FOLD_OR(Extents == dynamic_extent),
    __false_zero_delay,
    __has_constexpr_req_size_delay
  >::type::template __apply<
    typename LP::template mapping<std::experimental::extents<Extents...>>
  >;

  // TODO delay instantiation of policy until conditional is resolved?
  using type = typename std::conditional<
    __use_array::value,
    array_container_policy<T, __use_array::size>,
    vector_container_policy<T>
  >::type;
};

// </editor-fold> end container policy select implementation }}}1
//==============================================================================

} // end namespace __detail

} // end inline namespace __mdarray_version_0
} // end namespace experimental
} // end namespace std

#endif //MDARRAY_INCLUDE_EXPERIMENTAL_BITS_CONTAINER_POLICY_BASIC_HPP_
