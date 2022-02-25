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

#ifndef MDARRAY_INCLUDE_EXPERIMENTAL___P1684_BITS_CONST_WRAPPED_ACCESSOR_POLICY_HPP_
#define MDARRAY_INCLUDE_EXPERIMENTAL___P1684_BITS_CONST_WRAPPED_ACCESSOR_POLICY_HPP_

#include <experimental/__p0009_bits/macros.hpp>
#include <experimental/__p0009_bits/dynamic_extent.hpp>
#include <experimental/__p0009_bits/extents.hpp>

namespace std {
namespace experimental {
inline namespace __mdarray_version_0 {
namespace __detail {

template <class CP>
struct __const_wrapped_accessor_policy {
 private:

  _MDSPAN_NO_UNIQUE_ADDRESS CP __underlying_cp;

 public:

  using element_type = add_const_t<typename CP::element_type>;
  using pointer = typename CP::const_pointer;
  using reference = typename CP::const_reference;
  // TODO @proposal-bug if we keep offset_policy in P0009, we need `const_offset_policy`...
  using offset_policy = typename CP::const_offset_policy;


  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __const_wrapped_accessor_policy() noexcept(noexcept(CP())) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __const_wrapped_accessor_policy(__const_wrapped_accessor_policy const&)
  noexcept(is_nothrow_copy_constructible<CP>::value) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  constexpr __const_wrapped_accessor_policy(__const_wrapped_accessor_policy&&)
  noexcept(is_nothrow_move_constructible<CP>::value) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  _MDSPAN_CONSTEXPR_14_DEFAULTED __const_wrapped_accessor_policy&
  operator=(__const_wrapped_accessor_policy const&)
  noexcept(is_nothrow_copy_assignable<CP>::value) = default;
  _MDSPAN_CONSTEXPR_14_DEFAULTED __const_wrapped_accessor_policy&
  operator=(__const_wrapped_accessor_policy&&)
  noexcept(is_nothrow_move_assignable<CP>::value) = default;
  MDSPAN_INLINE_FUNCTION_DEFAULTED
  ~__const_wrapped_accessor_policy() noexcept = default;

  MDSPAN_INLINE_FUNCTION
  explicit constexpr
  __const_wrapped_accessor_policy(CP&& underlying)
  noexcept(is_nothrow_move_constructible<CP>::value)
    : __underlying_cp(std::move(underlying))
  { }

  MDSPAN_INLINE_FUNCTION
  explicit constexpr
  __const_wrapped_accessor_policy(CP const& underlying)
  noexcept(is_nothrow_copy_constructible<CP>::value)
    : __underlying_cp(underlying)
  { }

  MDSPAN_INLINE_FUNCTION
  constexpr reference access(pointer ptr, ptrdiff_t i) const
  noexcept(noexcept(__underlying_cp.access(ptr, i)))
  {
    return __underlying_cp.access(ptr, i);
  }

  MDSPAN_INLINE_FUNCTION
  constexpr pointer offset(pointer p, ptrdiff_t i) const
  noexcept(noexcept(__underlying_cp.offset(p, i)))
  {
    return __underlying_cp.offset(p, i);
  }
};

} // end namespace __detail
} // end inline namespace __mdarray_version_0
} // end namespace experimental
} // end namespace std

#endif //MDARRAY_INCLUDE_EXPERIMENTAL___P1684_BITS_CONST_WRAPPED_ACCESSOR_POLICY_HPP_
