
module sll_m_bsplines
  use iso_C_binding, only: c_ptr, c_null_ptr, c_bool, c_double, c_int
  use sll_m_bsplines_types

  implicit none

  public :: sll_bsplines_t, &
      sll_bsplines__new, &
      sll_bsplines__free, &
      sll_bsplines__get_knot,                &
      sll_bsplines__eval_basis,              &
      sll_bsplines__eval_deriv,              &
      sll_bsplines__eval_basis_and_n_derivs, &
      sll_bsplines__integrals

  private
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  interface
      function get_new_bspline_uniform( &
      degree  , &
      periodic, &
      xmin    , &
      xmax    , &
      ncells  )  bind(C)

         use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
         type(c_ptr) :: get_new_bspline_uniform
         integer(c_int) , value , intent(in) :: degree
         logical(c_bool), value , intent(in) :: periodic
         real(c_double) , value , intent(in) :: xmin
         real(c_double) , value , intent(in) :: xmax
         integer(c_int) , value , intent(in) :: ncells
      end function get_new_bspline_uniform
  end interface

  interface
      function get_new_bspline_non_uniform( &
      degree  , &
      periodic, &
      xmin    , &
      xmax    , &
      ncells  , &
      breaks  , &
      nbreaks   )  bind(C)

         use iso_C_binding, only: c_ptr, c_double, c_int, c_bool
         type(c_ptr) :: get_new_bspline_non_uniform
         integer(c_int) , value , intent(in) :: degree
         logical(c_bool), value , intent(in) :: periodic
         real(c_double) , value , intent(in) :: xmin
         real(c_double) , value , intent(in) :: xmax
         integer(c_int) , value , intent(in) :: ncells
         integer(c_int) , value , intent(in) :: nbreaks
         real(c_double)         , intent(in) :: breaks(0:nbreaks-1)
      end function get_new_bspline_non_uniform
  end interface

  interface
      subroutine bsplines_free( bspline  )  bind(C)
         use iso_C_binding, only: c_ptr
         type(c_ptr), value, intent(in   ) :: bspline
      end subroutine bsplines_free
  end interface

  interface
      subroutine bsplines_eval_basis( &
              bspline, &
              x,       &
              nvals,   &
              vals,    &
              jmin) bind(c)
         use iso_C_binding, only: c_ptr, c_int, c_double
         type(c_ptr),    value, intent(in   ) :: bspline
         real(c_double), value, intent(in   ) :: x
         integer(c_int), value, intent(in   ) :: nvals
         real(c_double)       , intent(  out) :: vals(0:nvals-1)
         integer(c_int)       , intent(  out) :: jmin
      end subroutine bsplines_eval_basis
  end interface

  interface
      subroutine bsplines_eval_deriv( &
              bspline, &
              x,       &
              nvals,   &
              vals,    &
              jmin) bind(c)
         use iso_C_binding, only: c_ptr, c_int, c_double
         type(c_ptr),    value, intent(in   ) :: bspline
         real(c_double), value, intent(in   ) :: x
         integer(c_int), value, intent(in   ) :: nvals
         real(c_double)       , intent(  out) :: vals(0:nvals-1)
         integer(c_int)       , intent(  out) :: jmin
      end subroutine bsplines_eval_deriv
  end interface

  interface
      subroutine bsplines_eval_basis_and_n_derivs( &
              bspline, &
              x,       &
              n,       &
              nvals1,  &
              nvals2,  &
              vals,    &
              jmin) bind(c)
         use iso_C_binding, only: c_ptr, c_int, c_double
         type(c_ptr),    value, intent(in   ) :: bspline
         real(c_double), value, intent(in   ) :: x
         integer(c_int), value, intent(in   ) :: n
         integer(c_int), value, intent(in   ) :: nvals1
         integer(c_int), value, intent(in   ) :: nvals2
         real(c_double)       , intent(  out) :: vals(0:nvals1-1,0:nvals2-1)
         integer(c_int)       , intent(  out) :: jmin
      end subroutine bsplines_eval_basis_and_n_derivs
  end interface

  interface
      subroutine bsplines_integrals( &
              bspline, &
              n,       &
              int_vals) bind(C)
         use iso_C_binding, only: c_ptr, c_int, c_double
         type(c_ptr),    value, intent(in   ) :: bspline
         integer(c_int), value, intent(in   ) :: n
         real(c_double)       , intent(  out) :: int_vals(0:n-1)
       end subroutine bsplines_integrals
   end interface

  interface
      function bsplines_get_knot( &
              bspline, &
              idx    ) bind(C) result(knot)
         use iso_C_binding, only: c_ptr, c_int, c_double
         type(c_ptr),    value, intent(in   ) :: bspline
         integer(c_int), value, intent(in   ) :: idx
         real(c_double)                       :: knot
       end function bsplines_get_knot
   end interface

contains

  !-----------------------------------------------------------------------------
  !> @brief      Allocate and initialize uniform or non-uniform B-splines
  !> @param[out] bspline  allocatable polymorphic object
  !> @param[in]  degree    spline degree
  !> @param[in]  periodic  .true. if domain is periodic
  !> @param[in]  xmin      x coordinate of left  boundary of domain
  !> @param[in]  xmax      x coordinate of right boundary of domain
  !> @param[in]  ncells    number of cells in domain (one polynomial per cell)
  !> @param[in]  breaks    list of breakpoints (only for non-uniform B-splines)
  !-----------------------------------------------------------------------------
  function sll_bsplines__new( &
      degree  , &
      periodic, &
      xmin    , &
      xmax    , &
      ncells  , &
      breaks  ) result(bspline)

    type(sll_bsplines_t)                              :: bspline
    integer(c_int)                    , intent(in   ) :: degree
    logical                           , intent(in   ) :: periodic
    real(c_double)                    , intent(in   ) :: xmin
    real(c_double)                    , intent(in   ) :: xmax
    integer(c_int)                    , intent(in   ) :: ncells
    real(c_double),           optional, intent(in   ) :: breaks(:)

    if ( present(breaks) ) then
        bspline % bspline_obj = get_new_bspline_non_uniform( &
            degree, LOGICAL(periodic, c_bool), xmin, &
            xmax, ncells, breaks, size(breaks,1))
    else
        bspline % bspline_obj = get_new_bspline_uniform( &
            degree, LOGICAL(periodic, c_bool), xmin, &
            xmax, ncells)
    endif

    bspline % degree = degree
    bspline % ncells = ncells
    bspline % xmin   = xmin
    bspline % xmax   = xmax
    bspline % nbasis = merge( ncells  , ncells+degree, periodic )
    bspline % periodic = periodic
    bspline % dx     = ( xmax - xmin ) / real(ncells, c_double)

  end function sll_bsplines__new

  !---------------------------------------------------------------------------
  !> @brief        Free storage
  !> @param[inout] self  B-splines object
  !---------------------------------------------------------------------------
  subroutine sll_bsplines__free( self )
   type(sll_bsplines_t), intent(inout) :: self

    call bsplines_free(self % bspline_obj)
    self % bspline_obj = c_null_ptr

  end subroutine sll_bsplines__free

  !-----------------------------------------------------------------------------
  !> @brief     Find the position of the idx-th knot
  !> @param[in] self  uniform B-splines
  !> @param[in] idx   index of interest
  !> results          knot position
  !-----------------------------------------------------------------------------
  function sll_bsplines__get_knot( self, idx ) result(knot)
    type(sll_bsplines_t), intent(in   ) :: self
    integer             , intent(in   ) :: idx
    real(c_double)                      :: knot

    knot = bsplines_get_knot( self % bspline_obj, idx - 1 )

  end function sll_bsplines__get_knot

  !-----------------------------------------------------------------------------
  !> Evaluate value at x of all basis functions with support in local cell
  !> values[j] = B_j(x) for jmin <= j <= jmin+degree
  !>
  !> @param[in]  self    B-splines
  !> @param[in]  x       evaluation point
  !> @param[out] values  array of B-splines' values
  !> @param[out] jmin    index of first non-zero B-spline
  !-----------------------------------------------------------------------------
  subroutine sll_bsplines__eval_basis( self, x, values, jmin )
    type(sll_bsplines_t), intent(in   ) :: self
    real(c_double)      , intent(in   ) :: x
    real(c_double)      , intent(  out) :: values(0:)
    integer(c_int)      , intent(  out) :: jmin

    call bsplines_eval_basis( self % bspline_obj, x, size(values, 1), values, jmin )
    jmin = jmin + 1
  end subroutine sll_bsplines__eval_basis

  !-----------------------------------------------------------------------------
  !> Evaluate derivative at x of all basis functions with support in local cell
  !> derivs[j] = B_j'(x) for jmin <= j <= jmin+degree
  !>
  !> @param[in]  self    B-splines
  !> @param[in]  x       evaluation point
  !> @param[out] derivs  array of B-splines' derivatives
  !> @param[out] jmin    index of first non-zero B-spline
  !-----------------------------------------------------------------------------
  subroutine sll_bsplines__eval_deriv( self, x, derivs, jmin )
    type(sll_bsplines_t), intent(in   ) :: self
    real(c_double)                    , intent(in   ) :: x
    real(c_double)                    , intent(  out) :: derivs(0:)
    integer(c_int)                    , intent(  out) :: jmin

    call bsplines_eval_deriv( self % bspline_obj, x, size(derivs, 1), derivs, jmin )
    jmin = jmin + 1

  end subroutine sll_bsplines__eval_deriv

  !-----------------------------------------------------------------------------
  !> Evaluate value and n derivatives at x of all basis functions with support in local cell
  !> derivs[i,j] = (d/dx)^i B_j(x) for 0 <= i <= n and jmin <= j <= jmin+degree
  !>
  !> @param[in]  self    B-splines
  !> @param[in]  x       evaluation point
  !> @param[in]  n       number of required derivatives
  !> @param[out] derivs  array of B-splines' (multiple) derivatives
  !> @param[out] jmin    index of first non-zero B-spline
  !-----------------------------------------------------------------------------
  subroutine sll_bsplines__eval_basis_and_n_derivs( self, x, n, derivs, jmin )
    type(sll_bsplines_t), intent(in   ) :: self
    real(c_double)      , intent(in   ) :: x
    integer(c_int)      , intent(in   ) :: n
    real(c_double)      , intent(  out) :: derivs(0:,0:)
    integer(c_int)      , intent(  out) :: jmin

    call bsplines_eval_basis_and_n_derivs( self % bspline_obj, x, n, size(derivs, 1), size(derivs, 2), derivs, jmin)
    jmin = jmin + 1

  end subroutine sll_bsplines__eval_basis_and_n_derivs

  subroutine sll_bsplines__integrals( self, int_vals )
     type(sll_bsplines_t), intent(in   ) :: self
     real(c_double)    , intent(  out) :: int_vals(:)

     call bsplines_integrals( self % bspline_obj, size(int_vals,1), int_vals )
  end subroutine

end module sll_m_bsplines
