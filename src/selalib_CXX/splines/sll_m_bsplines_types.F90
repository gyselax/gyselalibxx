module sll_m_bsplines_types
  use iso_C_binding, only: c_ptr, c_double
  use prec_const, only: F64

  implicit none

  public :: sll_bsplines_t

  private

  type :: sll_bsplines_t
     type(c_ptr) :: bspline_obj
     integer :: degree
     integer :: ncells
     integer :: nbasis
     real(F64) :: dx
     real(F64) :: xmin
     real(F64) :: xmax
     logical :: periodic
  end type
end module sll_m_bsplines_types

