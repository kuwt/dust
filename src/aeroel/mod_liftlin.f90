!!=====================================================================
!!
!! Copyright (C) 2018 Politecnico di Milano
!!
!! This file is part of DUST, an aerodynamic solver for complex
!! configurations.
!! 
!! Permission is hereby granted, free of charge, to any person
!! obtaining a copy of this software and associated documentation
!! files (the "Software"), to deal in the Software without
!! restriction, including without limitation the rights to use,
!! copy, modify, merge, publish, distribute, sublicense, and/or sell
!! copies of the Software, and to permit persons to whom the
!! Software is furnished to do so, subject to the following
!! conditions:
!! 
!! The above copyright notice and this permission notice shall be
!! included in all copies or substantial portions of the Software.
!! 
!! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!! EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
!! OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!! NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
!! HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
!! WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
!! FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
!! OTHER DEALINGS IN THE SOFTWARE.
!! 
!! Authors: 
!!          Federico Fonte             <federico.fonte@polimi.it>
!!          Davide Montagnani       <davide.montagnani@polimi.it>
!!          Matteo Tugnoli             <matteo.tugnoli@polimi.it>
!!=====================================================================

!> Module containing the specific subroutines for the lifting line 
!! type of aerodynamic elements
module mod_liftlin

use mod_param, only: &
  wp, pi

use mod_linsys_vars, only: &
  t_linsys

use mod_aero_elements, only: &
  c_elem, t_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_liftlin


!----------------------------------------------------------------------

type, extends(c_elem) :: t_liftlin
  real(wp), allocatable :: norm_coord_p
contains

  procedure, pass(this) :: build_row        => build_row_liftlin
  procedure, pass(this) :: build_row_static => build_row_static_liftlin
  procedure, pass(this) :: add_wake         => add_wake_liftlin
  procedure, pass(this) :: add_liftlin      => add_liftlin_liftlin
  procedure, pass(this) :: compute_pot      => compute_pot_liftlin
  procedure, pass(this) :: compute_vel      => compute_vel_liftlin
  procedure, pass(this) :: compute_psi      => compute_psi_liftlin
  procedure, pass(this) :: compute_cp       => compute_cp_liftlin
end type

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------
subroutine build_row_liftlin (this, elems, linsys, uinf, ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)       :: elems(:)
 type(t_linsys), intent(inout)    :: linsys
 real(wp), intent(in)             :: uinf(:)
 integer, intent(in)              :: ie
 integer, intent(in)              :: ista, iend

end subroutine build_row_liftlin

!----------------------------------------------------------------------
subroutine build_row_static_liftlin (this, elems, ll_elems, linsys, uinf, ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)       :: elems(:)
 type(t_elem_p), intent(in)       :: ll_elems(:)
 type(t_linsys), intent(inout)    :: linsys
 real(wp), intent(in)             :: uinf(:)
 integer, intent(in)              :: ie
 integer, intent(in)              :: ista, iend

end subroutine build_row_static_liftlin

!----------------------------------------------------------------------
subroutine add_wake_liftlin (this, wake_elems, impl_wake_ind, linsys, uinf, &
                             ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)      :: wake_elems(:)
 integer, intent(in)             :: impl_wake_ind(:,:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend

end subroutine add_wake_liftlin

!----------------------------------------------------------------------
subroutine add_liftlin_liftlin (this, ll_elems, linsys, uinf, &
                             ie, ista, iend)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in)      :: ll_elems(:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend

end subroutine add_liftlin_liftlin

!----------------------------------------------------------------------
subroutine compute_pot_liftlin (this, A, b, pos,i,j)
  class(t_liftlin), intent(inout) :: this
  real(wp), intent(out) :: A
  real(wp), intent(out) :: b(3)
  real(wp), intent(in) :: pos(:)
  integer , intent(in) :: i,j

end subroutine compute_pot_liftlin

!----------------------------------------------------------------------
subroutine compute_psi_liftlin (this, A, b, pos, nor, i, j )
  class(t_liftlin), intent(inout) :: this
  real(wp), intent(out) :: A
  real(wp), intent(out) :: b(3)
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: nor(:)
  integer , intent(in) :: i , j

end subroutine compute_psi_liftlin

!----------------------------------------------------------------------
subroutine compute_vel_liftlin (this, pos, uinf, vel)
  class(t_liftlin), intent(inout) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: uinf(3)
  real(wp), intent(out) :: vel(3)

end subroutine compute_vel_liftlin

!----------------------------------------------------------------------
subroutine compute_cp_liftlin (this, elems, uinf)
  class(t_liftlin), intent(inout) :: this
  type(t_elem_p), intent(in) :: elems(:)
  real(wp), intent(in) :: uinf(:)

end subroutine compute_cp_liftlin 

!----------------------------------------------------------------------

end module mod_liftlin
