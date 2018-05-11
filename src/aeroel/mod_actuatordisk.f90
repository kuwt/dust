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

!> Module containing the specific subroutines for the actuator disks
!! elements
module mod_actuatordisk

use mod_param, only: &
  wp, pi

use mod_handling, only: &
  error

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

use mod_aero_elements, only: &
  c_elem, t_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_actdisk


!----------------------------------------------------------------------

type, extends(c_elem) :: t_actdisk
  real(wp), allocatable :: traction
contains

  procedure, pass(this) :: build_row        => build_row_actdisk
  procedure, pass(this) :: build_row_static => build_row_static_actdisk
  procedure, pass(this) :: add_wake         => add_wake_actdisk
  procedure, pass(this) :: add_liftlin      => add_liftlin_actdisk
  procedure, pass(this) :: compute_pot      => compute_pot_actdisk
  procedure, pass(this) :: compute_vel      => compute_vel_actdisk
  procedure, pass(this) :: compute_psi      => compute_psi_actdisk
  procedure, pass(this) :: compute_cp       => compute_cp_actdisk
end type

character(len=*), parameter :: this_mod_name='mod_actuatordisk'

integer :: it=0

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

subroutine build_row_actdisk (this, elems, linsys, uinf, ie, ista, iend)
 class(t_actdisk), intent(inout) :: this
 type(t_elem_p), intent(in)       :: elems(:)
 type(t_linsys), intent(inout)    :: linsys
 real(wp), intent(in)             :: uinf(:)
 integer, intent(in)              :: ie
 integer, intent(in)              :: ista, iend
 character(len=*), parameter      :: this_sub_name='build_row_actdisk'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')
 
end subroutine build_row_actdisk

!----------------------------------------------------------------------

subroutine build_row_static_actdisk (this, elems, ll_elems, linsys, uinf, ie, ista, iend)
 class(t_actdisk), intent(inout) :: this
 type(t_elem_p), intent(in)       :: elems(:)
 type(t_elem_p), intent(in)       :: ll_elems(:)
 type(t_linsys), intent(inout)    :: linsys
 real(wp), intent(in)             :: uinf(:)
 integer, intent(in)              :: ie
 integer, intent(in)              :: ista, iend
 character(len=*), parameter      :: this_sub_name='build_row_static_actdisk'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine build_row_static_actdisk

!----------------------------------------------------------------------

subroutine add_wake_actdisk (this, wake_elems, impl_wake_ind, linsys, uinf, &
                             ie, ista, iend)
 class(t_actdisk), intent(inout) :: this
 type(t_elem_p), intent(in)      :: wake_elems(:)
 integer, intent(in)             :: impl_wake_ind(:,:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend
 character(len=*), parameter      :: this_sub_name='add_wake_actdisk'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine add_wake_actdisk

!----------------------------------------------------------------------

subroutine add_liftlin_actdisk (this, ll_elems, linsys, uinf, &
                             ie, ista, iend)
 class(t_actdisk), intent(inout) :: this
 type(t_elem_p), intent(in)      :: ll_elems(:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend
 character(len=*), parameter      :: this_sub_name='add_liftlin_actdisk'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine add_liftlin_actdisk

!----------------------------------------------------------------------

subroutine compute_pot_actdisk (this, A, b, pos,i,j)
 class(t_actdisk), intent(inout) :: this
 real(wp), intent(out) :: A
 real(wp), intent(out) :: b(3)
 real(wp), intent(in) :: pos(:)
 integer , intent(in) :: i,j

 real(wp) :: dou

  if ( i .ne. j ) then
    !TODO: this is not going to work, need to adapt to n sided
    call potential_calc_doublet(this, dou, pos)
  else
!   AIC (doublets) = 0.0   -> dou = 0
    dou = -2.0_wp*pi
  end if

  A = -dou

  b=0.0_wp

end subroutine compute_pot_actdisk

!----------------------------------------------------------------------

subroutine compute_psi_actdisk (this, A, b, pos, nor, i, j )
 class(t_actdisk), intent(inout) :: this
 real(wp), intent(out) :: A
 real(wp), intent(out) :: b(3)
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: nor(:)
 integer , intent(in) :: i , j

 real(wp) :: vdou(3) 

  call velocity_calc_doublet(this, vdou, pos)

  A = sum(vdou * nor)


  !  b = ... (from boundary conditions)
  !TODO: consider moving this outside
  if ( i .eq. j ) then
    b =  4.0_wp*pi*this%nor
  else
    b = 0.0_wp
  end if

end subroutine compute_psi_actdisk

!----------------------------------------------------------------------
subroutine compute_vel_actdisk (this, pos, uinf, vel)
 class(t_actdisk), intent(inout) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: uinf(3)
 real(wp), intent(out) :: vel(3)

 real(wp) :: vdou(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)
 

  vel = vdou*this%idou


end subroutine compute_vel_actdisk

!----------------------------------------------------------------------
subroutine compute_cp_actdisk (this, elems, uinf)
 class(t_actdisk), intent(inout) :: this
 type(t_elem_p), intent(in) :: elems(:)
 real(wp), intent(in) :: uinf(:)

 character(len=*), parameter      :: this_sub_name='compute_cp_actdisk'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine compute_cp_actdisk 

!----------------------------------------------------------------------

subroutine update_actdisk(elems_ll, linsys)
 type(t_elem_p), intent(inout) :: elems_ll(:)
 type(t_linsys), intent(inout) :: linsys

 real(wp), allocatable :: res_temp(:)


end subroutine update_actdisk

!----------------------------------------------------------------------

subroutine solve_actdisk(elems_ll, elems_tot, elems_wake,  uinf, airfoil_data)
 type(t_elem_p), intent(inout) :: elems_ll(:)
 type(t_elem_p), intent(in)    :: elems_tot(:)
 type(t_elem_p), intent(in)    :: elems_wake(:)
 real(wp), intent(in)          :: uinf(3)
 type(t_aero_tab),  intent(in) :: airfoil_data(:)

end subroutine solve_actdisk

!----------------------------------------------------------------------

end module mod_actuatordisk
