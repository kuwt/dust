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

use mod_handling, only: &
  error

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_aero_elements, only: &
  c_elem, t_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_liftlin, update_liftlin, solve_liftlin


!----------------------------------------------------------------------

type, extends(c_elem) :: t_liftlin
  real(wp)              :: norm_coord_p
  real(wp), allocatable :: tang_cen(:)
  real(wp), allocatable :: bnorm_cen(:)
  real(wp)              :: csi
  integer               :: i_airfoil(2)
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

character(len=*), parameter :: this_mod_name='mod_vortring'

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
 character(len=*), parameter      :: this_sub_name='build_row_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')
 
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
 character(len=*), parameter      :: this_sub_name='build_row_static_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

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
 character(len=*), parameter      :: this_sub_name='add_wake_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

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
 character(len=*), parameter      :: this_sub_name='add_liftlin_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine add_liftlin_liftlin

!----------------------------------------------------------------------

subroutine compute_pot_liftlin (this, A, b, pos,i,j)
 class(t_liftlin), intent(inout) :: this
 real(wp), intent(out) :: A
 real(wp), intent(out) :: b(3)
 real(wp), intent(in) :: pos(:)
 integer , intent(in) :: i,j

 real(wp) :: dou

  if ( i .ne. j ) then
    call potential_calc_doublet(this, dou, pos)
  else
!   AIC (doublets) = 0.0   -> dou = 0
    dou = -2.0_wp*pi
  end if

  A = -dou

  b=0.0_wp

end subroutine compute_pot_liftlin

!----------------------------------------------------------------------

subroutine compute_psi_liftlin (this, A, b, pos, nor, i, j )
 class(t_liftlin), intent(inout) :: this
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

end subroutine compute_psi_liftlin

!----------------------------------------------------------------------
subroutine compute_vel_liftlin (this, pos, uinf, vel)
 class(t_liftlin), intent(inout) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: uinf(3)
 real(wp), intent(out) :: vel(3)

 real(wp) :: vdou(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)

  vel = vdou*this%idou

end subroutine compute_vel_liftlin

!----------------------------------------------------------------------
subroutine compute_cp_liftlin (this, elems, uinf)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in) :: elems(:)
 real(wp), intent(in) :: uinf(:)

 character(len=*), parameter      :: this_sub_name='add_liftlin_liftlin'
 
  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

end subroutine compute_cp_liftlin 

!----------------------------------------------------------------------

subroutine update_liftlin(elems_ll, linsys)
 type(t_elem_p), intent(inout) :: elems_ll(:)
 type(t_linsys), intent(inout) :: linsys

 real(wp), allocatable :: res_temp(:)

  !HERE extrapolate the solution before the linear system
  allocate(res_temp(size(linsys%res_expl,1)))
  res_temp = linsys%res_expl(:,1)
  linsys%res_expl(:,1) = 2.0_wp*res_temp - linsys%res_expl(:,2)
  linsys%res_expl(:,2) = res_temp
  deallocate(res_temp)

end subroutine update_liftlin

!----------------------------------------------------------------------

subroutine solve_liftlin(elems_ll, elems_tot, uinf)
 type(t_elem_p), intent(inout) :: elems_ll(:)
 type(t_elem_p), intent(in)    :: elems_tot(:)
 real(wp), intent(in)          :: uinf(3)

 integer :: i_l, j
 real(wp) :: vel(3), v(3), up(3)
 real(wp) :: unorm, alpha, mach, re
 real(wp) :: aero_coeff(3)

 !Calculate the induced velocity on the airfoil
 do i_l = 1,size(elems_ll)
   vel = 0.0_wp
   do j = 1,size(elems_tot)
     call compute_vel(elems_ll(i_l)%p%cen,uinf,vel)
     vel = vel + v
   enddo
     vel = vel/(4.0_wp*pi) + uinf - elems_ll(i_l)%p%ub
     up = vel-elems_ll(i_l)%p%bnorm_cen*sum(elems_ll(i_l)%p%bnorm_cen*vel)
     unorm = norm2(up)
     alpha = atan2(sum(up*elem_ll(i_l)%p%norm), sum(up*elem_ll(i_l)%p%tang_cen))
     alpha = alpha * 180.0_wp/pi
     !TODO: fix these parameters which are still hard-coded
     mach = 0.0_wp
     re = 1000000.0_wp
     call interp_aero_coeff ( airfoil_data ,  &
                    csi , airfoil_id , (/alpha, mach, re/) , aero_coeff )
 enddo


 !Get the angle of attack, as well as the other parameters

 !Get into the tables to obtain the loads

 !Get the vorticity of the element
end subroutine solve_liftlin

!----------------------------------------------------------------------

end module mod_liftlin
