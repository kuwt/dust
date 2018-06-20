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

use mod_sim_param, only: &
  t_sim_param

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_actdisk, update_actdisk


!----------------------------------------------------------------------

type, extends(c_expl_elem) :: t_actdisk
  real(wp), allocatable :: traction
  real(wp), allocatable :: radius
contains

  procedure, pass(this) :: compute_pot      => compute_pot_actdisk
  procedure, pass(this) :: compute_vel      => compute_vel_actdisk
  procedure, pass(this) :: compute_psi      => compute_psi_actdisk
  procedure, pass(this) :: compute_pres     => compute_pres_actdisk
  procedure, pass(this) :: compute_dforce   => compute_dforce_actdisk
end type

character(len=*), parameter :: this_mod_name='mod_actuatordisk'

integer :: it=0

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Compute the potential due to an actuator disk
!!
!! this subroutine employs doublets  to calculate
!! the AIC of an actuator disk on a surface panel, adding the contribution
!! to an equation for the potential.
subroutine compute_pot_actdisk (this, A, b, pos,i,j)
 class(t_actdisk), intent(inout) :: this
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

end subroutine compute_pot_actdisk

!----------------------------------------------------------------------

!> Compute the velocity due to an actuator disk
!!
!! This subroutine employs doublets basic subroutines to calculate
!! the AIC coefficients of an actuator disk  to a vortex ring, adding
!! the contribution to an equation for the velocity
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

!> Compute the velocity induced by an actuator disk in a prescribed position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_actdisk (this, pos, uinf, vel)
 class(t_actdisk), intent(inout) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: uinf(3)
 real(wp), intent(out) :: vel(3)

 real(wp) :: vdou(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)


  vel = vdou*this%mag


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

!> The computation of the pressure in the actuator disk is not meant to
!! happen, loads are retrieved from the tables
subroutine compute_pres_actdisk (this, sim_param)
 class(t_actdisk), intent(inout) :: this
 !type(t_elem_p), intent(in) :: elems(:)
 type(t_sim_param), intent(in) :: sim_param

 character(len=*), parameter      :: this_sub_name='compute_pres_actdisk'

  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

!! only steady loads: steady data from table: L -> gam -> p_equiv
!this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
!        norm2(uinf - this%ub) * this%dy / this%area * &
!             elems(this%id)%p%idou

end subroutine compute_pres_actdisk

!----------------------------------------------------------------------

subroutine compute_dforce_actdisk (this, sim_param)
 class(t_actdisk), intent(inout) :: this
 !type(t_elem_p), intent(in) :: elems(:)
 type(t_sim_param), intent(in) :: sim_param

 character(len=*), parameter      :: this_sub_name='compute_dforce_actdisk'

  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

!! only steady loads: steady data from table: L -> gam -> p_equiv
!this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
!        norm2(uinf - this%ub) * this%dy / this%area * &
!             elems(this%id)%p%idou

end subroutine compute_dforce_actdisk

!----------------------------------------------------------------------

subroutine update_actdisk(elems_ad, linsys, sim_param)
 type(t_expl_elem_p), intent(inout) :: elems_ad(:)
 type(t_linsys), intent(inout) :: linsys
 type(t_sim_param), intent(in) :: sim_param

 real(wp), allocatable :: res_temp(:)

 integer :: ie

 do ie=1,size(elems_ad)
   select type(el => elems_ad(ie)%p)
   type is(t_actdisk)

     el%mag = -el%traction*sim_param%dt/(sim_param%rho_inf*pi*el%radius**2)

   end select
 enddo


end subroutine update_actdisk

!----------------------------------------------------------------------

!TODO: is this really necessary?
subroutine solve_actdisk(elems_ll, elems_tot, elems_wake,  uinf, airfoil_data)
 type(t_expl_elem_p), intent(inout) :: elems_ll(:)
 type(t_elem_p), intent(in)    :: elems_tot(:)
 type(t_elem_p), intent(in)    :: elems_wake(:)
 real(wp), intent(in)          :: uinf(3)
 type(t_aero_tab),  intent(in) :: airfoil_data(:)

end subroutine solve_actdisk

!----------------------------------------------------------------------

end module mod_actuatordisk
