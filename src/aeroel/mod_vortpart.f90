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

!> Module containing the specific subroutines for a single vortex line
module mod_vortpart

use mod_param, only: &
  wp, pi, max_char_len, prev_tri, next_tri, prev_qua, next_qua

use mod_handling, only: &
  error, printout

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_sim_param, only: &
  t_sim_param

use mod_math, only: &
  cross

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_vortpart, t_vortpart_p, initialize_vortpart


!----------------------------------------------------------------------

type, extends(c_vort_elem) :: t_vortpart
  !> Orientation of the vorticity vector
  real(wp) :: dir(3)
  real(wp) :: vel(3)
  real(wp), pointer :: stretch(:)
  logical :: free=.true.
contains

  procedure, pass(this) :: compute_vel       => compute_vel_vortpart
  procedure, pass(this) :: compute_stretch   => compute_stretch_vortpart  
  procedure, pass(this) :: compute_diffusion => compute_diffusion_vortpart  
  procedure, pass(this) :: calc_geo_data     => calc_geo_data_vortpart

end type

type :: t_vortpart_p

  type(t_vortpart), pointer :: p

end type

character(len=*), parameter :: this_mod_name='mod_vortpart'

real(wp) :: r_Vortex
real(wp) :: r_cutoff

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Initialize vortex line 
subroutine initialize_vortpart(r_Vortex_in, r_cutoff_in)
 real(wp), intent(in) :: r_Vortex_in, r_cutoff_in

  r_Vortex = r_Vortex_in
  r_cutoff  = r_cutoff_in

end subroutine initialize_vortpart

!----------------------------------------------------------------------

!> Compute the velocity induced by a vortex particle in a prescribed position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_vortpart (this, pos, uinf, vel)
 class(t_vortpart), intent(in) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: uinf(3)
 real(wp), intent(out) :: vel(3)

 real(wp) :: vvort(3)
 real(wp) :: dist(3)!, distn


  dist = pos-this%cen
  !Rosenhead kernel regularized velocity
  vvort =  cross(this%dir,dist) / (sqrt(sum(dist**2)+r_Vortex**2))**3
  vel = vvort*this%mag
  
  !Rankine velocity
  !distn = norm2(dist)
  !if ( distn .gt. r_Vortex ) then
  !  vvort =  cross(this%dir,dist) / distn**3
  !else
  !  vvort =  cross(this%dir,dist)  / r_Vortex**3
  !end if
  !vel = vvort*this%mag



end subroutine compute_vel_vortpart

!----------------------------------------------------------------------

!> Compute the vortex stretching induced by a vortex particle 
!! in a prescribed position with a prescribed vorticity (i.e. another particle)
!!
!! WARNING: the calculated term, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_stretch_vortpart (this, pos, alpha, stretch)
 class(t_vortpart), intent(in) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: alpha(3)
 real(wp), intent(out) :: stretch(3)

 real(wp) :: dist(3), distn

  !TODO: add far field approximations

  dist = pos-this%cen
  distn = sqrt(sum(dist**2)+r_Vortex**2)
  !"original"
  !stretch = -cross(alpha, this%dir*this%mag)/(distn)**3 &
  !     +3.0_wp/(distn)**5 * cross(dist, this%mag*this%dir) * &
  !     sum(alpha*dist)
  !"original" fixed sign
  !stretch = -cross(alpha, this%dir*this%mag)/(distn)**3 &
  !     -3.0_wp/(distn)**5 * cross(dist, this%mag*this%dir) * &
  !     sum(alpha*dist)
  !"transpose"
! stretch = -cross(this%dir*this%mag, alpha)/(distn)**3 &
!      +1.0_wp/(distn)**5 * dist * sum(dist*cross(this%dir*this%mag, alpha))

  stretch = -cross(this%dir*this%mag, alpha)/(distn)**3 &
       +1.0_wp/(distn)**5 * dist * sum(dist*cross(this%dir*this%mag, alpha))


end subroutine compute_stretch_vortpart

!----------------------------------------------------------------------

!> Compute the vorticity diffusion induced by a vortex particle 
!! in a prescribed position with a prescribed vorticity (i.e. another particle)
!!
subroutine compute_diffusion_vortpart (this, pos, alpha, diff)
 class(t_vortpart), intent(in) :: this
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: alpha(3)
 real(wp), intent(out) :: diff(3)

 real(wp) :: dist(3), distn
 real(wp) :: volp, volq

  dist = pos-this%cen
  distn = norm2(dist)

  volp = 4.0_wp/3.0_wp*pi*r_Vortex**3
  volq = 4.0_wp/3.0_wp*pi*r_Vortex**3
  diff = 1/(r_Vortex**2)*(volp*this%dir*this%mag - volq*alpha) &
                                                *etaeps(distn,r_Vortex)

end subroutine compute_diffusion_vortpart

!----------------------------------------------------------------------

function etaeps(dist, eps) result(eta)
 real(wp), intent(in) :: dist
 real(wp), intent(in) :: eps
 real(wp) :: eta

  eta = 105.0_wp/(8.0_wp*pi) / ((dist/eps)**2+1)**(9.0_wp/2.0_wp)
  eta = eta/(eps**3)

end function etaeps
!----------------------------------------------------------------------
subroutine calc_geo_data_vortpart(this, vert)
 class(t_vortpart), intent(inout) :: this
 real(wp), intent(in) :: vert(:)

  ! center, it is the only coordinate available
  this%cen = vert

end subroutine calc_geo_data_vortpart

!----------------------------------------------------------------------

end module mod_vortpart
