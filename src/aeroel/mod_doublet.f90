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


module mod_doublet

use mod_param, only: &
  wp, &
  prev_tri , next_tri , &
  prev_qua , next_qua , &
  pi

use mod_aero_elements, only: &
  c_elem, c_elem_pan, t_elem_p

use mod_math, only: &
  cross

implicit none

contains

!----------------------------------------------------------------------
!> compute AIC of panel <this>, on the control point <pos>
!! Compute Omega_ik: the solid angle seen from the point <pos>_i (passive)
!!   looking to the panel <this>_k (active). This can be interpreted as the
!!   velocity potential in <pos>_i induced by a (-4*pi)-intensity surface 
!!   doublet on the panel <this>_k.               ^---------
!!
!! Omega_ik = int_{S_k} { n \cdot (r_i - r) / |r_i - r|^3 }
!!
!! The relation with unitary surface doublet D_ik is:
!!   D_ik = -(1/(4*pi)) * Omega_ik
!!
subroutine potential_calc_doublet(this, dou, pos)
 class(c_elem_pan), intent(inout) :: this
 real(wp), intent(out) :: dou
 real(wp), intent(in) :: pos(:)

 real(wp), dimension(3) :: e3
 !real(wp), dimension(3,4) :: ver_p ! Projected vertices on the mean plane

 real(wp) :: zQ
 real(wp), dimension(3) ::Qp

 real(wp), dimension(3) :: ei , ap1 , am1
 real(wp) :: ap1n , am1n , sinB , cosB , beta

 integer :: indm1 , indp1

 real(wp), parameter :: eps_dou  = 1e-6_wp
 real(wp), parameter :: ff_ratio = 10.0_wp
 real(wp) :: radius

 integer :: i1
 
 radius = norm2(pos-this%cen)

 ! unit normal 
 e3 = this%nor 
 
 ! Control point (Q): distance (normal proj) of the point <pos> from the panel <this>
 zQ = sum( (pos-this%cen) * e3 )

 if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then ! far-field approximation (1) 

   dou = zQ * this%area / radius**3.0_wp 

 else

   ! initialisation
   dou = 0.0_wp
  
   Qp = pos - zQ * e3

   do i1 = 1 , this%n_ver
      
     if ( this%n_ver .eq. 3 ) then
       indm1 = prev_tri(i1)
       indp1 = next_tri(i1)
     else if ( this%n_ver .eq. 4 ) then
       indm1 = prev_qua(i1)
       indp1 = next_qua(i1)
     end if
   
     ! doublet  -----
     ! it is possible to use ver, instead of ver_p for the doublet
!    ei = - ( pos - this%verp(:,i1) ) ; ei = ei / norm2(ei)
     ei = - ( pos - this%ver (:,i1) ) ; ei = ei / norm2(ei)
     ap1 =   this%edge_vec(:,i1   ) - ei * sum( ei * this%edge_vec(:,i1   ) )
     am1 = - this%edge_vec(:,indm1) + ei * sum( ei * this%edge_vec(:,indm1) )
     ap1n= norm2(ap1)
     am1n= norm2(am1)
     sinB = sum ( ei * cross(am1,ap1) ) / ( ap1n * am1n )
     cosB = sum ( am1 * ap1 ) / ( ap1n * am1n )
     beta = atan2( sinB , cosB )
     dou = dou + beta
   
   end do

   ! Correct the result to obtain the solid angle (from Gauss-Bonnet theorem)
   if     ( dou .lt. -(this%n_ver-2)*pi - eps_dou ) then
     dou = dou + (this%n_ver-2) * pi
   elseif ( dou .gt. +(this%n_ver-2)*pi + eps_dou ) then
     dou = dou - (this%n_ver-2) * pi 
   else
     dou = 0.0_wp
   end if
   ! Check on the proj of r_i-r normal to the panel (useless?)
!  if ( abs(zQ) .lt. eps_dou ) then
!    dou = 0.0_wp
!  end if

 end if

end subroutine potential_calc_doublet

!----------------------------------------------------------------------
!> compute velocity AIC of panel <this>, on the control point <pos>
!! Biot-Savart law. Regular kernel: linear core (Rankine vortex)
!! Compute the velocity induced by a vortex ring (equivalent to a constant
!!  intentsity surface doublet) with intensity 4*pi. <------
!!
!! V^{vr}_ik = grad_{r_i} { - int_{S_k} { n \cdot (r_i - r) / |r_i - r|^3 } }
!!
subroutine velocity_calc_doublet(this, v_dou, pos)
 class(c_elem_pan), intent(inout) :: this
 real(wp), intent(out) :: v_dou(3)
 real(wp), intent(in) :: pos(:)

 real(wp) :: phix , phiy , pdou
 integer  :: indp1 , indm1
 real(wp) :: av(3) , hv(3)
 real(wp) :: ai    , hi   
 real(wp) :: R1 , R2

 real(wp), parameter :: ff_ratio = 10.0_wp
 real(wp) :: radius_v(3)
 real(wp) :: radius , rati

 real(wp), parameter :: r_Rankine = 0.0000005_wp
 real(wp), parameter :: r_cutoff  = 0.0000001_wp
 real(wp) :: r_Ran

 integer :: i1

 ! WARNING:
 ! to introduce the FLAT-PANEL APPROXIMATION, the projection of the nodes
 ! on the mean palen should be used ...
 ! ti = this%edge_uni(i1)
 ! si = this%edge_len(i1)

 v_dou = 0.0_wp

 radius_v = pos-this%cen
 radius   = norm2(radius_v)

  if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then ! far-field approximation (1) 
    
    rati = 3.0_wp * this%area * sum( radius_v * this%nor ) / radius**5.0_wp
 
    phix =   rati * sum( radius_v*this%tang(:,1) ) 
    phiy =   rati * sum( radius_v*this%tang(:,2) ) 
    pdou =   this%area * ( - radius**2.0_wp + 3.0_wp*sum(radius_v*this%nor)**2.0_wp ) / radius**5.0_wp 
 
    v_dou(1) = this%tang(1,1)*phix + this%tang(1,2)*phiy + this%nor(1)* pdou 
    v_dou(2) = this%tang(2,1)*phix + this%tang(2,2)*phiy + this%nor(2)* pdou 
    v_dou(3) = this%tang(3,1)*phix + this%tang(3,2)*phiy + this%nor(3)* pdou
  
  else

   do i1 = 1 , this%n_ver
   
     if ( this%n_ver .eq. 3 ) then
       indm1 = prev_tri(i1)
       indp1 = next_tri(i1)
     else if ( this%n_ver .eq. 4 ) then
       indm1 = prev_qua(i1)
       indp1 = next_qua(i1)
     end if

     ! use this%ver instead of its projection this%verp   
     !av = pos-this%verp(:,i1)
     av = pos-this%ver(:,i1)
     ai = sum(av*this%edge_uni(:,i1))
     R1 = norm2(av)
     !R2 = norm2(pos-this%verp(:,indp1))
     R2 = norm2(pos-this%ver(:,indp1))
     hv = av - ai*this%edge_uni(:,i1)
     hi = norm2(hv)
     if ( hi .gt. this%edge_len(i1)*r_Rankine ) then
       ! (a/r+(s-a)/r)/h 
       v_dou = v_dou + ( (this%edge_len(i1)-ai)/r2 + ai/r1 )/(hi**2.0_wp) * &
                       cross(this%edge_uni(:,i1),hv)
     else
!      ! (a/r+(s-a)/r)* h/r_Rankine^2.
       if ( ( R1 .gt. this%edge_len(i1)*r_cutoff ) .and. &     ! avoid singularity ...
            ( R2 .gt. this%edge_len(i1)*r_cutoff )   ) then
         r_Ran = r_Rankine * this%edge_len(i1)
         v_dou = v_dou + ( (this%edge_len(i1)-ai)/R2 + ai/R1 )/(r_Ran**2.0_wp) * &
                          cross(this%edge_uni(:,i1),hv)
!      else
       end if
     end if
   
   end do

  end if

end subroutine velocity_calc_doublet

!----------------------------------------------------------------------

end module mod_doublet
