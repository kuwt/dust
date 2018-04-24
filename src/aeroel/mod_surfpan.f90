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


!> Module to treat surface doublet + source panels 
!!TODO: consider a massive cleanup of unused subroutines

module mod_surfpan

use mod_param, only: &
  wp , & 
  prev_tri , next_tri , &
  prev_qua , next_qua , &
  pi

use mod_aero_elements, only: &
  c_elem, t_elem_p

use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_math, only: &
  cross

!----------------------------------------------------------------------

implicit none

public :: t_surfpan

!----------------------------------------------------------------------

!> Surface panel (Morino kind) with uniform distribution of doublets
!! and sources
!!
!! This type inherits most of its members from the general \ref c_elem class
!! however implements its own subroutines to calculate the coefficients of the
!! linear system and the 
type, extends(c_elem) :: t_surfpan

contains

  procedure, pass(this) :: build_row        => build_row_surfpan
  procedure, pass(this) :: build_row_static => build_row_static_surfpan
  procedure, pass(this) :: add_wake         => add_wake_surfpan
  procedure, pass(this) :: add_liftlin      => add_liftlin_surfpan
  procedure, pass(this) :: compute_pot      => compute_pot_surfpan
  procedure, pass(this) :: compute_vel      => compute_vel_surfpan
  procedure, pass(this) :: compute_psi      => compute_psi_surfpan
  procedure, pass(this) :: compute_cp       => compute_cp_surfpan

end type

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------


!----------------------------------------------------------------------
!> compute AIC of panel <this>, on the control point <pos>
!! Compute I_ik : velocity potential in <pos>_i induced by a (-4*pi)
!!   -intensity surface source  on the panel <this>_k.         ^---------
!!
!! I_ik = int_{S_k} { 1 / |r_i - r| }
!!
!! The relation with unitary surface doublet D_ik is:
!!   S_ik = -(1/(4*pi)) * I_ik
!!
subroutine potential_calc_sou_surfpan(this, sou, dou, pos)
 class(t_surfpan), intent(inout) :: this
 real(wp), intent(out) :: sou
 real(wp), intent(in)  :: dou
 real(wp), intent(in) :: pos(:)

 real(wp) :: radius

 real(wp), dimension(3) :: e3
 real(wp) :: zQ , souLog , vi , R1 , R2
 real(wp), dimension(3) :: Qp 
 integer :: indm1 , indp1

 real(wp), parameter :: eps_sou  = 1.0e-6_wp
 real(wp), parameter :: ff_ratio = 10.0_wp

 integer :: i1

 radius = norm2(pos-this%cen)

 if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then ! far-field approximation (1) 

   sou = this%area / radius

   if ( isnan(sou) ) then
     write(*,*) '  stop. NaN in potential_calc_sou_surfpan. FF '
   end if

 else

   ! initialisation
   sou = 0.0_wp
  
   ! unit normal 
   e3 = this%nor 
  
   ! Control point (Q)
   zQ = sum( (pos-this%cen) * e3 )
   Qp = pos - zQ * e3
  
   do i1 = 1 , this%n_ver
  
     if ( this%n_ver .eq. 3 ) then
       indm1 = prev_tri(i1)
       indp1 = next_tri(i1)
     else if ( this%n_ver .eq. 4 ) then
       indm1 = prev_qua(i1)
       indp1 = next_qua(i1)
     end if
  
     R1 = norm2( pos - this%verp(:,i1) )
     R2 = norm2( pos - this%verp(:,indp1) )
     ! si = this%edge_len(i1)
     souLog = log( (R1+R2+this%edge_len(i1)) / (R1+R2-this%edge_len(i1)) )
 
 
  !  WARNING ---
     if ( abs(R1+R2-this%edge_len(i1)) < 1e-6 ) then
       write(*,*) ' Warning: in potential_calc_sou_surfpan '
       write(*,*) '   abs(R1+R2-this%edge_len(i1)) < 1e-6 '
       write(*,*) ' this%n_ver , this%verp(:,i1) , this%verp(:,indp1) '
       write(*,*) this%n_ver , this%verp(:,i1) , this%verp(:,indp1)
       write(*,*) ' i1 , indp1 , R1 , R2 , this%edge_len(i1) '
       write(*,*) i1 , indp1 , R1 , R2 , this%edge_len(i1)
       stop
     end if
  
     vi = - sum( cross( Qp-this%verp(:,i1), this%edge_vec(:,i1) ) * e3 ) / this%edge_len(i1)
     sou = sou + vi * souLog
  
   end do
  
   sou = sou - zQ * dou

 end if

end subroutine potential_calc_sou_surfpan

!----------------------------------------------------------------------

!> compute velocity AIC of panel <this>, on the control point <pos>
!! Compute the velocity induced in <pos>_i by the surface source <this>_k 
!!  with intentsity intensity -4*pi. <------
!!
!! V^{sou}_ik = grad_{r_i} { int_{S_k} { 1 / |r_i - r| } }
!!
subroutine velocity_calc_sou_surfpan(this, vel, pos)
 class(t_surfpan), intent(inout) :: this
 real(wp), intent(out) :: vel(3)
 real(wp), intent(in) :: pos(:)

 real(wp) :: phix , phiy , pdou
 real(wp) :: R1 , R2 , souLog

 real(wp), parameter :: ff_ratio = 10.0_wp
 real(wp) :: radius_v(3)
 real(wp) :: radius

 integer :: indm1 , indp1
 integer :: i1

 radius_v = pos-this%cen
 radius   = norm2(radius_v)

 if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then ! far-field approximation (1) 

   phix = - this%area * sum( radius_v*this%tang(:,1) ) / radius**3.0_wp 
   phiy = - this%area * sum( radius_v*this%tang(:,2) ) / radius**3.0_wp 
   pdou = - this%area * sum( radius_v*this%nor       ) / radius**3.0_wp 

 else

   phix = 0.0_wp
   phiy = 0.0_wp

   do i1 = 1 , this%n_ver
  
     if ( this%n_ver .eq. 3 ) then
       indm1 = prev_tri(i1)
       indp1 = next_tri(i1)
     else if ( this%n_ver .eq. 4 ) then
       indm1 = prev_qua(i1)
       indp1 = next_qua(i1)
     end if
  
     R1 = norm2( pos - this%verp(:,i1) )
     R2 = norm2( pos - this%verp(:,indp1) )
     ! si = this%edge_len(i1)
     souLog = log( (R1+R2+this%edge_len(i1)) / (R1+R2-this%edge_len(i1)) )
  
     phix = phix + this%sinTi(i1) * souLog
     phiy = phiy - this%cosTi(i1) * souLog
  
   end do
  
   call potential_calc_doublet(this, pdou, pos)

   ! debug : it seems that ONLY the near-field formulas had the wrong sign !!! CHECK it again !!!
   phix = - phix
   phiy = - phiy
   pdou = - pdou

 end if 
 
 ! vsou = (/ phix , phiy , pdou /)
 vel(1) = this%tang(1,1)*phix + this%tang(1,2)*phiy + this%nor(1)* pdou 
 vel(2) = this%tang(2,1)*phix + this%tang(2,2)*phiy + this%nor(2)* pdou 
 vel(3) = this%tang(3,1)*phix + this%tang(3,2)*phiy + this%nor(3)* pdou 

!!!!! ! debug ---- check if the formulas were inverted
!!!!! ! vel = - vel

end subroutine velocity_calc_sou_surfpan

!----------------------------------------------------------------------

!> Build a row of the linear system for a surface panel
!!
!! Only the dynamic part of the linear system is actually built here:
!! the rest of the system was already built in the \ref build_row_static 
!! subroutine. 
subroutine build_row_surfpan(this, elems, linsys, uinf, ie, ista, iend)
 class(t_surfpan), intent(inout) :: this
 type(t_elem_p), intent(in)      :: elems(:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista, iend

 integer :: j1
 real(wp) :: b1(3)
 
  !Components not moving, no body velocity in the boundary condition
  linsys%b(ie) = sum(linsys%b_static(:,ie) * (-uinf))
 


  ! ista and iend will be the end of the unknowns vector, containing
  ! the moving elements
  do j1 = ista , iend
    
    !Compute the AIC and the contribution of j1 to the rhs of ie
    call elems(j1)%p%compute_pot( linsys%A(ie,j1), b1, &
                                  this%cen, ie, j1 )
    
    !Add the contribution to the rhs with the 
    linsys%b(ie) = linsys%b(ie) + sum(b1*(elems(j1)%p%ub-uinf))

  end do

end subroutine build_row_surfpan

!----------------------------------------------------------------------

!> Build a static row of the linear system for a surface panel
!!
!! In this subroutine only the static part of the equations is built. It is
!! called just once at the beginning of the simulation, and saves the AIC 
!! coefficients for te static part and the static contribution to the rhs
subroutine build_row_static_surfpan(this, elems, ll_elems, linsys, uinf, ie, ista, iend)
 class(t_surfpan), intent(inout) :: this
 type(t_elem_p), intent(in)      :: elems(:)
 type(t_elem_p), intent(in)      :: ll_elems(:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista, iend

 integer :: j1
 real(wp) :: b1(3)
 
  linsys%b(ie) = 0.0_wp
  linsys%b_static(:,ie) = 0.0_wp

  !Cycle just all the static elements, ista and iend will be the beginning of 
  !the result vector. Then save the rhs in b_static
  do j1 = ista , iend
 
    call elems(j1)%p%compute_pot( linsys%A(ie,j1), b1,  &
                                  this%cen, ie, j1 )

    linsys%b_static(:,ie) = linsys%b_static(:,ie) + b1
 
  end do

  !Now build the static contribution from the lifting line elements
  do j1 = 1,linsys%nstatic_ll
    call ll_elems(j1)%p%compute_pot( linsys%L_static(ie,j1), b1,  &
                                  this%cen, 1, 2 )
  enddo
  
  
  !The rest of the dynamic part will be completed during the first 
  ! iteration of the assembling

end subroutine build_row_static_surfpan

!----------------------------------------------------------------------

!> Add the contribution of the wake to one equation for a surface panel
!!
!! The rhs of the equation for a surface panel is updated  adding the 
!! the contribution of potential due to the wake
subroutine add_wake_surfpan(this, wake_elems, impl_wake_ind, linsys, uinf, &
                            ie,ista, iend)
 class(t_surfpan), intent(inout) :: this
 type(t_elem_p), intent(in)      :: wake_elems(:)
 integer, intent(in)             :: impl_wake_ind(:,:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend

 integer :: j1, ind1, ind2
 real(wp) :: a, b(3)
 integer :: n_impl
  
  !Count the number of implicit wake contributions
  n_impl = size(impl_wake_ind,2)

  !Add the contribution of the implicit wake panels to the linear system
  do j1 = 1 , n_impl
    ind1 = impl_wake_ind(1,j1); ind2 = impl_wake_ind(2,j1)
    if ((ind1.ge.ista .and. ind1.le.iend) .and. &
        (ind2.ge.ista .and. ind2.le.iend)) then
    
  
      !todo: find a more elegant solution to avoid i=j
      call wake_elems(j1)%p%compute_pot( a, b, this%cen, 1, 2 )
      
      linsys%A(ie,ind1) = linsys%A(ie,ind1) + a
      linsys%A(ie,ind2) = linsys%A(ie,ind2) - a

    endif
    
  end do
  
  ! Add the explicit vortex panel wake contribution to the rhs
  do j1 = n_impl+1 , size(wake_elems)
  
    !todo: find a more elegant solution to avoid i=j
    call wake_elems(j1)%p%compute_pot( a, b, this%cen, 1, 2 )
    
    linsys%b(ie) = linsys%b(ie) - a*wake_elems(j1)%p%idou

  end do

end subroutine add_wake_surfpan

!----------------------------------------------------------------------

!> Add the contribution of the lifing lines to one equation for a surface panel
!!
!! The rhs of the equation for a surface panel is updated  adding the 
!! the contribution of potential due to the lifting lines
subroutine add_liftlin_surfpan(this, ll_elems, linsys, uinf, &
                            ie,ista, iend)
 class(t_surfpan), intent(inout) :: this
 type(t_elem_p), intent(in)      :: ll_elems(:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend

 integer :: j1, ind1, ind2
 real(wp) :: a, b(3)
 integer :: n_impl
  
  !Static part: take what was already computed
  do  j1 = 1, ista-1
    linsys%b(ie) = linsys%b(ie) - linsys%L_static(ie,j1)*ll_elems(j1)%p%idou
  enddo

  !Dynamic part: compute the things now
  do j1 = ista , iend
  
    !todo: find a more elegant solution to avoid i=j
    call ll_elems(j1)%p%compute_pot( a, b, this%cen, 1, 2 )
    
    linsys%b(ie) = linsys%b(ie) - a*ll_elems(j1)%p%idou

  end do

end subroutine add_liftlin_surfpan

!----------------------------------------------------------------------

!> Compute the potential due to a surface panel
!!
!! this subroutine employs doublets and sources basic subroutines to calculate
!! the AIC of a suface panel on another surface panel, and the contribution 
!! to its rhs
subroutine compute_pot_surfpan(this, A, b, pos , i , j )
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(out) :: A
  real(wp), intent(out) :: b(3)
  real(wp), intent(in) :: pos(:)
  integer , intent(in) :: i , j

  real(wp) :: dou , sou

  if ( i .ne. j ) then
    call potential_calc_doublet(this, dou, pos)
  else
!   AIC (doublets) = 0.0   -> dou = 0
    dou = -2.0_wp*pi
  end if

  ! TODO: check coefficients 1/4*pi, ...
  A = -dou


  call potential_calc_sou_surfpan(this, sou, dou, pos)

! b = ... (sources from doublets)
  b =  sou * this%nor

end subroutine compute_pot_surfpan

!----------------------------------------------------------------------

!> Compute the velocity due to a surface panel
!!
!! This subroutine employs doublets and sources basic subroutines to calculate
!! the AIC coefficients of a surface panel to a vortex ring and the 
!! contribution to its rhs
subroutine compute_psi_surfpan(this, A, b, pos, nor, i , j )
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(out) :: A
  real(wp), intent(out) :: b(3)
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: nor(:)
  integer , intent(in) :: i , j

  real(wp) :: vdou(3) , vsou(3)


  call velocity_calc_doublet(this, vdou, pos)

  ! TODO: check coefficients 1/4*pi, ...
  A = sum(vdou * nor)

  call velocity_calc_sou_surfpan(this, vsou, pos)

! b = ... (sources from doublets)
  b =   sum(-vsou * nor ) * this%nor

end subroutine compute_psi_surfpan

!----------------------------------------------------------------------

!> Compute the velocity induced by a surface panel in a prescribed position
!!
!! The velocity in the position is calculated considering the influece of
!! both the doublets and the sources
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of 
!! the equations is multiplied by 4*pi, to obtain the actual velocity the 
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_surfpan(this, pos , uinf, vel )
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: uinf(3)
  real(wp), intent(out) :: vel(3)

  real(wp) :: vdou(3) , vsou(3)

  ! WARNING:
  ! vdou: velocity induced by a doublet of intensity  4*pi ---> v = vdou / 4*pi * idou
  ! vsou: velocity induced by a source  of intensity -4*pi ---> v =-vsou / 4*pi * isou

  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos) 

  ! source ----
  call velocity_calc_sou_surfpan(this, vsou, pos)
 
  !TODO: FIX ALL THIS
  ! In the outer routines vel is used to update the induced velocity
  ! velocity = velocity + vel / ( 4*pi ) 
  ! vel =  vdou * idou - vsou * isou, where isou = psi = this%nor * ( this%b - U_infy )
  !                                   with  this%b the velocity of the collocation point
  ! TODO: pass this%psi and uinf ( up to now uinf = (/-1,0,0/) )
  !vel = vdou*this%idou - vsou*( sum(this%nor*(this%ub+(/-1.0_wp, 0.0_wp, 0.0_wp/))) )
  vel = vdou*this%idou - vsou*( sum(this%nor*(this%ub-uinf)) )

end subroutine compute_vel_surfpan

!----------------------------------------------------------------------

subroutine compute_cp_surfpan(this, elems, uinf) !, uinf)
  class(t_surfpan), intent(inout) :: this
  type(t_elem_p), intent(in) :: elems(:)
  real(wp), intent(in) :: uinf(:)

  real(wp) :: vel_phi(3)

  integer :: i_e

  ! perturbation velocity, u ---------------------------------
  ! Compute velocity from the potential (mu = -phi), exploiting the stencil
  ! contained in pot_vel_stencil
  vel_phi = 0.0_wp
  do i_e = 1 , this%n_ver
    if ( this%i_neigh(i_e) .ne. 0 ) then
      vel_phi = vel_phi + &
        this%pot_vel_stencil(:,i_e) * (elems(this%i_neigh(i_e))%p%idou - this%idou)
    ! else
    end if
  end do

  if (.not. allocated(this%vel)) then !
    allocate(this%vel(3)) ; this%vel = 0.0_wp
    write(*,*) 'allocating this%vel'
  end if

  vel_phi  = - vel_phi    ! mu = - phi

  ! velocity, U = u_t \hat{t} + u_n \hat{n} + U_inf ----------
  this%vel = vel_phi - sum(vel_phi*this%nor)*this%nor    +  &
             this%nor * sum(this%nor * (-uinf+this%ub) ) +  &
             uinf

  ! pressure coefficient, cp ---------------------------------
  ! steady problems  : cp = 1 - (V/V_inf)^2 - dphi/dt * 2/(V_inf^2)
  this%cp  = 1.0_wp - ( norm2(this%vel) / norm2(uinf) )**2.0_wp  &
           + 2.0_wp * this%didou_dt / norm2(uinf)**2.0_wp
 

end subroutine compute_cp_surfpan

!----------------------------------------------------------------------

end module mod_surfpan
