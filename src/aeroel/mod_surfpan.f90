!./\\\\\\\\\\\...../\\\......./\\\..../\\\\\\\\\..../\\\\\\\\\\\\\.
!.\/\\\///////\\\..\/\\\......\/\\\../\\\///////\\\.\//////\\\////..
!..\/\\\.....\//\\\.\/\\\......\/\\\.\//\\\....\///.......\/\\\......
!...\/\\\......\/\\\.\/\\\......\/\\\..\////\\.............\/\\\......
!....\/\\\......\/\\\.\/\\\......\/\\\.....\///\\...........\/\\\......
!.....\/\\\......\/\\\.\/\\\......\/\\\.......\///\\\........\/\\\......
!......\/\\\....../\\\..\//\\\...../\\\../\\\....\//\\\.......\/\\\......
!.......\/\\\\\\\\\\\/....\///\\\\\\\\/..\///\\\\\\\\\/........\/\\\......
!........\///////////........\////////......\/////////..........\///.......
!!=========================================================================
!!
!! Copyright (C) 2018-2020 Davide   Montagnani,
!!                         Matteo   Tugnoli,
!!                         Federico Fonte
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
!!          Federico Fonte             <federico.fonte@outlook.com>
!!          Davide Montagnani       <davide.montagnani@gmail.com>
!!          Matteo Tugnoli                <tugnoli.teo@gmail.com>
!!=========================================================================


!> Module to treat surface doublet + source panels

module mod_surfpan

use mod_param, only: &
  wp , nl, &
  prev_tri , next_tri , &
  prev_qua , next_qua , &
  pi, max_char_len

use mod_handling, only: &
  error, warning, printout

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p


use mod_doublet, only: &
  potential_calc_doublet , &
  velocity_calc_doublet  , &
  gradient_calc_doublet

use mod_linsys_vars, only: &
  t_linsys

use mod_sim_param, only: &
  t_sim_param, sim_param

use mod_math, only: &
  cross , compute_qr

!----------------------------------------------------------------------

implicit none

public :: t_surfpan, initialize_surfpan

private

!----------------------------------------------------------------------

!> Surface panel (Morino kind) with uniform distribution of doublets
!! and sources
!!
!! This type inherits most of its members from the general \ref c_elem class
!! however implements its own subroutines to calculate the coefficients of the
!! linear system and the
type, extends(c_impl_elem) :: t_surfpan

  real(wp), allocatable :: pot_vel_stencil(:,:)
  real(wp), allocatable :: cosTi(:) , sinTi(:)
  real(wp), allocatable :: verp(:,:)
  real(wp)              :: surf_vel(3)
  !> stencil for the Constrained Hermite Taylor Series
  ! Least Square computation of derivatives
  real(wp), allocatable ::   chtls_stencil(:,:)

  !> surface quantities for Bernoulli integral equation
  real(wp) :: dUn_dt
  real(wp) :: bernoulli_source
  real(wp), pointer :: pres_sol

  !> boundary layer and flow separation
  real(wp) :: h_bl      ! height of the surface (boundary??) layer
  real(wp) :: al_free   ! ratio: shed vorticity / total vorticity
  real(wp) :: surf_vort(3) ! free vorticity
  real(wp) :: free_vort(3) ! free vorticity

contains

  procedure, pass(this) :: build_row        => build_row_surfpan
  procedure, pass(this) :: build_row_static => build_row_static_surfpan
  procedure, pass(this) :: add_wake         => add_wake_surfpan
  procedure, pass(this) :: add_expl         => add_expl_surfpan
  procedure, pass(this) :: compute_pot      => compute_pot_surfpan
  procedure, pass(this) :: compute_vel      => compute_vel_surfpan
  procedure, pass(this) :: compute_grad     => compute_grad_surfpan
  procedure, pass(this) :: compute_psi      => compute_psi_surfpan
  procedure, pass(this) :: compute_pres     => compute_pres_surfpan
  procedure, pass(this) :: compute_dforce   => compute_dforce_surfpan
  procedure, pass(this) :: calc_geo_data    => calc_geo_data_surfpan
  procedure, pass(this) :: get_vort_vel     => get_vort_vel_surfpan
  procedure, pass(this) :: correct_pressure_kutta => &
                           correct_pressure_kutta_surfpan

  procedure, pass(this) :: create_local_velocity_stencil => &
                           create_local_velocity_stencil_surfpan
  procedure, pass(this) :: create_chtls_stencil => &
                           create_chtls_stencil_surfpan

  procedure, pass(this) :: get_bernoulli_source => get_bernoulli_source_surfpan

end type

real(wp) :: ff_ratio

character(len=max_char_len) :: msg(3)

character(len=*), parameter :: this_mod_name = 'mod_surfpan'
!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Subroutine to populate the module variables from input
!!
subroutine initialize_surfpan()

  ff_ratio = sim_param%FarFieldRatioSource

end subroutine initialize_surfpan

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
 real(wp) :: den

 real(wp), parameter :: eps_sou  = 1.0e-6_wp
 real(wp), parameter :: ff_ratio = 10.0_wp

 integer :: i1, i2
 character(len=*), parameter :: this_sub_name = 'potential_calc_sou_surfpan'

 radius = norm2(pos-this%cen)

 if ( radius .gt. ff_ratio * maxval(this%edge_len) ) then
   ! far-field approximation (1)

   sou = this%area / radius

   !if ( isnan(sou) ) then
   !workaround to use only standard fortran: by definition a NaN is not equal
   !even to itself
   if ( sou .ne. sou ) then
     call error(this_sub_name, this_mod_name, 'Divide by zero in sources &
      &far field approximation')
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

     den = R1+R2-this%edge_len(i1)
     if ( den < 1e-6_wp ) then
       write(msg(1),'(A)') 'Too small denominator in &
       &source computation with point projection, using actual &
       &points instead.'
       write(msg(2),'(A,F12.6,F12.6,F12.6)') 'Computing sources on point: ',&
       pos(1),pos(2),pos(3)
       write(msg(3),'(A)')'This is most likely due to severely warped &
       &quadrilateral elements adjacent to small elements.'//nl//&
       &'      === CHECK MESH QUALITY! ==='

       call warning(this_sub_name, this_mod_name, msg)
       R1 = norm2( pos - this%ver(:,i1) )
       R2 = norm2( pos - this%ver(:,indp1) )
       den = R1+R2-this%edge_len(i1)
     end if

     souLog = log( (R1+R2+this%edge_len(i1)) / (den) )



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
 class(t_surfpan), intent(in) :: this
 real(wp), intent(out) :: vel(3)
 real(wp), intent(in) :: pos(:)

 real(wp) :: phix , phiy , pdou
 real(wp) :: R1 , R2 , souLog

 real(wp), parameter :: ff_ratio = 10.0_wp
 real(wp) :: radius_v(3)
 real(wp) :: radius

 integer :: indm1 , indp1
 integer :: i1
 character(len=max_char_len) :: message
 character(len=*), parameter :: this_sub_name='velocity_calc_sou_surfpan'

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
     if(sim_param%debug_level .ge.5) then
     if ( abs(R1+R2-this%edge_len(i1)) .lt. 1e-6_wp ) then
      call warning(this_sub_name, this_mod_name, &
        'too small denominator in calculation of velocity')
      write(message,*) ' R1,R2,this%edge_len,i1',R1,R2,this%edge_len(i1),i1
      call printout(message)
     end if
     if ( abs(this%edge_len(i1) ) .lt. 1e-6_wp ) then
      call warning(this_sub_name, this_mod_name, &
        'too small edge length in calculation of velocity')
      write(message,*) ' R1,R2,this%edge_len,i1',R1,R2,this%edge_len(i1),i1
      call printout(message)
     end if
     if ( abs(R1+R2+this%edge_len(i1)) .lt. 1e-6_wp ) then
      call warning(this_sub_name, this_mod_name, &
        'too small numerator in calculation of velocity')
      write(message,*) ' R1,R2,this%edge_len,i1',R1,R2,this%edge_len(i1),i1
      call printout(message)
     end if
     endif

     if ( R1+R2-this%edge_len(i1) .lt. 1e-12_wp ) then
       souLog = 0.0_wp
     else
       souLog = log( (R1+R2+this%edge_len(i1)) / (R1+R2-this%edge_len(i1)) )
     endif


     phix = phix + this%sinTi(i1) * souLog
     phiy = phiy - this%cosTi(i1) * souLog

   end do

   call potential_calc_doublet(this, pdou, pos)

   ! debug : it seems that ONLY the near-field formulas had the wrong sign !!! CHECK it again !!!
   phix = - phix ! * ( -1.0_wp )
   phiy = - phiy ! * ( -1.0_wp )
   pdou = - pdou ! * ( -1.0_wp )

 end if

 ! vsou = (/ phix , phiy , pdou /)
 vel(1) = this%tang(1,1)*phix + this%tang(1,2)*phiy + this%nor(1)* pdou
 vel(2) = this%tang(2,1)*phix + this%tang(2,2)*phiy + this%nor(2)* pdou
 vel(3) = this%tang(3,1)*phix + this%tang(3,2)*phiy + this%nor(3)* pdou

!!!!! ! debug ---- check if the formulas were inverted
!!!!! ! vel = - vel

end subroutine velocity_calc_sou_surfpan

!----------------------------------------------------------------------
!> TODO: compute the gradient of the velocity induced by a source and
!        write this routine
subroutine gradient_calc_sou_surfpan(this, grad, pos)
 class(t_surfpan), intent(in) :: this
 real(wp), intent(out) :: grad(3,3)
 real(wp), intent(in) :: pos(:)

 grad = 0.0_wp

end subroutine gradient_calc_sou_surfpan

!----------------------------------------------------------------------

!> Build a row of the linear system for a surface panel
!!
!! Only the dynamic part of the linear system is actually built here:
!! the rest of the system was already built in the \ref build_row_static
!! subroutine.
subroutine build_row_surfpan(this, elems, linsys, ie, ista, iend)
 class(t_surfpan), intent(inout) :: this
 type(t_impl_elem_p), intent(in)      :: elems(:)
 type(t_linsys), intent(inout)   :: linsys
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista, iend

 integer :: j1 , ipres
 real(wp) :: b1


! RHS for \phi equation --------------------------
  linsys%b(ie) = 0.0_wp
  !Components not moving, no body velocity in the boundary condition
  !linsys%b(ie) = sum(linsys%b_static(:,ie) * (-uinf))
! RHS for Bernoulli polynomial equation ----------
  ipres = linsys%idSurfPanG2L(ie)
  linsys%b_pres(ipres) = 0.0_wp


  ! Static part ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! + \phi equation ------------------------------
  do j1 = 1,ista-1

    linsys%b(ie) = linsys%b(ie) + &
          linsys%b_static(ie,j1) *sum(elems(j1)%p%nor*(-sim_param%u_inf-elems(j1)%p%uvort))
  enddo

  ! + Bernoulli polynomial equation --------------
  do j1 = 1, min(ista-1 , size(linsys%b_static_pres,1))
      linsys%b_pres( ipres ) = linsys%b_pres( ipres ) + &
                               linsys%b_static_pres( ipres , j1 ) * &
                        elems(linsys%idSurfPan(j1))%p%get_bernoulli_source()
  end do

  ! Moving part ++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! ista and iend will be the end of the unknowns vector, containing
  ! the moving elements
  do j1 = ista , iend

    !Compute the AIC and the contribution of j1 to the rhs of ie
    call elems(j1)%p%compute_pot( linsys%A(ie,j1), b1, &
                                  this%cen, ie, j1 )

    ! + \phi equation ----------------------------
    !Add the contribution to the rhs with the
    linsys%b(ie) = linsys%b(ie) &
              + b1* sum(elems(j1)%p%nor*(elems(j1)%p%ub-sim_param%u_inf-elems(j1)%p%uvort))

    ! + Bernoulli polynomial equation ------------
    !Add the contribution to the rhs

    select type( el => elems(j1)%p ) ; class is(t_surfpan)
      linsys%b_pres( ipres ) = &
               linsys%b_pres( ipres ) + &
               b1* el%bernoulli_source

    end select


  end do


end subroutine build_row_surfpan

!----------------------------------------------------------------------

!> Build a static row of the linear system for a surface panel
!!
!! In this subroutine only the static part of the equations is built. It is
!! called just once at the beginning of the simulation, and saves the AIC
!! coefficients for te static part and the static contribution to the rhs
subroutine build_row_static_surfpan(this, elems, expl_elems, linsys, &
                                    ie, ista, iend)
 class(t_surfpan), intent(inout) :: this
 type(t_impl_elem_p), intent(in)      :: elems(:)
 type(t_expl_elem_p), intent(in)      :: expl_elems(:)
 type(t_linsys), intent(inout)   :: linsys
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista, iend

 integer :: j1
 real(wp) :: b1

  linsys%b(ie) = 0.0_wp
  !linsys%b_static(:,ie) = 0.0_wp

  !Cycle just all the static elements, ista and iend will be the beginning of
  !the result vector. Then save the rhs in b_static
  do j1 = ista , iend

    call elems(j1)%p%compute_pot( linsys%A(ie,j1), b1,  &
                                  this%cen, ie, j1 )

    !linsys%b_static(:,ie) = linsys%b_static(:,ie) + b1
    linsys%b_static(ie,j1)  = b1

  end do

  !!Now build the static contribution from the lifting line elements
  !do j1 = 1,linsys%nstatic_ll
  !  call expl_elems(j1)%p%compute_pot( linsys%L_static(ie,j1), b1,  &
  !                                this%cen, 1, 2 )
  !enddo

  !!Now build the static contribution from the actuator disk elements
  !do j1 = 1,linsys%nstatic_ad
  !  call ad_elems(j1)%p%compute_pot( linsys%D_static(ie,j1), b1,  &
  !                                this%cen, 1, 2 )
  !enddo


  !Now build the static contribution from the explicit elements
  do j1 = 1,linsys%nstatic_expl
    call expl_elems(j1)%p%compute_pot( linsys%L_static(ie,j1), b1,  &
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
subroutine add_wake_surfpan(this, wake_elems, impl_wake_ind, linsys, &
                            ie,ista, iend)
 class(t_surfpan), intent(inout) :: this
 type(t_pot_elem_p), intent(in)      :: wake_elems(:)
 integer, intent(in)             :: impl_wake_ind(:,:)
 type(t_linsys), intent(inout)   :: linsys
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend

 integer :: j1, ind1, ind2
 real(wp) :: a, b
 integer :: n_impl

  !Count the number of implicit wake contributions
  n_impl = size(impl_wake_ind,2)

  !Add the contribution of the implicit wake panels to the linear system
  !Implicitly we assume that the first set of wake panels are the implicit
  !ones since are at the beginning of the list
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

    linsys%b(ie) = linsys%b(ie) - a*wake_elems(j1)%p%mag

!   ! debug ---
!   if ( ie .eq. 42 ) then
!     write(*,*) 'wake_elems(' , j1 , '): ' , a
!   end if
!   ! debug ---

  end do

end subroutine add_wake_surfpan

!----------------------------------------------------------------------

! Remove Kutta contribution to obtain the matrix of Bernoulli lin.sys. from
! the matrix of the \phi lin.sys.
! This routine removes the extra terms added in add_wake_surfpan.
! -> OBS: the input IE is already the panel id in the surfpan numeration
subroutine correct_pressure_kutta_surfpan(this, wake_elems, impl_wake_ind, &
                 linsys, uinf, ie,ista, iend)

 class(t_surfpan), intent(inout) :: this
 type(t_pot_elem_p), intent(in)      :: wake_elems(:)
 integer, intent(in)             :: impl_wake_ind(:,:)
 type(t_linsys), intent(inout)   :: linsys
 real(wp), intent(in)            :: uinf(:)
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend

 integer :: j1, ind1, ind2
 real(wp) :: a, b
 integer :: n_impl

  !Count the number of implicit wake contributions
  n_impl = size(impl_wake_ind,2)

  !Remove the contribution of the implicit wake panels to the linear system
  !Implicitly we assume that the first set of wake panels are the implicit
  !ones since are at the beginning of the list
  do j1 = 1 , n_impl
    ind1 = impl_wake_ind(1,j1); ind2 = impl_wake_ind(2,j1)
    if ((ind1.ge.ista .and. ind1.le.iend) .and. &
        (ind2.ge.ista .and. ind2.le.iend)) then

      ! in linsys%idSurfPanG2L, an element is different from zero if it identifies
      ! a surfpan elements. We need to remove only the TE contribution of wake elements
      ! originating from surfpan elements

!     if ( ( linsys%idSurfPanG2L(ind1) .ne. 0 ) .and. &
!          ( linsys%idSurfPanG2L(ind2) .ne. 0 )         ) then
        !todo: find a more elegant solution to avoid i=j
        call wake_elems(j1)%p%compute_pot( a, b, this%cen, 1, 2 )

        linsys%A_pres(ie,linsys%idSurfPanG2L(ind1)) = &
                         linsys%A_pres(ie,linsys%idSurfPanG2L(ind1)) - a
        linsys%A_pres(ie,linsys%idSurfPanG2L(ind2)) = &
                         linsys%A_pres(ie,linsys%idSurfPanG2L(ind2)) + a
!     endif

    endif

  end do

end subroutine correct_pressure_kutta_surfpan

!----------------------------------------------------------------------

!> Add the contribution of the lifing lines to one equation for a surface panel
!!
!! The rhs of the equation for a surface panel is updated  adding the
!! the contribution of potential due to the lifting lines
subroutine add_expl_surfpan(this, expl_elems, linsys, &
                            ie,ista, iend)
 class(t_surfpan), intent(inout) :: this
 type(t_expl_elem_p), intent(in)      :: expl_elems(:)
 type(t_linsys), intent(inout)   :: linsys
 integer, intent(in)             :: ie
 integer, intent(in)             :: ista
 integer, intent(in)             :: iend

 integer :: j1
 real(wp) :: a, b

  !Static part: take what was already computed
  do  j1 = 1, ista-1
    linsys%b(ie) = linsys%b(ie) - linsys%L_static(ie,j1)*expl_elems(j1)%p%mag
  enddo

  !Dynamic part: compute the things now
  do j1 = ista , iend

    !todo: find a more elegant solution to avoid i=j
    call expl_elems(j1)%p%compute_pot( a, b, this%cen, 1, 2 )

    linsys%b(ie) = linsys%b(ie) - a*expl_elems(j1)%p%mag

  end do

end subroutine add_expl_surfpan

!----------------------------------------------------------------------

!> Compute the potential due to a surface panel
!!
!! this subroutine employs doublets and sources basic subroutines to calculate
!! the AIC of a suface panel on another surface panel, and the contribution
!! to its rhs
subroutine compute_pot_surfpan(this, A, b, pos , i , j )
  class(t_surfpan), intent(inout) :: this
  real(wp), intent(out) :: A
  real(wp), intent(out) :: b
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
  !b =  sou * this%nor
  b =  sou

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
  real(wp), intent(out) :: b
  real(wp), intent(in) :: pos(:)
  real(wp), intent(in) :: nor(:)
  integer , intent(in) :: i , j

  real(wp) :: vdou(3) , vsou(3)


  call velocity_calc_doublet(this, vdou, pos)

  ! TODO: check coefficients 1/4*pi, ...
  A = sum(vdou * nor)

  call velocity_calc_sou_surfpan(this, vsou, pos)

! b = ... (sources from doublets)
  b =   sum(-vsou * nor )

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
subroutine compute_vel_surfpan(this, pos, vel )
  class(t_surfpan), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(out) :: vel(3)

  real(wp) :: vdou(3) , vsou(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)

  ! source ----
  call velocity_calc_sou_surfpan(this, vsou, pos)

  vel = vdou*this%mag - vsou*( sum(this%nor*(this%ub-sim_param%u_inf-this%uvort)) )
  ! vel = - vsou*( sum(this%nor*(this%ub-uinf-this%uvort)) )

end subroutine compute_vel_surfpan

!----------------------------------------------------------------------

!> Compute the gradient of the induced velocity in the position pos
!!
!! The velocity in the position is calculated considering the influece of
!! both the doublets and the sources
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_grad_surfpan(this, pos, grad )
  class(t_surfpan), intent(in) :: this
  real(wp), intent(in) :: pos(:)
  real(wp), intent(out) :: grad(3,3)

  real(wp) :: grad_dou(3,3) , grad_sou(3,3)

  ! doublet ---
  call gradient_calc_doublet(this, grad_dou, pos)

  ! source ----
  call gradient_calc_sou_surfpan(this, grad_sou, pos)

  grad = grad_dou*this%mag &
       - grad_sou*( sum(this%nor*(this%ub-sim_param%u_inf-this%uvort)) )

  ! vel = - vsou*( sum(this%nor*(this%ub-uinf-this%uvort)) )

end subroutine compute_grad_surfpan

!----------------------------------------------------------------------

!> Compute an approximate value of the mean pressure on the actual element
! TODO: move the velocity update outside this routine, in a dedicated routine
subroutine compute_pres_surfpan(this, R_g)
  class(t_surfpan) , intent(inout) :: this
  real(wp)         , intent(in)    :: R_g(3,3)
  !type(t_elem_p),   intent(in)    :: elems(:)

  real(wp) :: vel_phi(3), force_pres
  real(wp) :: f(5)    ! <- max n_ver of a surfpan = 4 ; +1 for the constraint eqn

  integer :: i_e , n_neigh
  real(wp) :: mach

! This routine contains the velocity update as well. TODO, move to a dedicated routine
! Two methods have been implemented for surface velocity computation:
! 1. stencil relying on a FV approx of the surface: pot_vel_stencil
! 2. stencil relying on a CHTLS* method
! *CHTLS: Constrained Hermite Taylor series Least Square method

!   ! 1. FV approx
!   ! perturbation velocity, u ---------------------------------
!   ! Compute velocity from the potential (mu = -phi), exploiting the stencil
!   ! contained in pot_vel_stencil.
!   ! ''Tangential component'' from the surface stencil
!   !   Normal component       from the boundary conditions U.n = b.n
!
!   ! tangential part
!   vel_phi_t = 0.0_wp
!   do i_e = 1 , this%n_ver
!     if ( associated(this%neigh(i_e)%p) ) then !  .and. &
!       vel_phi_t = vel_phi_t + &
!         this%pot_vel_stencil(:,i_e) * (this%neigh(i_e)%p%mag - this%mag)
!     end if
!   end do
!
! ! vel_phi_t = - vel_phi_t    ! mu = - phi
!   vel_phi_t = matmul( R_g , vel_phi_t ) ! transpose(R_g) ???
!   vel_phi_t = - vel_phi_t + sum(vel_phi_t*this%nor) * this%nor
!   vel_phi = vel_phi_t +  &
!         sum(this%nor*(this%ub-sim_param%u_inf-this%uvort)) * this%nor
!
!   ! velocity, U = u_t \hat{t} + u_n \hat{n} + U_inf ----------
!   this%surf_vel = sim_param%u_inf + vel_phi + this%uvort
! ! old and wrong
! ! this%surf_vel = vel_phi - sum(vel_phi*this%nor)*this%nor + &
! !      this%nor * sum(this%nor * (-sim_param%u_inf-this%uvort+this%ub) ) + &
! !            sim_param%u_inf + this%uvort

  ! 2. CHTLS method
  n_neigh = 0 ; f = 0.0_wp
  do i_e = 1 , this%n_ver
    if ( associated(this%neigh(i_e)%p) ) then
      n_neigh = n_neigh + 1
      f(n_neigh) = - ( this%neigh(i_e)%p%mag - this%mag )
    end if
  end do
  f(n_neigh+1) = sum(this%nor * (-sim_param%u_inf - this%uvort + this%ub) )

  vel_phi = matmul( this%chtls_stencil , f(1:n_neigh+1) )

  ! Rotation of the result =====================================
  ! - The stencil is computed in the local ref.sys. at the beginning of the code
  ! - The velocity is computed at each timestep
  
  vel_phi = matmul( (R_g) , vel_phi ) 

  !write(*,*) 'vel_phi', vel_phi

  ! Rotation of the result =====================================

  ! vel = u_inf + vel_phi + vel_rot
  this%surf_vel = sim_param%u_inf + vel_phi + this%uvort

  ! pressure -------------------------------------------------
  ! unsteady problems  : P = P_inf + 0.5*rho_inf*V_inf^2
  !                                - 0.5*rho_inf*V^2 - rho_inf*dphi/dt
  !                                     + rho * ub.u_phi
  ! with idou = -phi
!!!  this%pres  = sim_param%P_inf &
!!!! reduced equation after some manipulation
!!!!   - 0.5 * sim_param%rho_inf * norm2(vel_phi+this%uvort)**2.0_wp &
!!!!         - sim_param%rho_inf * sum( &
!!!!            (sim_param%u_inf-this%ub)*(vel_phi+this%uvort) ) &
!!!!         + sim_param%rho_inf * this%didou_dt
!!!! full equation
!!!    + 0.5_wp * sim_param%rho_inf * norm2(sim_param%u_inf)**2.0_wp &
!!!    - 0.5_wp * sim_param%rho_inf * norm2(this%surf_vel)**2.0_wp  &
!!!             + sim_param%rho_inf * sum(this%ub*(vel_phi+this%uvort)) &
!!!             + sim_param%rho_inf * this%didou_dt


  ! === Using result of the Uhlman's equation for pressure ===
  ! See linsys/mod_pressure_equation.f90:assemble_pressure_linsys
  ! as a reference for the 2 different implementations
  !> (a.1) original implementation ===
  ! this%pres = this%pres_sol - 0.5_wp*sim_param%rho_inf * &
  !             norm2(this%surf_vel)**2.0_wp
  !> (a.2) trick of setting B_inf = P_inf + 0.5 * rhoinf * uinf^2.0 ===
  this%pres = this%pres_sol &
            - 0.5_wp*sim_param%rho_inf * norm2(  this%surf_vel)**2.0_wp &
            + 0.5_wp*sim_param%rho_inf * norm2(sim_param%u_inf)**2.0_wp &
            + sim_param%P_inf

  if (this%moving) then
    ! old computation of the loads, expoliting Bernoulli
    ! TODO: to be updated and treated w/ Uhlman's equation formulation
    force_pres  = sim_param%P_inf &
      - 0.5_wp * sim_param%rho_inf * norm2(  this%surf_vel)**2.0_wp  &
      + 0.5_wp * sim_param%rho_inf * norm2(sim_param%u_inf)**2.0_wp &
               + sim_param%rho_inf * sum(this%ub*(vel_phi+this%uvort)) &
               + sim_param%rho_inf * this%didou_dt

    !write(*,*) 'this%surf_vel', this%surf_vel
  else
    force_pres = this%pres 
  endif

! ! OLD PRESSURE EVALUATION
! ! Bernoulli-based pressure for all the components
! force_pres  = sim_param%P_inf &
!   - 0.5_wp * sim_param%rho_inf * norm2(  this%surf_vel)**2.0_wp  &
!   + 0.5_wp * sim_param%rho_inf * norm2(sim_param%u_inf)**2.0_wp &
!            + sim_param%rho_inf * sum(this%ub*(vel_phi+this%uvort)) &
!            + sim_param%rho_inf * this%didou_dt
!
! this%pres = force_pres

  ! Prandt -- Glauert correction for compressibility effect
  mach = abs(norm2(sim_param%u_inf) / sim_param%a_inf)
  
  this%dforce = - (force_pres - sim_param%P_inf) * this%area * this%nor / sqrt(1 - mach**2)

end subroutine compute_pres_surfpan

!----------------------------------------------------------------------

!> Compute the elementary force on the on the actual element
!!
!! WARNING: at the moment this is completely bypassed, due to the
!! different treatment of stationary and moving elements.
!! Its functionalities were moved into compute_pres_surfpan
subroutine compute_dforce_surfpan(this)
  class(t_surfpan), intent(inout) :: this
  !type(t_elem_p), intent(in) :: elems(:)

  ! first rough approximation
  ! vec{F} = - this%pres * vec{n}

  !this%dforce = - this%pres * this%area * this%nor


end subroutine compute_dforce_surfpan

!----------------------------------------------------------------------

!> Compute the coefficients (pot_vel_stencil, for surfpan elements) for
!!  computing the velocity from the velocity potential (phi = -mu).
!! On-body analysis for 3dPanels (surfpan). This coefficients are constant
!!  in the local frame, associated with the component. In order to obtain
!!  the components of the velcoity in the base frame, the global rotation
!!  matrix is needed.
subroutine create_local_velocity_stencil_surfpan ( this , R_g )
  class(t_surfpan)  , intent(inout) :: this
! type(t_pot_elem_p), intent(in)    :: elems(:)
  real(wp)          , intent(in)    :: R_g(3,3)

 real(wp) :: bubble_surf
 integer  :: i_v

 !real(wp) :: n_vect(3)

  if ( .not. allocated(this%pot_vel_stencil) ) then
    allocate(this%pot_vel_stencil(3,this%n_ver) )
  end if

!  ! method #1: very sensitive to dimension gradient of neighbouring elements
!  bubble_surf = this%area
!
!  do i_v = 1 , this%n_ver
!
!    ! Update surf_bubble
!    !sum the contribuition only if the neighbour is really present
!    if(associated(this%neigh(i_v)%p)) then
!
!        bubble_surf = bubble_surf + &
!           this%neigh(i_v)%p%area / real(this%neigh(i_v)%p%n_ver,wp)
!
!    endif
!
!    this%pot_vel_stencil(:,i_v) = &
!             cross( this%edge_vec(:,i_v) , this%nor )
!
!  end do
!
!  this%pot_vel_stencil = this%pot_vel_stencil / bubble_surf

  ! method #2: w/o averaging on neighbouring elements ; 0.5 factor added
  bubble_surf = this%area

  do i_v = 1 , this%n_ver

    this%pot_vel_stencil(:,i_v) = &
     0.5_wp * cross( this%edge_vec(:,i_v) , this%nor )

  end do

!   ! method #3: w/o averaging on neighbouring elements ; 0.5 factor added ;
!   ! taking into account curvature !!!
!   bubble_surf = this%area
!
!   do i_v = 1 , this%n_ver
!
! !   write(*,*) ' this%id ' , this%id
! !   write(*,*) ' this%nor' , this%nor
!     if ( associated(this%neigh(i_v)%p) ) then
! !     write(*,*) ' shape(this%neigh(i_v)) : ' , shape(this%neigh(i_v))
!       n_vect = 0.5_wp * ( this%nor + this%neigh(i_v)%p%nor )
!       n_vect = n_vect / norm2(n_vect)
!     else
!       n_vect = this%nor
!     end if
!     this%pot_vel_stencil(:,i_v) = &
!      0.5_wp * cross( this%edge_vec(:,i_v) , n_vect )
!
!     this%pot_vel_stencil(:,i_v) = &
!                       this%pot_vel_stencil(:,i_v) &
!      - this%nor * sum(this%pot_vel_stencil(:,i_v)*this%nor)
!
!   end do

  this%pot_vel_stencil = this%pot_vel_stencil / bubble_surf

  ! chtls stencil need to be defined in the local ref.sys.
  ! it is built in the global ref.sys. at the first timestep
  !  -> rotation is needed
  this%pot_vel_stencil = matmul( transpose(R_g) , this%pot_vel_stencil )

end subroutine create_local_velocity_stencil_surfpan

!----------------------------------------------------------------------
!> create stencil in the "local" reference system. The components of the
!  velocity vector in the "global" ref.sys. are obtained with the rotation
!  matrix of the "local" reference system
subroutine create_chtls_stencil_surfpan( this , R_g )
  class(t_surfpan)  , intent(inout) :: this
! type(t_pot_elem_p), intent(in)    :: elems(:)
  real(wp)          , intent(in)    :: R_g(3,3)

  real(wp), allocatable :: A(:,:) , B(:,:) , W(:,:) , V(:,:)
  real(wp) :: dx(3)
  real(wp), allocatable :: C(:,:) , CQ(:,:)
  real(wp), allocatable :: Cls_tilde(:,:) , iCls_tilde(:,:) , chtls_tmp(:,:)
  real(wp) :: det_cls
  real(wp) :: r1

  real(wp), allocatable :: Q(:,:) , R(:,:)

  integer :: n_neigh
  integer :: i_n , i_nn

  ! write(*,*) 'R_g_chtls_stencil_surfpan', R_g
  ! # of neighbouring elements
  n_neigh = 0
  do i_n = 1 , this%n_ver
    if ( associated(this%neigh(i_n)%p) ) then ! the neighbouring element exists
      n_neigh = n_neigh + 1
    end if
  end do

  ! arrays A (differences), B (constraints), W (weights) -----------
  allocate( A(n_neigh,3) , W(n_neigh+1,n_neigh+1) ) ; W = 0.0_wp
  allocate( B(1,3) ) ; B(1,:) = this%nor
  i_n = 0
  do i_nn = 1 , this%n_ver
    if ( associated(this%neigh(i_nn)%p) ) then ! the neighbouring element exists
      i_n = i_n + 1
      dx = this%neigh(i_nn)%p%cen - this%cen
      A(i_n, : ) = dx
      W(i_n,i_n) = 1.0_wp / norm2(dx) * &
               max( 0.0_wp , abs( sum( this%nor * this%neigh(i_nn)%p%nor ) ) )
    end if
  end do
  W(n_neigh+1,n_neigh+1) = sum( W ) / real(n_neigh,wp)

  allocate( V(3,1) ) ; V(:,1) = B(1,:)

  call compute_qr( V , Q , R )

  allocate(         C(n_neigh+1,3) ) ; C(1:n_neigh,:) = A ; C(n_neigh+1,:) = B(1,:)
  allocate(        CQ(n_neigh+1,3) ) ; CQ = matmul( C , Q )
  allocate( Cls_tilde(        2,2) ) ; allocate( iCls_tilde(2,2) )
  Cls_tilde = matmul( transpose(CQ(:,2:3)) , matmul( W , CQ(:,2:3) ) )


  ! inverse Cls ----
  det_cls = Cls_tilde(1,1) * Cls_tilde(2,2) - Cls_tilde(1,2) * Cls_tilde(2,1)
  iCls_tilde(1,1) =  Cls_tilde(2,2) / det_cls ; iCls_tilde(1,2) = -Cls_tilde(1,2) / det_cls
  iCls_tilde(2,1) = -Cls_tilde(2,1) / det_cls ; iCls_tilde(2,2) =  Cls_tilde(1,1) / det_cls

  if ( .not. allocated(this%chtls_stencil) ) then
    allocate(this%chtls_stencil( 3 , n_neigh + 1 ) ) ; this%chtls_stencil = 0.0_wp
  end if

  r1 = R(1,1)
  allocate(chtls_tmp( 3 , n_neigh + 1 )) ; chtls_tmp = 0.0_wp
  chtls_tmp(  1,n_neigh + 1 ) = 1.0_wp / R(1,1)
  chtls_tmp(2:3,1 : n_neigh ) = matmul( iCls_tilde , &
        matmul(transpose(CQ(1:n_neigh,2:3)), W(1:n_neigh,1:n_neigh) ) )
  chtls_tmp(2:3, n_neigh+1:n_neigh+1 ) = matmul( iCls_tilde , &
        matmul(transpose(CQ( n_neigh+1:n_neigh+1,2:3)), W(n_neigh+1:n_neigh+1, n_neigh+1:n_neigh+1) ) &
      - matmul( matmul( transpose(CQ(1:n_neigh+1,2:3)), W(1:n_neigh+1,1:n_neigh+1) ) , &
                CQ(1:n_neigh+1,1:1) / r1 )  )

  this%chtls_stencil = matmul( Q , chtls_tmp )

  ! chtls stencil need to be defined in the local ref.sys.
  ! it is built in the global ref.sys. at the first timestep
  !  -> rotation is needed
  this%chtls_stencil = matmul( transpose(R_g) , this%chtls_stencil )


  deallocate(C,CQ,Cls_tilde,iCls_tilde,chtls_tmp)
  deallocate(A,B,V,W,Q,R)


end subroutine create_chtls_stencil_surfpan

!----------------------------------------------------------------------

!> Calculate the geometrical quantities of a surface panel
!!
!! The subroutine calculates all the relevant geometrical quantities of a
!! surface panel
subroutine calc_geo_data_surfpan(this,vert)
 class(t_surfpan), intent(inout) :: this
 real(wp), intent(in) :: vert(:,:)

 integer :: nsides, is
 real(wp):: nor(3), tanl(3)

  this%ver = vert
  nsides = this%n_ver

  ! center
  this%cen =  sum ( this%ver,2 ) / real(nsides,wp)

  ! unit normal and area
  if ( nsides .eq. 4 ) then
    nor = cross( this%ver(:,3) - this%ver(:,1) , &
                 this%ver(:,4) - this%ver(:,2)     )
  else if ( nSides .eq. 3 ) then
    nor = cross( this%ver(:,3) - this%ver(:,2) , &
                 this%ver(:,1) - this%ver(:,2)     )
  end if

  this%area = 0.5_wp * norm2(nor)
  this%nor = nor / norm2(nor)

  ! local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( this%ver(:,nsides) + this%ver(:,1) ) - this%cen

  this%tang(:,1) = tanl / norm2(tanl)
  this%tang(:,2) = cross( this%nor, this%tang(:,1)  )

  ! vector connecting two consecutive vertices:
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  if ( nsides .eq. 3 ) then
    do is = 1 , nsides
      this%edge_vec(:,is) = this%ver(:,next_tri(is)) - this%ver(:,is)
    end do
  else if ( nsides .eq. 4 ) then
    do is = 1 , nsides
      this%edge_vec(:,is) = this%ver(:,next_qua(is)) - this%ver(:,is)
    end do
  end if

  ! edge: edge_len(:)
  do is = 1 , nsides
    this%edge_len(is) = norm2(this%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nSides
    this%edge_uni(:,is) = this%edge_vec(:,is) / this%edge_len(is)
  end do

  !surface panels own fields
  do is = 1 , nsides
    ! cosTi , sinTi
    this%cosTi(is) = sum( this%edge_uni(:,is) * this%tang(:,1) )
    this%sinTi(is) = sum( this%edge_uni(:,is) * this%tang(:,2) )
    ! projection of the vertices on the mean plane
    this%verp(:,is) = this%ver(:,is) - this%nor * &
                    sum( (this%ver(:,is) - this%cen ) * this%nor )
  end do

  !TODO: is it necessary to initialize it here?
  this%dforce = 0.0_wp
  this%dmom   = 0.0_wp


end subroutine calc_geo_data_surfpan

!----------------------------------------------------------------------

!> Calculat the vorticity induced velocity from vortical elements
subroutine get_vort_vel_surfpan(this, vort_elems)
 class(t_surfpan), intent(inout)   :: this
 type(t_vort_elem_p), intent(in)    :: vort_elems(:)

 integer :: iv
 real(wp) :: vel(3)
 
 !> Initialize to zero, before accumulation
 this%uvort = 0.0_wp

 do iv=1,size(vort_elems)
   call vort_elems(iv)%p%compute_vel(this%cen, vel)
   this%uvort = this%uvort + vel/(4*pi)
 enddo

end subroutine

!----------------------------------------------------------------------

function get_bernoulli_source_surfpan(this) result(source)
  class(t_surfpan), intent(inout) :: this
  real(wp) :: source

  source = this%bernoulli_source

end function

!----------------------------------------------------------------------
end module mod_surfpan
