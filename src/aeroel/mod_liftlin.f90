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
!! Copyright (C) 2018-2019 Davide   Montagnani, 
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

!> Module containing the specific subroutines for the lifting line
!! type of aerodynamic elements
module mod_liftlin

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
  sim_param

use mod_math, only: &
  cross

use mod_c81, only: &
  t_aero_tab, interp_aero_coeff

!use mod_aero_elements, only: &
!  c_elem, t_elem_p

use mod_aeroel, only: &
  c_elem, c_pot_elem, c_vort_elem, c_impl_elem, c_expl_elem, &
  t_elem_p, t_pot_elem_p, t_vort_elem_p, t_impl_elem_p, t_expl_elem_p
!----------------------------------------------------------------------

implicit none

public :: t_liftlin, t_liftlin_p, update_liftlin, solve_liftlin


!----------------------------------------------------------------------

type :: t_liftlin_p
  type(t_liftlin), pointer :: p
end type

type, extends(c_expl_elem) :: t_liftlin
  real(wp), allocatable :: tang_cen(:)
  real(wp), allocatable :: bnorm_cen(:)
  real(wp)              :: csi_cen
  integer               :: i_airfoil(2)
  real(wp)              :: chord
  real(wp)              :: ctr_pt(3)
  real(wp)              :: d_2pi_coslambda
  real(wp)              :: nor_zeroLift(3)
  real(wp)              :: alpha
  real(wp)              :: vel_2d
  real(wp)              :: vel_outplane
  real(wp)              :: aero_coeff(3)
  real(wp)              :: alpha_isolated
  real(wp)              :: vel_2d_isolated
  real(wp)              :: vel_outplane_isolated
contains

  procedure, pass(this) :: compute_pot      => compute_pot_liftlin
  procedure, pass(this) :: compute_vel      => compute_vel_liftlin
  procedure, pass(this) :: compute_psi      => compute_psi_liftlin
  procedure, pass(this) :: compute_pres     => compute_pres_liftlin
  procedure, pass(this) :: compute_dforce   => compute_dforce_liftlin
  procedure, pass(this) :: calc_geo_data    => calc_geo_data_liftlin
end type

character(len=*), parameter :: this_mod_name='mod_liftlin'

integer :: it=0

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

!> Compute the potential due to a lifting line
!!
!! this subroutine employs doublets  to calculate
!! the AIC of a lifting line on a surface panel, adding the contribution
!! to an equation for the potential.
subroutine compute_pot_liftlin (this, A, b, pos,i,j)
 class(t_liftlin), intent(inout) :: this
 real(wp), intent(out) :: A
 real(wp), intent(out) :: b
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

!> Compute the velocity due to a lifting line
!!
!! This subroutine employs doublets basic subroutines to calculate
!! the AIC coefficients of a lifting line to a vortex ring, adding
!! the contribution to an equation for the velocity
subroutine compute_psi_liftlin (this, A, b, pos, nor, i, j )
 class(t_liftlin), intent(inout) :: this
 real(wp), intent(out) :: A
 real(wp), intent(out) :: b
 real(wp), intent(in) :: pos(:)
 real(wp), intent(in) :: nor(:)
 integer , intent(in) :: i , j

 real(wp) :: vdou(3)

  call velocity_calc_doublet(this, vdou, pos)

  A = sum(vdou * nor)


  !  b = ... (from boundary conditions)
  !TODO: consider moving this outside
  if ( i .eq. j ) then
    b =  4.0_wp*pi
  else
    b = 0.0_wp
  end if

end subroutine compute_psi_liftlin

!----------------------------------------------------------------------

!> Compute the velocity induced by a lifting line in a prescribed position
!!
!! WARNING: the velocity calculated, to be consistent with the formulation of
!! the equations is multiplied by 4*pi, to obtain the actual velocity the
!! result of the present subroutine MUST be DIVIDED by 4*pi
subroutine compute_vel_liftlin (this, pos, uinf, vel)
 class(t_liftlin), intent(in) :: this
 real(wp), intent(in) :: pos(:)

 real(wp), intent(in) :: uinf(3)
 real(wp), intent(out) :: vel(3)

 real(wp) :: vdou(3)


  ! doublet ---
  call velocity_calc_doublet(this, vdou, pos)

  vel = vdou*this%mag


end subroutine compute_vel_liftlin

!----------------------------------------------------------------------

subroutine compute_cp_liftlin (this, elems, uinf)
 class(t_liftlin), intent(inout) :: this
 type(t_elem_p), intent(in) :: elems(:)
 real(wp), intent(in) :: uinf(:)

 character(len=*), parameter      :: this_sub_name='compute_cp_liftlin'

  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

!! only steady loads: steady data from table: L -> gam -> p_equiv
!this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
!        norm2(uinf - this%ub) * this%dy / this%area * &
!             elems(this%id)%p%idou

end subroutine compute_cp_liftlin

!----------------------------------------------------------------------

!> The computation of the pressure in the lifting line is not meant to
!! happen, loads are retrieved from the tables
subroutine compute_pres_liftlin (this, R_g)
 class(t_liftlin) , intent(inout) :: this
 real(wp)         , intent(in)    :: R_g(3,3)
 !type(t_elem_p), intent(in) :: elems(:)

 character(len=*), parameter      :: this_sub_name='compute_pres_liftlin'

  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

!! only steady loads: steady data from table: L -> gam -> p_equiv
!this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
!        norm2(uinf - this%ub) * this%dy / this%area * &
!             elems(this%id)%p%idou

end subroutine compute_pres_liftlin

!----------------------------------------------------------------------

!> The computation of loads happens in solve_liftlin and not in
!! this subroutine, which is not used
subroutine compute_dforce_liftlin (this)
 class(t_liftlin), intent(inout) :: this
 !type(t_elem_p), intent(in) :: elems(:)

 character(len=*), parameter      :: this_sub_name='compute_dforce_liftlin'

  call error(this_sub_name, this_mod_name, 'This was not supposed to &
  &happen, a team of professionals is underway to remove the evidence')

!! only steady loads: steady data from table: L -> gam -> p_equiv
!this%cp =   2.0_wp / norm2(uinf)**2.0_wp * &
!        norm2(uinf - this%ub) * this%dy / this%area * &
!             elems(this%id)%p%idou

end subroutine compute_dforce_liftlin

!----------------------------------------------------------------------

!> Update of the solution of the lifting line.
!!
!! Here the extrapolation of the lifting line solution is performed
subroutine update_liftlin(elems_ll, linsys)
 type(t_liftlin_p), intent(inout) :: elems_ll(:)
 !type(t_expl_elem_p), intent(inout) :: elems_ll(:)
 type(t_linsys), intent(inout) :: linsys

 real(wp), allocatable :: res_temp(:)

  it = it + 1

  !HERE extrapolate the solution before the linear system
  if (it .gt. 2) then
    allocate(res_temp(size(linsys%res_expl,1)))
    res_temp = linsys%res_expl(:,1)
    linsys%res_expl(:,1) = 2.0_wp*res_temp - linsys%res_expl(:,2)
    linsys%res_expl(:,2) = res_temp
    deallocate(res_temp)
  else
    linsys%res_expl(:,2) = linsys%res_expl(:,1)
  endif

end subroutine update_liftlin

!----------------------------------------------------------------------

!> Solve the lifting line, in an iterative way
!!
!! The lifting line solution is not obtained from the solution of the linear
!! system. It is fully explicit, but by being nonlinear requires an
!! iterative solution.
subroutine solve_liftlin(elems_ll, elems_tot, &
                         elems_impl, elems_ad, &
                         elems_wake, elems_vort, &
                         airfoil_data, it)
 !type(t_expl_elem_p), intent(inout) :: elems_ll(:)
 type(t_liftlin_p), intent(inout) :: elems_ll(:)
 type(t_pot_elem_p),  intent(in)    :: elems_tot(:)
 type(t_impl_elem_p), intent(in)    :: elems_impl(:)
 type(t_expl_elem_p), intent(in)    :: elems_ad(:)
 type(t_pot_elem_p),  intent(in)    :: elems_wake(:)
 type(t_vort_elem_p), intent(in)    :: elems_vort(:)
 type(t_aero_tab),    intent(in)    :: airfoil_data(:)
 real(wp) :: uinf(3)
 integer  :: i_l, j, ic
 real(wp) :: vel(3), v(3), up(3)
 real(wp), allocatable :: vel_w(:,:)
 real(wp) :: unorm, alpha, alpha_2d
 real(wp) :: cl
 real(wp), allocatable :: aero_coeff(:)
 real(wp), allocatable :: dou_temp(:)
 ! debug ---
 real(wp), allocatable :: al2d_v(:)
 ! debug ---

 ! mach and reynolds number for each el
 real(wp) :: mach , reynolds
 ! arrays used for force projection
 real(wp) , allocatable :: a_v(:)   ! size(elems_ll)
 real(wp) , allocatable :: c_m(:,:) ! size(elems_ll) , 3
 real(wp) , allocatable :: u_v(:)   ! size(elems_ll)

 ! fixed point algorithm for ll ----
 real(wp) :: fp_tol , fp_damp , diff
 integer  :: fp_maxIter
 ! stall regularisation: params read as inputs
 logical  :: stall_regularisation
 real(wp):: al_stall 
 integer :: i_do , i , nn_stall , n_stall
 integer :: n_iter_reg

 real(wp) :: max_mag_ll

 type(t_liftlin), pointer :: el

 integer, intent(in) :: it

 
 character(len=max_char_len) :: msg
 character(len=*), parameter :: this_sub_name = 'solve_liftlin'

  ! params of the fixed point iterations
  fp_tol     = sim_param%llTol
  fp_damp    = sim_param%llDamp
  fp_maxIter = sim_param%llMaxIter
  ! params for stall regularisation
  stall_regularisation = sim_param%llStallRegularisation
  n_stall    = sim_param%llStallRegularisationNelems
  n_iter_reg = sim_param%llStallRegularisationNiters
  al_stall   = sim_param%llStallRegularisationAlphaStall * pi / 180.0_wp

  uinf = sim_param%u_inf

  allocate(dou_temp(size(elems_ll))) ; dou_temp = 0.0_wp

  !=== Compute the velocity from all the elements except for liftling elems ===
  ! and store it outside the loop, since it is constant
  allocate(vel_w(3,size(elems_ll))) ; vel_w = 0.0_wp 
!$omp parallel do private(i_l, j, v) schedule(dynamic)
  do i_l = 1,size(elems_ll)
    do j = 1,size(elems_impl) ! body panels: liftlin, vor che tenga contotlat 
      call elems_impl(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
      vel_w(:,i_l) = vel_w(:,i_l) + v
    enddo
    do j = 1,size(elems_ad) ! actuator disks 
      call elems_ad(j)%p%compute_vel(  elems_ll(i_l)%p%cen,uinf,v)
      vel_w(:,i_l) = vel_w(:,i_l) + v
    enddo
    do j = 1,size(elems_wake) ! wake panels
      call elems_wake(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
      vel_w(:,i_l) = vel_w(:,i_l) + v
    enddo
    do j = 1,size(elems_vort) ! wake vort
      call elems_vort(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
      vel_w(:,i_l) = vel_w(:,i_l) + v
    enddo
    if(sim_param%use_fmm) vel_w(:,i_l) = vel_w(:,i_l) + elems_ll(i_l)%p%uvort*4.0_wp*pi
  enddo
!$omp end parallel do

  vel_w = vel_w/(4.0_wp*pi)

  ! allocate arrasy containing aoa, aero coeffs and relative velocity
  allocate( a_v( size(elems_ll)  )) ;   a_v = 0.0_wp
  allocate( c_m( size(elems_ll),3)) ;   c_m = 0.0_wp
  allocate( u_v( size(elems_ll)  )) ;   u_v = 0.0_wp

  ! debug ---
  allocate( al2d_v( size(elems_ll) ) ) ; al2d_v = 0.0_wp 
  ! debug ---

  ! Remove the "out-of-plane" component of the relative velocity:
  ! 2d-velocity to enter the airfoil look-up-tables
  ! IS THIS LOOP USED (u_v) seems to be overwritten few lines down)
!$omp parallel do private(i_l, el) schedule(dynamic,4)
  do i_l=1,size(elems_ll)
   !select type(el => elems_ll(i_l)%p)
   !type is(t_liftlin)
     el => elems_ll(i_l)%p
     u_v(i_l) = norm2((uinf-el%ub) - &
         el%bnorm_cen*sum(el%bnorm_cen*(uinf-el%ub))) 
     el%vel_2d_isolated = norm2((uinf-el%ub) - &
                          el%bnorm_cen*sum(el%bnorm_cen*(uinf-el%ub)))
     el%vel_outplane_isolated = sum(el%bnorm_cen*(uinf-el%ub))
     el%alpha_isolated = atan2(sum((uinf-el%ub)*el%nor), &
                               sum((uinf-el%ub)*el%tang_cen))*180.0_wp/pi
   !end select
  end do
!$omp end parallel do

  do ic = 1, fp_maxIter
    diff = 0.0_wp             ! max diff ("norm \infty")
    max_mag_ll = 0.0_wp
!$omp parallel do private(i_l, el, j, v, vel, up, unorm, alpha, alpha_2d, mach, &
!$omp& reynolds, aero_coeff, cl) schedule(dynamic,4)
    do i_l = 1,size(elems_ll)

      ! compute velocity
      vel = 0.0_wp
      do j = 1,size(elems_ll)
        call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
        vel = vel + v
      enddo

      !select type(el => elems_ll(i_l)%p) ; type is(t_liftlin)
      el => elems_ll(i_l)%p

        ! overall relative velocity computed in the centre of the ll elem
        vel = vel/(4.0_wp*pi) + uinf - el%ub + vel_w(:,i_l)
        ! "effective" velocity = proj. of vel in the n-t plane
        up =  el%nor*sum(el%nor*vel) + el%tang_cen*sum(el%tang_cen*vel)
        u_v(i_l) = norm2(up) 
        unorm = u_v(i_l)      ! velocity w/o induced velocity
       
        ! Angle of incidence (full velocity)
        alpha = atan2(sum(up*el%nor), sum(up*el%tang_cen))
        alpha = alpha * 180.0_wp/pi  ! .c81 tables defined with angles in [deg]

        ! === Piszkin, Lewinski (1976) LL model for swept wings ===
        ! the control point is approximately at 3/4 of the chord, but the induced
        ! angle of incidence needs to be modified, introducing a "2D correction"
        !
        !> "2D correction" of the induced angle
        alpha_2d = el%mag / ( pi * el%chord * unorm ) *180.0_wp/pi
        alpha = alpha - alpha_2d 

        al2d_v(i_l) = alpha_2d

        ! === Piszkin, Lewinski (1976) LL model for swept wings ===

        ! compute local Reynolds and Mach numbers for the section
        ! needed to enter the LUT (.c81) of aerodynamic loads (2d airfoil)
        mach     = unorm / sim_param%a_inf      
        reynolds = sim_param%rho_inf * unorm * & 
                   el%chord / sim_param%mu_inf    

        ! Read the aero coeff from .c81 tables
        call interp_aero_coeff ( airfoil_data,  el%csi_cen, el%i_airfoil , &
                                   (/alpha, mach, reynolds/) , sim_param , &
                                                              aero_coeff )
        cl = aero_coeff(1)   ! cl needed for the iterative process

        ! Compute the "equivalent" intensity of the vortex line 
        dou_temp(i_l) = - 0.5_wp * unorm * cl * el%chord
!$omp atomic
        diff = max(diff,abs(elems_ll(i_l)%p%mag-dou_temp(i_l)))
!$omp end atomic

        c_m(i_l,:) = aero_coeff
        a_v(i_l)   = alpha * pi/180.0_wp ! [rad]

        el%vel_outplane = sum(el%bnorm_cen*vel)
      !end select

    enddo  ! i_l
!$omp end parallel do


    ! === Stall regularisation ===
    ! avoid "unphysical" stall on an elem between 2 elems without stall. Check
    ! ll section of nn_stall elems, with nn_stall \in (1,llRegularisationNelems)
    ! todo:
    !   - read neighbouring elements and not the prev and next elem in ll numbering
    !   - read al_stall from tables (?), now it is an input from the user 
    if ( stall_regularisation ) then ! *** stall regularisation ***

      if ( ( mod( ic , n_iter_reg ) .eq. 0 ) .or. &
           ( ic .eq. fp_maxIter ) ) then

        al_stall = al_stall * 180.0_wp / pi
        a_v = a_v * 180.0_wp / pi
        i_l = 2 

        nn_stall = 1
        do while( i_l .lt. size(elems_ll) )

          nn_stall = 1 ; i_do = 1

          do while ( ( i_do .eq. 1 ) .and. ( nn_stall .le. n_stall ) .and. &
                    ( i_l + nn_stall .le. size(elems_ll) ) )

            if ( ( all(a_v(i_l:i_l+nn_stall-1) .ge. al_stall ) ) .and. &
                 ( a_v(i_l-1       ) .lt. al_stall )             .and. &
                 ( a_v(i_l+nn_stall) .lt. al_stall ) ) then ! correct

! ! todo: assign a proper debug_level to this screen output
! ! debug ----
!               write(*,*) ' stall_regularisation '
!               write(*,*) ' i_l , i_l+nn_stall-1 : ' , i_l , i_l+nn_stall-1
!               write(*,*) (a_v(i_l:i_l+nn_stall-1) .ge. al_stall ) , '       ' , &
!                          (a_v(i_l-1) .lt. al_stall ) , &
!                          (a_v(i_l+nn_stall) .lt. al_stall)
!               write(*,*) ' a_v(i_l-1,i_l+nn_stall): ' , a_v(i_l-1) , a_v(i_l+nn_stall)
!               write(*,*) ' a_v(i_l:i_l+nn_stall-1)    : ' , a_v(i_l:i_l+nn_stall-1)
! ! debug ---- 
              do i = 1 , nn_stall
                a_v(     i_l+i-1) =  real(i,wp)/real(nn_stall+1,wp) &
                        * a_v(i_l+nn_stall) + (real(nn_stall+1-i,wp))/&
                                       real(nn_stall+1,wp) * a_v(i_l-1)

                dou_temp(i_l+i-1) = real(i,wp)/real(nn_stall+1,wp) * &
                     dou_temp(i_l+nn_stall) + (real(nn_stall+1-i,wp))/&
                                  real(nn_stall+1,wp) * dou_temp(i_l-1)
              end do

! ! todo: assign a proper debug_level to this screen output
! ! debug ---- 
!               write(*,*) ' a_v(i_l:i_l+nn_stall-1) new: ' , a_v(i_l:i_l+nn_stall-1)
! ! debug ---- 

              i_do = 0 ;  ! update variable for the do while loops

            end if

            nn_stall = nn_stall + 1 ! look for larger ll sections

          end do        

          ! update the index of the next elems to analyse for regularisation
          if ( i_do .eq. 0 ) then ;  i_l = i_l + nn_stall-1
          else ;                     i_l = i_l + 1
          end if

        end do ! *** loop over ll elems ***

        a_v = a_v * pi / 180.0_wp
        al_stall = al_stall * pi / 180.0_wp

      end if ! *** if ic/n_iter_reg integer, and ... ***

    end if ! *** stall regularisation ***

    ! === Update ll intensity ===
    do i_l = 1,size(elems_ll)
      elems_ll(i_l)%p%mag = ( dou_temp(i_l)+ fp_damp*elems_ll(i_l)%p%mag )&
                             /(1.0_wp+fp_damp)
      max_mag_ll = max(max_mag_ll,abs(elems_ll(i_l)%p%mag))
    enddo

    if ( diff/max_mag_ll .le. fp_tol ) exit ! convergence

  enddo !solver iterations

! ! todo: assign a proper debug_level to this screen output
  ! debug ----
    write(*,*) '     i_l      ,        ll_alpha  ,          ll_alpha_2d '
    do i_l = 1 , size(elems_ll)
      write(*,*) i_l , a_v(i_l)*180.0_wp/pi , al2d_v(i_l)
    end do
  ! debug ----

! ! debug ---
! write(*,*) ' i_l , u_v , al , cl , cd , chord , area ' 
! ! debug ---

  ! === Loads computation ===
  do i_l = 1,size(elems_ll)

!   ! debug ---
!   write(*,*) i_l , u_v(i_l) , a_v(i_l) *180.0_wp/pi , &
!              c_m(i_l,1) , c_m(i_l,2) , &
!              elems_ll(i_l)%p%chord , &
!              elems_ll(i_l)%p%area
!   ! debug ---

   !select type(el => elems_ll(i_l)%p)
   !type is(t_liftlin)
   el => elems_ll(i_l)%p
    ! avg delta_p = \vec{F}.\vec{n} / A = ( L*cos(al)+D*sin(al) ) / A
    el%pres   = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
               ( c_m(i_l,1) * cos(a_v(i_l)) +  c_m(i_l,2) * sin(a_v(i_l)) )
    ! elementary force = p*n + tangential contribution from L,D
    el%dforce = ( el%nor * el%pres + &
                  el%tang_cen * &
                  0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * ( &
                 -c_m(i_l,1) * sin(a_v(i_l)) + c_m(i_l,2) * cos(a_v(i_l)) &
                 ) ) * el%area
    ! elementary moment = 0.5 * rho * v^2 * A * c * cm, 
    ! - around bnorm_cen (always? TODO: check)
    ! - referred to the ref.point of the elem,
    !   ( here, cen of the elem = cen of the liftlin (for liftlin elems) )
    el%dmom = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
                   el%chord * el%area * c_m(i_l,3)

    el%alpha = a_v(i_l) * 180_wp/pi
    el%vel_2d = u_v(i_l)
    el%aero_coeff = c_m(i_l,:)


   !end select
  end do

  if(sim_param%debug_level .ge. 3) then
    write(msg,*) 'iterations: ',ic; call printout(trim(msg))
    write(msg,*) 'diff',diff; call printout(trim(msg))
    write(msg,*) 'diff/max_mag_ll:',diff/max_mag_ll; call printout(trim(msg))
  endif

  ! useful arrays ---
  deallocate(dou_temp, vel_w)
  deallocate(a_v,c_m,u_v)
  
! ! debug ----
! deallocate(re_v, ma_vm chord_v, cl_v, el_i_airfoil_m, el_csi_cen_v)
! ! debug ----

end subroutine solve_liftlin

!----------------------------------------------------------------------

subroutine calc_geo_data_liftlin(this, vert)
 class(t_liftlin), intent(inout) :: this
 real(wp), intent(in) :: vert(:,:)

 integer  :: is, nsides
 real(wp) :: nor(3), tanl(3) , cen(3)
 real(wp) :: cos_lambda

  this%ver = vert
  nsides = this%n_ver

! ! correct here and at the end of the subroutine
! this%ver(:,3) = this%ver(:,3) + ( this%ver(:,3) - this%ver(:,2) ) / 3.0_wp
! this%ver(:,4) = this%ver(:,4) + ( this%ver(:,4) - this%ver(:,1) ) / 3.0_wp


  ! center, for the lifting line is the mid-point
  this%cen =  sum ( this%ver(:,1:2),2 ) / 2.0_wp
! this%cen =  sum ( this%ver,2 ) / real(nsides,wp) ! <<<< NO ! ............
! ... %cen coinces with control point. %cen must be c/2 far from the ll ...
! ... if the non penetration b.c. is used ( u_rel.n = 0 )               ...

  ! unit normal and area, ll should always have 4 sides
  nor = cross( this%ver(:,3) - this%ver(:,1) , &
               this%ver(:,4) - this%ver(:,2)     )

  ! -- 0.75 chord -- look for other "0.75 chord" tag
  this%area = 0.5_wp * norm2(nor) / 0.75_wp 
  this%nor = nor / norm2(nor)   ! then overwritten

  ! *** 2019-02-27 ***
  ! computation of tangent vectors %tang, moved at the end of the routine...
  ! ...
  ! *** 2019-02-27 ***

  ! vector connecting two consecutive vertices:
  ! edge_vec(:,1) =  ver(:,2) - ver(:,1)
  ! ll should always have 4 sides
  do is = 1 , nsides
    this%edge_vec(:,is) = this%ver(:,next_qua(is)) - this%ver(:,is)
  end do

  ! edge: edge_len(:)
  do is = 1 , nsides
    this%edge_len(is) = norm2(this%edge_vec(:,is))
  end do

  ! unit vector
  do is = 1 , nSides
    this%edge_uni(:,is) = this%edge_vec(:,is) / this%edge_len(is)
  end do

  ! ll-specific fields
  this%tang_cen = this%edge_vec(:,2) - this%edge_vec(:,4)
  this%tang_cen = this%tang_cen / norm2(this%tang_cen)

  this%bnorm_cen = cross(this%tang_cen, this%nor)  ! old
! this%bnorm_cen = this%ver(:,2) - this%ver(:,1)
  this%bnorm_cen = this%bnorm_cen / norm2(this%bnorm_cen)

  ! -- 0.75 chord -- look for other "0.75 chord" tag
  ! correct the chord value ----
  this%chord = sum(this%edge_len((/2,4/)))*0.5_wp
  this%chord = this%chord / 0.75_wp

  ! === Piszkin, Lewinski (1976) LL model for swept wings ===
  ! - cos_lambda = cos(lambda) , where lambda is the local sweep angle
  ! - %cen is overwritten and set = to %cen, because the fmm routines
  !   compute velocity on the %cen of the elems, so far.
  ! - %d_2pi_coslambda = x_CP - x_{1/4*c} * 2*pi * cos(lambda)
  !
  !> sweep angle
  cos_lambda  = norm2( cross( this%tang_cen , -this%edge_uni(:,1) ) )

! ! debug ---
! write(*,*) ' cos_lambda : ' , cos_lambda
! ! debug --- 
 
  !> modified ~3/4 control point
  this%ctr_pt = this%cen + this%tang_cen * this%chord / 2.0_wp 
  ! this%ctr_pt = this%cen + this%tang_cen * this%chord / ( 2.0_wp * cos_lambda )

  !> 2 * pi * | x_CP - x{1/4*c} | * cos(lambda)
  this%d_2pi_coslambda = norm2( this%cen - this%ctr_pt ) * &
                         2.0_wp * pi * cos_lambda
  ! !> overwrite centre
  this%cen    = this%ctr_pt
  !
  ! === Piszkin, Lewinski (1976) LL model for swept wings ===

  ! overwrite nor
  this%nor = cross( this%bnorm_cen , this%tang_cen )
  this%nor = this%nor / norm2(this%nor)

  ! *** 2019-02-27 ***
  ! ...before 2019-02-27 at the beginning of the routine
  ! local tangent unit vector as in PANAIR
  cen =  sum ( this%ver,2 ) / real(nsides,wp)
  tanl = 0.5_wp * ( this%ver(:,nsides) + this%ver(:,1) ) - cen

  ! they should not be used, but ...
  this%tang(:,1) = tanl / norm2(tanl)                   ! this%tang_cen    ! 
  this%tang(:,2) = cross( this%nor, this%tang(:,1)  )   ! this%bnorm_cen   ! 
  ! *** 2019-02-27 ***

  !TODO: is it necessary to initialize it here?
  this%dforce = 0.0_wp
  this%dmom   = 0.0_wp

! ! debug ----
! write(*,*) ' this%cen       :' , this%cen   
! write(*,*) ' this%ctr_pt    :' , this%ctr_pt
! write(*,*) ' this%chord     :' , this%chord 
! write(*,*) ' this%nor       :' , this%nor   
! write(*,*) ' this%tang_cen  :' , this%tang_cen  
! write(*,*) ' this%bnorm_cen :' , this%bnorm_cen 
! write(*,*)
! ! debug ----

! ! correct here and at the beginning of the subroutine
! this%ver(:,3) = this%ver(:,3) - ( this%ver(:,3) - this%ver(:,2) ) / 4.0_wp
! this%ver(:,4) = this%ver(:,4) - ( this%ver(:,4) - this%ver(:,1) ) / 4.0_wp

end subroutine calc_geo_data_liftlin

!----------------------------------------------------------------------

end module mod_liftlin
