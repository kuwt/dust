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

public :: t_liftlin, t_liftlin_p, update_liftlin, solve_liftlin ,  &
          build_ll_array_optim, solve_liftlin_optim , &
          solve_liftlin_optim_regul , &
          solve_liftlin_newton, solve_liftlin_piszkin


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

  !> new field for load computations
  real(wp) :: vel_ctr_pt(3)
  real(wp) ::  al_ctr_pt
  real(wp) :: dn_dt(3)

  !> time derivative of the LL intensity (for unsteady loads)
  real(wp) :: dGamma_dt

contains

  procedure, pass(this) :: compute_pot      => compute_pot_liftlin
  procedure, pass(this) :: compute_vel      => compute_vel_liftlin
  procedure, pass(this) :: compute_psi      => compute_psi_liftlin
  procedure, pass(this) :: compute_pres     => compute_pres_liftlin
  procedure, pass(this) :: compute_dforce   => compute_dforce_liftlin
  procedure, pass(this) :: calc_geo_data    => calc_geo_data_liftlin

  !> new routines for load computations
  procedure, pass(this) :: get_vel_ctr_pt   => get_vel_ctr_pt_liftlin
  procedure, pass(this) :: compute_dforce_jukowski => &
                           compute_dforce_jukowski_liftlin

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
    ! linear extrapolation ---
    linsys%res_expl(:,1) = 2.0_wp*res_temp - linsys%res_expl(:,2)
    linsys%res_expl(:,2) = res_temp
!   ! no extrapolation ---
!   linsys%res_expl(:,1) = res_temp
!   linsys%res_expl(:,2) = res_temp
    deallocate(res_temp)
  else
    linsys%res_expl(:,2) = linsys%res_expl(:,1)
  endif

end subroutine update_liftlin

!----------------------------------------------------------------------
subroutine solve_liftlin_newton( &
                         elems_ll, elems_tot, &
                         elems_impl, elems_ad, &
                         elems_wake, elems_vort, &
                         airfoil_data, it , &
                         M_array )

  type(t_liftlin_p)  , intent(inout) :: elems_ll(:)
  type(t_pot_elem_p) , intent(in)    :: elems_tot(:)
  type(t_impl_elem_p), intent(in)    :: elems_impl(:)
  type(t_expl_elem_p), intent(in)    :: elems_ad(:)
  type(t_pot_elem_p) , intent(in)    :: elems_wake(:)
  type(t_vort_elem_p), intent(in)    :: elems_vort(:)
  type(t_aero_tab)   , intent(in)    :: airfoil_data(:)
  integer            , intent(in)    :: it
  real(wp),allocatable,intent(in)    :: M_array(:,:)

  real(wp) :: artificial_mu = 0.10_wp
  real(wp), allocatable :: mu_v(:)

  real(wp) :: vel(3) , v(3) , up(3)
  real(wp), allocatable :: vel_w(:,:)
  real(wp) :: unorm , alpha , alpha_2d
  !> Gam, Lam, grad_Gam, grad_Lam (for gradient descent)
  real(wp), allocatable :: Gam(:) , Lam(:) , Gamma_old(:) , d_Gam(:)
  real(wp), allocatable :: grad_Gam(:) , grad_Lam(:) , df_dGam(:,:) , f_Gam(:)
  real(wp), allocatable :: a_ik(:,:,:)
  !> uinf
  real(wp) :: uinf(3)
  real(wp) :: mach , reynolds ! mach and reynolds number for each el
  ! arrays used to store aoa, aero coeff and vel
  real(wp) , allocatable :: a_v(:)   ! size(elems_ll)
  real(wp) , allocatable :: c_m(:,:) ! size(elems_ll) , 3
  real(wp) , allocatable :: u_v(:)   ! size(elems_ll)
  real(wp) , allocatable :: dcl_v(:) ! size(elems_ll)
  real(wp) , allocatable :: aero_coeff(:)
  real(wp) :: dcl_da , Jopt
  real(wp) :: diff , max_mag_ll
  !> sim_param: fixed point algorithm for ll
  real(wp) :: fp_tol , fp_damp
  integer  :: fp_maxIter
  !> sim_param: load computation
  logical :: load_avl
  real(wp) :: e_l(3) , e_d(3)

  ! lapack
  integer , allocatable :: ipiv(:)
  integer :: info

  type(t_liftlin), pointer :: el
  integer :: i_l , j , ic

  real(wp) , allocatable :: diff_v(:)
  character(len=4) :: msg_4
  character(len=max_char_len) :: msg
  character(len=*), parameter :: this_sub_name = 'solve_liftlin_newton'
 
  allocate( diff_v(sim_param%llMaxIter+1) ) ; diff_v = 0.0_wp

  ! params of the fixed point iterations
  fp_tol     = sim_param%llTol
  fp_damp    = sim_param%llDamp
  fp_maxIter = sim_param%llMaxIter
  ! param for load computation
  load_avl   = sim_param%llLoadsAVL
 
  uinf = sim_param%u_inf

  !> allocate and fill Gamma_old array of the ll intensity at previous dt
  allocate(Gamma_old(size(elems_ll)))
  do i_l = 1 , size(elems_ll)
    Gamma_old(i_l) = elems_ll(i_l)%p%mag
  end do

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

  ! allocate array containing aoa, aero coeffs and relative velocity
  allocate(   a_v( size(elems_ll)  )) ;   a_v = 0.0_wp
  allocate(   c_m( size(elems_ll),3)) ;   c_m = 0.0_wp
  allocate(   u_v( size(elems_ll)  )) ;   u_v = 0.0_wp
  allocate( dcl_v( size(elems_ll)  )) ; dcl_v = 0.0_wp

  ! Remove the "out-of-plane" component of the relative velocity:
  ! 2d-velocity to enter the airfoil look-up-tables
  ! IS THIS LOOP USED? (u_v) seems to be overwritten few lines down)
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

  ! === Allocate and initialize Gam and Lam arrays ===
  allocate(   d_Gam(size(elems_ll))) ;    d_Gam = 0.0_wp
  allocate(     Gam(size(elems_ll))) ;      Gam = Gamma_old 
  allocate(     Lam(size(elems_ll))) ;      Lam = 1.0_wp
  allocate(grad_Gam(size(elems_ll))) ; grad_Gam = 0.0_wp
  allocate(grad_Lam(size(elems_ll))) ; grad_Lam = 0.0_wp
  allocate(   f_Gam(size(elems_ll))) ;    f_Gam = 0.0_wp
  allocate( df_dGam(size(elems_ll),size(elems_ll))) ;  df_dGam = 0.0_wp
  allocate(    mu_v(size(elems_ll))) ;     mu_v = 0.0_wp

  ! === Save AIC for derivatives d(alpha)/d(gamma) ===
  allocate( a_ik( size(elems_ll) , size(elems_ll) , 3 ) ) ; a_ik = 0.0_wp
  !> set intensity of the ll elems to 1, to compute AICs
  do i_l = 1 , size(elems_ll)
    elems_ll(i_l)%p%mag = 1.0_wp
  end do 
  !> compute AICs for derivatives
  do i_l = 1 , size(elems_ll)
    do j = 1 , size(elems_ll)
      call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf, a_ik(i_l,j,:) )
      a_ik(i_l,j,:) = a_ik(i_l,j,:) / ( 4.0_wp * pi )
    end do
  end do
  !> set intensity of the ll elems to their actual value
  do i_l = 1 , size(elems_ll)
    elems_ll(i_l)%p%mag = Gamma_old(i_l)
  end do 

  ! === Iterative algorithm ===
  do ic = 1 , fp_maxIter

    !> Reset diff and max_mag_ll variables
    diff = 0.0_wp ; max_mag_ll = 0.0_wp

!$omp parallel do private(i_l, el, j, v, vel, up, unorm, alpha, alpha_2d, mach, &
!$omp& reynolds, aero_coeff, dcl_da) schedule(dynamic,4)
    do i_l = 1 , size(elems_ll)
 
      ! === Compute velocity ===
      !> LL influence ( to be divided by 4*pi )
      vel = 0.0_wp
      do j = 1 , size(elems_ll)
        call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
        vel = vel + v
      enddo

      el => elems_ll(i_l)%p

      !> overall relative velocity computed in the centre of the ll elem
      vel = vel/(4.0_wp*pi) + uinf - el%ub + vel_w(:,i_l)
      !> "effective" velocity = proj. of vel in the n-t plane
      up =  el%nor*sum(el%nor*vel) + el%tang_cen*sum(el%tang_cen*vel)
      unorm = norm2(up)      ! velocity w/o induced velocity
      
      ! === Angle of incidence (full velocity) ===
      alpha = atan2(sum(up*el%nor), sum(up*el%tang_cen))
      alpha = alpha * 180.0_wp/pi  ! .c81 tables defined with angles in [deg]

      ! === Piszkin, Lewinski (1976) LL model for swept wings ===
      ! the control point is approximately at 3/4 of the chord, but the induced
      ! angle of incidence needs to be modified, introducing a "2D correction"
      !
      !> "2D correction" of the induced angle
      alpha_2d = el%mag / ( pi * el%chord * unorm ) *180.0_wp/pi
      alpha = alpha - alpha_2d 
      ! =========================================================

      ! === Mach and Reynolds numbers ===
      mach     = unorm / sim_param%a_inf      
      reynolds = sim_param%rho_inf * unorm * el%chord / sim_param%mu_inf    

      ! === Read the aero coeff from .c81 tables ===
      call interp_aero_coeff ( airfoil_data,  el%csi_cen, el%i_airfoil , &
                                 (/alpha, mach, reynolds/) , sim_param , &
                                                  aero_coeff ,  dcl_da )

      ! === Fill arrays for iterative method ===
      u_v(  i_l)   = unorm
      a_v(  i_l)   = alpha * pi / 180.0_wp
      c_m(  i_l,:) = aero_coeff
      dcl_v(i_l)   = dcl_da

    end do
!$omp end parallel do

!   ! debug ---
!   if ( ic .eq. 1 ) then
!     do i_l = 1 , size(elems_ll)
!       write(*,*) ' alpha, dcl_da(' , i_l ,') : ' , a_v(i_l) * 180.0_wp/pi , dcl_v(i_l)
!     end do
!   end if
!   ! debug ---

    !> Compute nonlinear function F(Gam)
    f_Gam   = 0.0_wp ! compute_f_lam()
    do i_l = 1 , size(elems_ll)
      f_Gam(i_l) = 0.5_wp * u_v(i_l) * c_m(i_l,1) * elems_ll(i_l)%p%chord + Gam(i_l) 
    end do
    
    !> Stopping criterion
    diff_v(ic) = maxval( abs( f_Gam - mu_v * matmul( M_array , Gam ) ) )
    if ( diff_v(ic) .lt. fp_tol ) exit ! convergence

    !> Jacobian
    df_dGam = 0.0_wp ! compute_df_dlam()
    do i_l = 1 , size(elems_ll)
      do j = 1 , size(elems_ll)
        df_dGam(i_l,j) = 0.5_wp * elems_ll(i_l)%p%chord * &
                         dcl_v(i_l) * &
                  sum( a_ik(i_l,j,:) * & 
                       ( elems_ll(i_l)%p%nor      * cos(a_v(i_l)) &
                       - elems_ll(i_l)%p%tang_cen * sin(a_v(i_l)) ) )
      end do
      df_dGam(i_l,i_l) = df_dGam(i_l,i_l) + 1.0_wp
    end do

    ! === Artificial viscosity regularisation ===
    do i_l = 1 , size(elems_ll)
      !> artificial viscosity following Chattot / Gallay&Laurendeau
      mu_v(i_l) = - min( 0.25_wp * elems_ll(i_l)%p%chord * &
                         dcl_v(i_l) * &
                         sum( a_ik(i_l,i_l,:) * & 
                         ( elems_ll(i_l)%p%nor      * cos(a_v(i_l)) &
                         - elems_ll(i_l)%p%tang_cen * sin(a_v(i_l)) ) ) &
                       , 0.0_wp ) * 1.5_wp

!     ! debug ---
!     write(*,*) ' a_v , mu_v(',i_l,') : ' , a_v(i_l)*180_wp/pi , mu_v(i_l)
!     ! debug ---

      !> Update Jacobian dF/dGam
      df_dGam(i_l,:) = df_dGam(i_l,:) &
                     - mu_v(i_l) * M_array(i_l,:)
      !> Update Nonlinear function F(Gam)
      f_Gam(i_l)     = f_Gam(i_l) &
                     - mu_v(i_l) * sum( M_array(i_l,:) * Gam )

    end do

    !> Solve Newton linear system
    d_Gam = f_Gam        ! lapack overwrites the solution to the rhs

    allocate(ipiv(size(elems_ll)))
    call dgesv( size(elems_ll) , 1 , df_dGam , size(elems_ll) , &
                          ipiv ,       d_Gam , size(elems_ll) , info )
    deallocate(ipiv)

    !> Update ll intensity
    do i_l = 1 , size(elems_ll)

      Gam(i_l) = Gam(i_l) - 0.001_wp * d_Gam(i_l)
      elems_ll(i_l)%p%mag = Gam(i_l)

    enddo

    !> Stopping criterion
    if ( ic .ge. fp_maxIter ) exit

  end do

  ! debug ---
  write(*,*) '         i_l ,     artif_visc ,              alpha[deg] ,                  cl ,  &
  &                   dcl_da ,                   Gam ,               &
  &     f_gam ,                  Gam+f_Gam ,              grad_Gam'
  do i_l = 1 , size(elems_ll)
    write(*,*) i_l , mu_v(i_l), a_v(i_l)*180.0_wp/pi , c_m(i_l,1) , dcl_v(i_l) , &
           Gam(i_l) , f_gam(i_l) , Gam(i_l)+f_Gam(i_l) , grad_Gam(i_l)
  end do
  ! debug ---

  write(msg_4,'(I4.4)') it
  open(unit=21,file=trim(sim_param%basename_debug)//'_conv_'//msg_4//'.dat')
  do i_l = 1 , ic-1
    write(21,*) diff_v(i_l)
  end do
  close(21)

  ! === Update dGamma_dt field ===
  do i_l = 1,size(elems_ll)
    elems_ll(i_l)%p%dGamma_dt = ( elems_ll(i_l)%p%mag - Gamma_old(i_l) ) / &
                                  sim_param%dt 
  end do

  ! === Load computations === ( same as solve_liftlin() routine )
  ! compute LL sectional loads from singularity intensities. 
  do i_l = 1,size(elems_ll)

   !select type(el => elems_ll(i_l)%p)
   !type is(t_liftlin)
   el => elems_ll(i_l)%p
    ! avg delta_p = \vec{F}.\vec{n} / A = ( L*cos(al)+D*sin(al) ) / A
    !> steady pressure contribution ( overwritten below )
    el%pres   = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
               ( c_m(i_l,1) * cos(a_v(i_l)) +  c_m(i_l,2) * sin(a_v(i_l)) )

    if ( load_avl ) then

      ! === Kutta-Joukowski theorem ~ AVL ===
      !> Inviscid contribution ~ AVL
      call el % compute_dforce_jukowski ( elems_tot , elems_wake )

      !> Add viscous contribution
      el % dforce = el % dforce + &
           0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * el%area * &
           c_m(i_l,2) * ( sin(el%al_ctr_pt) * el%nor      +  &
                          cos(el%al_ctr_pt) * el%tang_cen )

      ! === Update AOAs and aerodynamic coefficients ===
      !> unit vector in the direction of the relative velocity (_d for drag)
      !  and in the direction of the lift (_l for lift)
      e_l = el%nor * cos( el%al_ctr_pt ) - el%tang_cen * sin( el%al_ctr_pt )
      e_d = el%nor * sin( el%al_ctr_pt ) + el%tang_cen * cos( el%al_ctr_pt )
      e_l = e_l / norm2(e_l)
      e_d = e_d / norm2(e_d)

      !> Unsteady contribution
      el%dforce = el%dforce &
                  - sim_param%rho_inf * el%area * el%dGamma_dt &
                  * e_l ! lift direction

      !> Update aerodynamic coefficients and AOA, and pressure
      c_m(i_l,1) = sum( el % dforce * e_l ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      c_m(i_l,2) = sum( el % dforce * e_d ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      a_v(i_l) = el % al_ctr_pt

      !> overwrite pressure field to take into account unsteady contributions
      el%pres = sum( el%dforce * el%nor ) / el%area

    else

      ! === elementary force = p*n + tangential contribution from L,D ===
      el%dforce = ( el%nor * el%pres + &
                    el%tang_cen * &
                    0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * ( &
                   -c_m(i_l,1) * sin(a_v(i_l)) + c_m(i_l,2) * cos(a_v(i_l)) &
                   ) ) * el%area

      ! === Update AOAs and aerodynamic coefficients ===
      !> unit vector in the direction of the relative velocity (_d for drag)
      !  and in the direction of the lift (_l for lift)
      e_l = el%nor * cos( a_v(i_l) ) - el%tang_cen * sin( a_v(i_l) )
      e_d = el%nor * sin( a_v(i_l) ) + el%tang_cen * cos( a_v(i_l) )
      e_l = e_l / norm2(e_l)
      e_d = e_d / norm2(e_d)

      !> Unsteady contribution
      el%dforce = el%dforce &
                  - sim_param%rho_inf * el%area * el%dGamma_dt &
                  * e_l ! lift direction

      !> Update aerodynamic coefficients and AOA, and pressure
      c_m(i_l,1) = sum( el % dforce * e_l ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      c_m(i_l,2) = sum( el % dforce * e_d ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )

      !> overwrite pressure field to take into account unsteady contributions
      el%pres = sum( el%dforce * el%nor ) / el%area

    end if

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

  diff = maxval( abs( f_Gam ) ) 
  if(sim_param%debug_level .ge. 3) then
    write(msg,*) 'iterations: ',ic; call printout(trim(msg))
    write(msg,*) 'diff',diff; call printout(trim(msg))
    ! write(msg,*) 'diff/max_mag_ll:',diff/max_mag_ll; call printout(trim(msg))
  endif

  ! === Deallocations ==
  deallocate(Gam, Lam, grad_Gam, grad_Lam, df_dGam, f_Gam, Gamma_old,d_Gam)
  deallocate(a_v, c_m, u_v, dcl_v, mu_v)
  deallocate(a_ik)


end subroutine solve_liftlin_newton

!----------------------------------------------------------------------
subroutine update_kernel( elems_ll , kernel_in , dcl_v , kernel_out , f_chord_v )
 type(t_liftlin_p)     , intent(in)  :: elems_ll(:)
 real(wp), allocatable , intent(in)  :: dcl_v(:)
 real(wp), allocatable , intent(in)  :: kernel_in(:,:)
 real(wp), allocatable , intent(out) :: kernel_out(:,:)

 integer , allocatable :: global_to_local_ll(:)
 real(wp), allocatable :: p_v(:)

 real(wp), allocatable , intent(out) :: f_chord_v(:)
 real(wp) :: sigma

 integer :: i_l , j_l , n_l , j , n_neigh , max_ll_id

 real(wp) :: p_tilde = 0.5_wp

 n_l = size(elems_ll)

 allocate( f_chord_v(n_l) ) ; f_chord_v = 0.0_wp
 allocate(p_v(n_l))
 p_v = matmul( kernel_in , dcl_v )
 p_v = abs( dcl_v - p_v )

 p_v = p_v / p_tilde
 
 do i_l = 1 , n_l ! cubic function
   if ( p_v(i_l) .lt. 1.0_wp ) then
     f_chord_v(i_l) = sim_param%llArtificialViscosity * ( &
                   3.0_wp * p_v(i_l)**2.0_wp - 2.0_wp * p_v(i_l)**3.0_wp ) 
   else
     f_chord_v(i_l) = sim_param%llArtificialViscosity
   end if 
 end do
 
 deallocate(p_v)

 !> array linking global id of ll to ll-indexing
 max_ll_id = elems_ll(1)%p%id
 do i_l = 1 , n_l
   max_ll_id = max( max_ll_id , elems_ll(i_l)%p%id )
 end do 

 allocate(global_to_local_ll(max_ll_id)) ; global_to_local_ll = 0
 do i_l = 1 , n_l
   global_to_local_ll( elems_ll(i_l)%p%id ) = i_l
 end do

 ! === Kernel ===
 allocate(kernel_out(n_l,n_l)) ; kernel_out = 0.0_wp 

 !> Normal distribution
 do i_l = 1 , n_l

   if ( f_chord_v(i_l) .eq. 0.0_wp ) then

     kernel_out(i_l,i_l) = 1.0_wp 

   else

     sigma = 2.0_wp * f_chord_v(i_l) * elems_ll(i_l)%p%chord
     do j_l = 1 , n_l
    
       kernel_out(i_l,j_l) = exp( - 0.5_wp * norm2( elems_ll(j_l)%p%cen &
                                                  - elems_ll(i_l)%p%cen )**2.0_wp / &
                              sigma**2.0_wp )
    
     end do

   end if

   !> normalisation
   kernel_out(i_l,:) = kernel_out(i_l,:) / sum(kernel_out(i_l,:))

 end do

 ! Deallocations
 deallocate(global_to_local_ll) !,f_chord_v)

end subroutine update_kernel

!----------------------------------------------------------------------

subroutine build_ll_array_optim( elems_ll , M , kernel )
 type(t_liftlin_p)     , intent(in)  :: elems_ll(:)
 real(wp), allocatable , intent(out) :: M(:,:)
 real(wp), allocatable , intent(out) :: kernel(:,:)

 integer , allocatable :: global_to_local_ll(:)
 real(wp), allocatable :: p_v(:)

 real(wp), allocatable :: f_chord_v(:)
 real(wp) :: sigma

 integer :: i_l , j_l , n_l , j , n_neigh , max_ll_id


 n_l = size(elems_ll)

 allocate(     M(n_l,n_l) ) ;         M = 0.0_wp 
 allocate( f_chord_v(n_l) ) ; f_chord_v = 0.0_wp

 f_chord_v = sim_param%llArtificialViscosity

 !> array linking global id of ll to ll-indexing
 max_ll_id = elems_ll(1)%p%id
 do i_l = 1 , n_l
   max_ll_id = max( max_ll_id , elems_ll(i_l)%p%id )
 end do 

 allocate(global_to_local_ll(max_ll_id)) ; global_to_local_ll = 0
 do i_l = 1 , n_l
   global_to_local_ll( elems_ll(i_l)%p%id ) = i_l
 end do

 do i_l = 1 , n_l

   n_neigh = 0

   !> Out-Of-Diagonal elements
   do j = 1 , size( elems_ll(i_l)%p%neigh )
     if ( associated(elems_ll(i_l)%p%neigh(j)%p) ) then
       
       M(i_l, &
         global_to_local_ll( elems_ll(i_l)%p%neigh(j)%p%id ) ) = 1.0_wp
       n_neigh = n_neigh + 1

     end if
   end do 

   !> Diagonal elements
   if (     n_neigh .eq. 1 ) then ; M(i_l,i_l) = -1.0_wp 
   elseif ( n_neigh .eq. 2 ) then ; M(i_l,i_l) = -2.0_wp
   else ;  write(*,*) ' Warning: elems_ll(', i_l , ') has no neighbor.'
   end if

 end do

 ! M = matmul( M , M )


 ! === Kernel ===
 allocate(kernel(n_l,n_l)) ; kernel = 0.0_wp 

! !> Identity
! do i_l = 1 , size(elems_ll) ; kernel(i_l,i_l) = 1.0_wp ; end do
!
! !> Rough kernel [ 0.25 , 0.5 , 0.25 ] (or  [ 0.5 , 0.5])
! do i_l = 1 , n_l
!
!   n_neigh = 0
!
!   !> Out-Of-Diagonal elements
!   do j = 1 , size( elems_ll(i_l)%p%neigh )
!     if ( associated(elems_ll(i_l)%p%neigh(j)%p) ) then
!       
!       kernel( i_l , &
!         global_to_local_ll( elems_ll(i_l)%p%neigh(j)%p%id ) ) = 0.25_wp
!       n_neigh = n_neigh + 1
!
!     end if
!   end do 
!
!   kernel(i_l,i_l) = 1.0_wp
!
!   !> Diagonal element
!   if ( n_neigh .eq. 1 ) then
!
!       kernel(i_l,i_l) = 0.25_wp
!       kernel(i_l,:) = 2.0_wp * kernel(i_l,:)
!
!   elseif ( n_neigh .eq. 2 ) then
!
!       kernel(i_l,i_l) = 0.50_wp
!
!   else
!     write(*,*) ' Warning: elems_ll(', i_l , ') has no neighbor.'
!   end if
!
! end do

 !> Normal distribution
 do i_l = 1 , n_l

   if ( f_chord_v(i_l) .eq. 0.0_wp ) then

     kernel(i_l,i_l) = 1.0_wp 

   else

     sigma = 2.0_wp * f_chord_v(i_l) * elems_ll(i_l)%p%chord
     do j_l = 1 , n_l
    
       kernel(i_l,j_l) = exp( - 0.5_wp * norm2( elems_ll(j_l)%p%cen &
                                              - elems_ll(i_l)%p%cen )**2.0_wp / &
                              sigma**2.0_wp )
    
     end do

   end if

   !> normalisation
   kernel(i_l,:) = kernel(i_l,:) / sum(kernel(i_l,:))

 end do

 ! Deallocations
 deallocate(global_to_local_ll,f_chord_v)


end subroutine build_ll_array_optim

!----------------------------------------------------------------------
subroutine solve_liftlin_optim_regul( &
                         elems_ll, elems_tot, &
                         elems_impl, elems_ad, &
                         elems_wake, elems_vort, &
                         airfoil_data, it, M_array )

  type(t_liftlin_p)  , intent(inout) :: elems_ll(:)
  type(t_pot_elem_p) , intent(in)    :: elems_tot(:)
  type(t_impl_elem_p), intent(in)    :: elems_impl(:)
  type(t_expl_elem_p), intent(in)    :: elems_ad(:)
  type(t_pot_elem_p) , intent(in)    :: elems_wake(:)
  type(t_vort_elem_p), intent(in)    :: elems_vort(:)
  type(t_aero_tab)   , intent(in)    :: airfoil_data(:)
  integer            , intent(in)    :: it
  real(wp)           , intent(in)    :: M_array(:,:)

  real(wp) :: vel(3) , v(3) , up(3)
  real(wp), allocatable :: vel_w(:,:)
  real(wp) :: unorm , alpha , alpha_2d
  !> Gam, Lam, grad_Gam, grad_Lam (for gradient descent)
  real(wp), allocatable :: Gam(:) , Gamma_old(:)
  real(wp), allocatable :: grad_Gam(:) , df_dGam(:,:) , f_Gam(:)
  real(wp), allocatable :: a_ik(:,:,:)
  !> uinf
  real(wp) :: uinf(3)
  real(wp) :: mach , reynolds ! mach and reynolds number for each el
  ! arrays used to store aoa, aero coeff and vel
  real(wp) , allocatable :: a_v(:)   ! size(elems_ll)
  real(wp) , allocatable :: c_m(:,:) ! size(elems_ll) , 3
  real(wp) , allocatable :: u_v(:)   ! size(elems_ll)
  real(wp) , allocatable :: dcl_v(:) ! size(elems_ll)
  real(wp) , allocatable :: aero_coeff(:)
  real(wp) :: dcl_da , Jopt
  real(wp) :: diff , max_mag_ll
  !> sim_param: fixed point algorithm for ll
  real(wp) :: fp_tol , fp_damp
  integer  :: fp_maxIter
  !> sim_param: load computation
  logical :: load_avl
  real(wp) :: e_l(3) , e_d(3)

  real(wp) :: lambda   = 10.0_wp     ! regularisation parameter
  real(wp) :: grad_tol = 0.00001_wp  ! convergence parameter
  ! TODO: provide them as an input
  real(wp) :: Jopt_old , dJopt

  type(t_liftlin), pointer :: el
  integer :: i_l , j , ic

  real(wp) , allocatable :: diff_v(:)
  character(len=4) :: msg_4
  character(len=max_char_len) :: msg
  character(len=*), parameter :: this_sub_name = 'solve_liftlin_optim_regul'
 
  allocate( diff_v(sim_param%llMaxIter+1) ) ; diff_v = 0.0_wp

  ! params of the fixed point iterations
  fp_tol     = sim_param%llTol
  fp_damp    = sim_param%llDamp
  fp_maxIter = sim_param%llMaxIter
  ! param for load computation
  load_avl   = sim_param%llLoadsAVL
 
  uinf = sim_param%u_inf

  !> allocate and fill Gamma_old array of the ll intensity at previous dt
  allocate(Gamma_old(size(elems_ll)))
  do i_l = 1 , size(elems_ll)
    Gamma_old(i_l) = elems_ll(i_l)%p%mag
  end do

  ! debug --- 
  write(*,*) ' Gamma_old '
  do i_l = 1 , size(elems_ll)
    write(*,*) Gamma_old(i_l)
  end do
  ! debug --- 

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

  ! allocate array containing aoa, aero coeffs and relative velocity
  allocate(   a_v( size(elems_ll)  )) ;   a_v = 0.0_wp
  allocate(   c_m( size(elems_ll),3)) ;   c_m = 0.0_wp
  allocate(   u_v( size(elems_ll)  )) ;   u_v = 0.0_wp
  allocate( dcl_v( size(elems_ll)  )) ; dcl_v = 0.0_wp

  ! Remove the "out-of-plane" component of the relative velocity:
  ! 2d-velocity to enter the airfoil look-up-tables
  ! IS THIS LOOP USED? (u_v) seems to be overwritten few lines down)
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

  ! === Allocate and initialize Gam and Lam arrays ===
  allocate(     Gam(size(elems_ll))) ;      Gam = Gamma_old 
  allocate(grad_Gam(size(elems_ll))) ; grad_Gam = 0.0_wp
  allocate(   f_Gam(size(elems_ll))) ;    f_Gam = 0.0_wp
  allocate( df_dGam(size(elems_ll),size(elems_ll))) ;  df_dGam = 0.0_wp

  ! === Save AIC for derivatives d(alpha)/d(gamma) ===
  allocate( a_ik( size(elems_ll) , size(elems_ll) , 3 ) ) ; a_ik = 0.0_wp
  !> set intensity of the ll elems to 1, to compute AICs
  do i_l = 1 , size(elems_ll)
    elems_ll(i_l)%p%mag = 1.0_wp
  end do 
  !> compute AICs for derivatives
  do i_l = 1 , size(elems_ll)
    do j = 1 , size(elems_ll)
      call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf, a_ik(i_l,j,:) )
    end do
  end do
  a_ik = a_ik / ( 4.0_wp * pi )
  !> set intensity of the ll elems to their actual value
  do i_l = 1 , size(elems_ll)
    elems_ll(i_l)%p%mag = Gamma_old(i_l)
  end do 
 
  ! === Iterative algorithm ===
  do ic = 1 , fp_maxIter

    !> Reset diff and max_mag_ll variables
    diff = 0.0_wp ; max_mag_ll = 0.0_wp

!$omp parallel do private(i_l, el, j, v, vel, up, unorm, alpha, alpha_2d, mach, &
!$omp& reynolds, aero_coeff, dcl_da) schedule(dynamic,4)
    do i_l = 1 , size(elems_ll)
 
      ! === Compute velocity ===
      !> LL influence ( to be divided by 4*pi )
      vel = 0.0_wp
      do j = 1 , size(elems_ll)
        call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
        vel = vel + v
      enddo

      el => elems_ll(i_l)%p

      !> overall relative velocity computed in the centre of the ll elem
      vel = vel/(4.0_wp*pi) + uinf - el%ub + vel_w(:,i_l)
      !> "effective" velocity = proj. of vel in the n-t plane
      up =  el%nor*sum(el%nor*vel) + el%tang_cen*sum(el%tang_cen*vel)
      unorm = norm2(up)      ! velocity w/o induced velocity
      
      ! === Angle of incidence (full velocity) ===
      alpha = atan2(sum(up*el%nor), sum(up*el%tang_cen))
      alpha = alpha * 180.0_wp/pi  ! .c81 tables defined with angles in [deg]

      ! === Piszkin, Lewinski (1976) LL model for swept wings ===
      ! the control point is approximately at 3/4 of the chord, but the induced
      ! angle of incidence needs to be modified, introducing a "2D correction"
      !
      !> "2D correction" of the induced angle
      alpha_2d = el%mag / ( pi * el%chord * unorm ) *180.0_wp/pi
      alpha = alpha - alpha_2d 
      ! =========================================================

      ! === Mach and Reynolds numbers ===
      mach     = unorm / sim_param%a_inf      
      reynolds = sim_param%rho_inf * unorm * el%chord / sim_param%mu_inf    

      ! === Read the aero coeff from .c81 tables ===
      call interp_aero_coeff ( airfoil_data,  el%csi_cen, el%i_airfoil , &
                                 (/alpha, mach, reynolds/) , sim_param , &
                                                  aero_coeff ,  dcl_da )

      ! === Fill arrays for iterative method ===
      u_v(  i_l)   = unorm
      a_v(  i_l)   = alpha * pi / 180.0_wp
      c_m(  i_l,:) = aero_coeff
      dcl_v(i_l)   = dcl_da

    end do
!$omp end parallel do

    !> Compute gradient
    df_dGam = 0.0_wp ! compute_df_dlam()
    f_Gam   = 0.0_wp ! compute_f_lam()
    do i_l = 1 , size(elems_ll)
      f_Gam(i_l) = 0.5_wp * u_v(i_l) * c_m(i_l,1) * elems_ll(i_l)%p%chord  
    end do
    do i_l = 1 , size(elems_ll)
      do j = 1 , size(elems_ll)
        df_dGam(i_l,j) = 0.5_wp * elems_ll(i_l)%p%chord * &
                         dcl_v(i_l) * &
                         sum( a_ik(i_l,j,:) * ( &
                         elems_ll(i_l)%p%nor      * cos(a_v(i_l)) &
                       - elems_ll(i_l)%p%tang_cen * sin(a_v(i_l)) ) )
      end do
      ! Update diagonal component
      df_dGam(i_l,i_l) = 1.0_wp + df_dGam(i_l,i_l)
    end do

    ! Identity weight matrix for the nonlinear system error
    grad_Gam = matmul( transpose(df_dGam) , Gam+f_Gam ) + &
       lambda* matmul( transpose(M_array) , matmul( M_array , Gam ) )

!   ! debug ---
!   write(*,*) ' minval , maxval : '
!   write(*,*) ' M_array : ' , minval(M_array ) , maxval(M_array )
!   write(*,*) ' Gam     : ' , minval(Gam     ) , maxval(Gam     )
!   write(*,*) ' Lam     : ' , minval(Lam     ) , maxval(Lam     )
!   write(*,*) ' grad_Gam: ' , minval(grad_Gam) , maxval(grad_Gam)
!   write(*,*) ' grad_Lam: ' , minval(grad_Lam) , maxval(grad_Lam)
!   ! debug ---

!   ! debug ---
!   write(*,*) ' i_l , alpha[rad] , alpha[deg] , cl , dcl_da '
!   do i_l = 1 , size(elems_ll)
!     write(*,*) i_l , a_v(i_l) , a_v(i_l)*180.0_wp/pi , c_m(i_l,1) , dcl_v(i_l) 
!   end do
!   write(*,*) ' stop in mod_liftlin.f90, l.563' ; stop
!   ! debug ---

    ! debug ---
    write(*,*) '         i_l ,     alpha[deg] ,                  cl ,  &
    &                   dcl_da ,                   Gam ,               &
    &     f_gam ,                  Gam+f_Gam ,              grad_Gam'
    do i_l = 1 , size(elems_ll)
      write(*,*) i_l , a_v(i_l)*180.0_wp/pi , c_m(i_l,1) , dcl_v(i_l) , &
             Gam(i_l) , f_gam(i_l) , Gam(i_l)+f_Gam(i_l) , grad_Gam(i_l)
    end do
    ! debug ---

    diff = 0.0_wp
    !> Update ll intensity
    do i_l = 1 , size(elems_ll)

      !>> Gradient descent
      !> Update ll intensity
      Gam(i_l) = Gam(i_l) - 0.001_wp * grad_Gam(i_l)
      elems_ll(i_l)%p%mag = Gam(i_l)
      
      !>> Other optimization algorithms ???
      ! ...

      !> Update max_mag_ll 
      max_mag_ll = max(max_mag_ll,abs(elems_ll(i_l)%p%mag))

    enddo

    ! === Objective function and Stopping Criterion === 
    Jopt = 0.5_wp * ( norm2( Gam + f_Gam ) ) **2.0_wp + &
           0.5_wp * ( sum( Gam * matmul( M_array , Gam ) ) ) * lambda

    diff_v(ic) = Jopt

!   ! debug ---
!   write(*,*) ' J(',ic,'): ' , Jopt
!   ! debug ---

    if ( ic .gt. 1 ) then

      !> Update
      dJopt = abs( Jopt - Jopt_old )

!     write(*,*) ' Jopt, Jopt_old : ' , Jopt , Jopt_old

      !> Check convergence
      if ( dJopt .lt. grad_tol ) exit ! convergence

    end if
    
    Jopt_old = Jopt

  end do
 
  write(msg_4,'(I4.4)') it
  open(unit=21,file=trim(sim_param%basename_debug)//'_conv_'//msg_4//'.dat')
  do i_l = 1 , ic-1
    write(21,*) diff_v(i_l)
  end do
  close(21)

! write(*,*) ' stop in mod_liftlin.f90, l.624. ' ; stop
  ! === Update dGamma_dt field ===
  do i_l = 1,size(elems_ll)
    elems_ll(i_l)%p%dGamma_dt = ( elems_ll(i_l)%p%mag - Gamma_old(i_l) ) / &
                                  sim_param%dt 
  end do

  ! === Loads computation ===
  ! compute LL sectional loads from singularity intensities. 
  do i_l = 1,size(elems_ll)

   !select type(el => elems_ll(i_l)%p)
   !type is(t_liftlin)
   el => elems_ll(i_l)%p
    ! avg delta_p = \vec{F}.\vec{n} / A = ( L*cos(al)+D*sin(al) ) / A
    !> steady pressure contribution ( overwritten below )
    el%pres   = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
               ( c_m(i_l,1) * cos(a_v(i_l)) +  c_m(i_l,2) * sin(a_v(i_l)) )

    if ( load_avl ) then

      ! === Kutta-Joukowski theorem ~ AVL ===
      !> Inviscid contribution ~ AVL
      call el % compute_dforce_jukowski ( elems_tot , elems_wake )

      !> Add viscous contribution
      el % dforce = el % dforce + &
           0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * el%area * &
           c_m(i_l,2) * ( sin(el%al_ctr_pt) * el%nor      +  &
                          cos(el%al_ctr_pt) * el%tang_cen )

      ! === Update AOAs and aerodynamic coefficients ===
      !> unit vector in the direction of the relative velocity (_d for drag)
      !  and in the direction of the lift (_l for lift)
      e_l = el%nor * cos( el%al_ctr_pt ) - el%tang_cen * sin( el%al_ctr_pt )
      e_d = el%nor * sin( el%al_ctr_pt ) + el%tang_cen * cos( el%al_ctr_pt )
      e_l = e_l / norm2(e_l)
      e_d = e_d / norm2(e_d)

      !> Unsteady contribution
      el%dforce = el%dforce &
                  - sim_param%rho_inf * el%area * el%dGamma_dt &
                  * e_l ! lift direction

      !> Update aerodynamic coefficients and AOA, and pressure
      c_m(i_l,1) = sum( el % dforce * e_l ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      c_m(i_l,2) = sum( el % dforce * e_d ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      a_v(i_l) = el % al_ctr_pt

      !> overwrite pressure field to take into account unsteady contributions
      el%pres = sum( el%dforce * el%nor ) / el%area

    else

      ! === elementary force = p*n + tangential contribution from L,D ===
      el%dforce = ( el%nor * el%pres + &
                    el%tang_cen * &
                    0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * ( &
                   -c_m(i_l,1) * sin(a_v(i_l)) + c_m(i_l,2) * cos(a_v(i_l)) &
                   ) ) * el%area

      ! === Update AOAs and aerodynamic coefficients ===
      !> unit vector in the direction of the relative velocity (_d for drag)
      !  and in the direction of the lift (_l for lift)
      e_l = el%nor * cos( a_v(i_l) ) - el%tang_cen * sin( a_v(i_l) )
      e_d = el%nor * sin( a_v(i_l) ) + el%tang_cen * cos( a_v(i_l) )
      e_l = e_l / norm2(e_l)
      e_d = e_d / norm2(e_d)

      !> Unsteady contribution
      el%dforce = el%dforce &
                  - sim_param%rho_inf * el%area * el%dGamma_dt &
                  * e_l ! lift direction

      !> Update aerodynamic coefficients and AOA, and pressure
      c_m(i_l,1) = sum( el % dforce * e_l ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      c_m(i_l,2) = sum( el % dforce * e_d ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )

      !> overwrite pressure field to take into account unsteady contributions
      el%pres = sum( el%dforce * el%nor ) / el%area

    end if

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

  ! === Output Statistics ===
  if(sim_param%debug_level .ge. 3) then
    write(msg,*) 'iterations: ',ic; call printout(trim(msg))
    write(msg,*) ' J = ', Jopt; call printout(trim(msg))
  endif

  ! debug --- 
  write(*,*) ' Gam '
  do i_l = 1 , size(elems_ll)
    write(*,*) Gam(i_l) , elems_ll(i_l)%p%mag
  end do
  ! debug --- 

! write(*,*) ' stop in mod_liftlin.f90, l 1130. ' ; stop

  ! === Deallocations ==
  deallocate(Gam, grad_Gam, df_dGam, f_Gam, Gamma_old)
  deallocate(a_v, c_m, u_v, dcl_v)
  deallocate(a_ik)


end subroutine solve_liftlin_optim_regul

!----------------------------------------------------------------------

subroutine solve_liftlin_optim( &
                         elems_ll, elems_tot, &
                         elems_impl, elems_ad, &
                         elems_wake, elems_vort, &
                         airfoil_data, it, M_array )

  type(t_liftlin_p)  , intent(inout) :: elems_ll(:)
  type(t_pot_elem_p) , intent(in)    :: elems_tot(:)
  type(t_impl_elem_p), intent(in)    :: elems_impl(:)
  type(t_expl_elem_p), intent(in)    :: elems_ad(:)
  type(t_pot_elem_p) , intent(in)    :: elems_wake(:)
  type(t_vort_elem_p), intent(in)    :: elems_vort(:)
  type(t_aero_tab)   , intent(in)    :: airfoil_data(:)
  integer            , intent(in)    :: it
  real(wp)           , intent(in)    :: M_array(:,:)

  real(wp) :: vel(3) , v(3) , up(3)
  real(wp), allocatable :: vel_w(:,:)
  real(wp) :: unorm , alpha , alpha_2d
  !> Gam, Lam, grad_Gam, grad_Lam (for gradient descent)
  real(wp), allocatable :: Gam(:) , Lam(:) , Gamma_old(:)
  real(wp), allocatable :: grad_Gam(:) , grad_Lam(:) , df_dGam(:,:) , f_Gam(:)
  real(wp), allocatable :: a_ik(:,:,:)
  !> uinf
  real(wp) :: uinf(3)
  real(wp) :: mach , reynolds ! mach and reynolds number for each el
  ! arrays used to store aoa, aero coeff and vel
  real(wp) , allocatable :: a_v(:)   ! size(elems_ll)
  real(wp) , allocatable :: c_m(:,:) ! size(elems_ll) , 3
  real(wp) , allocatable :: u_v(:)   ! size(elems_ll)
  real(wp) , allocatable :: dcl_v(:) ! size(elems_ll)
  real(wp) , allocatable :: aero_coeff(:)
  real(wp) :: dcl_da , Jopt
  real(wp) :: diff , max_mag_ll
  !> sim_param: fixed point algorithm for ll
  real(wp) :: fp_tol , fp_damp
  integer  :: fp_maxIter
  !> sim_param: load computation
  logical :: load_avl

  type(t_liftlin), pointer :: el
  integer :: i_l , j , ic

  character(len=max_char_len) :: msg
  character(len=*), parameter :: this_sub_name = 'solve_liftlin_optim'
 
  ! params of the fixed point iterations
  fp_tol     = sim_param%llTol
  fp_damp    = sim_param%llDamp
  fp_maxIter = sim_param%llMaxIter
  ! param for load computation
  load_avl   = sim_param%llLoadsAVL
 
  uinf = sim_param%u_inf

  !> allocate and fill Gamma_old array of the ll intensity at previous dt
  allocate(Gamma_old(size(elems_ll)))
  do i_l = 1 , size(elems_ll)
    Gamma_old(i_l) = elems_ll(i_l)%p%mag
  end do

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

  ! allocate array containing aoa, aero coeffs and relative velocity
  allocate(   a_v( size(elems_ll)  )) ;   a_v = 0.0_wp
  allocate(   c_m( size(elems_ll),3)) ;   c_m = 0.0_wp
  allocate(   u_v( size(elems_ll)  )) ;   u_v = 0.0_wp
  allocate( dcl_v( size(elems_ll)  )) ; dcl_v = 0.0_wp

  ! Remove the "out-of-plane" component of the relative velocity:
  ! 2d-velocity to enter the airfoil look-up-tables
  ! IS THIS LOOP USED? (u_v) seems to be overwritten few lines down)
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

  ! === Allocate and initialize Gam and Lam arrays ===
  allocate(     Gam(size(elems_ll))) ;      Gam = Gamma_old 
  allocate(     Lam(size(elems_ll))) ;      Lam = 1.0_wp
  allocate(grad_Gam(size(elems_ll))) ; grad_Gam = 0.0_wp
  allocate(grad_Lam(size(elems_ll))) ; grad_Lam = 0.0_wp
  allocate(   f_Gam(size(elems_ll))) ;    f_Gam = 0.0_wp
  allocate( df_dGam(size(elems_ll),size(elems_ll))) ;  df_dGam = 0.0_wp

  ! === Save AIC for derivatives d(alpha)/d(gamma) ===
  allocate( a_ik( size(elems_ll) , size(elems_ll) , 3 ) ) ; a_ik = 0.0_wp
  !> set intensity of the ll elems to 1, to compute AICs
  do i_l = 1 , size(elems_ll)
    elems_ll(i_l)%p%mag = 1.0_wp
  end do 
  !> compute AICs for derivatives
  do i_l = 1 , size(elems_ll)
    do j = 1 , size(elems_ll)
      call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf, a_ik(i_l,j,:) )
    end do
  end do
  a_ik = a_ik / ( 4.0_wp * pi )
  !> set intensity of the ll elems to their actual value
  do i_l = 1 , size(elems_ll)
    elems_ll(i_l)%p%mag = Gamma_old(i_l)
  end do 

  ! === Iterative algorithm ===
  do ic = 1 , fp_maxIter

    !> Reset diff and max_mag_ll variables
    diff = 0.0_wp ; max_mag_ll = 0.0_wp

!$omp parallel do private(i_l, el, j, v, vel, up, unorm, alpha, alpha_2d, mach, &
!$omp& reynolds, aero_coeff, dcl_da) schedule(dynamic,4)
    do i_l = 1 , size(elems_ll)
 
      ! === Compute velocity ===
      !> LL influence ( to be divided by 4*pi )
      vel = 0.0_wp
      do j = 1 , size(elems_ll)
        call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
        vel = vel + v
      enddo

      el => elems_ll(i_l)%p

      !> overall relative velocity computed in the centre of the ll elem
      vel = vel/(4.0_wp*pi) + uinf - el%ub + vel_w(:,i_l)
      !> "effective" velocity = proj. of vel in the n-t plane
      up =  el%nor*sum(el%nor*vel) + el%tang_cen*sum(el%tang_cen*vel)
      unorm = norm2(up)      ! velocity w/o induced velocity
      
      ! === Angle of incidence (full velocity) ===
      alpha = atan2(sum(up*el%nor), sum(up*el%tang_cen))
      alpha = alpha * 180.0_wp/pi  ! .c81 tables defined with angles in [deg]

      ! === Piszkin, Lewinski (1976) LL model for swept wings ===
      ! the control point is approximately at 3/4 of the chord, but the induced
      ! angle of incidence needs to be modified, introducing a "2D correction"
      !
      !> "2D correction" of the induced angle
      alpha_2d = el%mag / ( pi * el%chord * unorm ) *180.0_wp/pi
      alpha = alpha - alpha_2d 
      ! =========================================================

      ! === Mach and Reynolds numbers ===
      mach     = unorm / sim_param%a_inf      
      reynolds = sim_param%rho_inf * unorm * el%chord / sim_param%mu_inf    

      ! === Read the aero coeff from .c81 tables ===
      call interp_aero_coeff ( airfoil_data,  el%csi_cen, el%i_airfoil , &
                                 (/alpha, mach, reynolds/) , sim_param , &
                                                  aero_coeff ,  dcl_da )

      ! === Fill arrays for iterative method ===
      u_v(  i_l)   = unorm
      a_v(  i_l)   = alpha * pi / 180.0_wp
      c_m(  i_l,:) = aero_coeff
      dcl_v(i_l)   = dcl_da

    end do
!$omp end parallel do

!   ! debug ---
!   write(*,*) ' i_l , alpha[rad] , alpha[deg] , cl , dcl_da '
!   do i_l = 1 , size(elems_ll)
!     write(*,*) i_l , a_v(i_l) , a_v(i_l)*180.0_wp/pi , c_m(i_l,1) , dcl_v(i_l) 
!   end do
!   write(*,*) ' stop in mod_liftlin.f90, l.563' ; stop
!   ! debug ---

    !> Compute gradient
    df_dGam = 0.0_wp ! compute_df_dlam()
    f_Gam   = 0.0_wp ! compute_f_lam()
    do i_l = 1 , size(elems_ll)
      f_Gam(i_l) = 0.5_wp * u_v(i_l) * c_m(i_l,1) * elems_ll(i_l)%p%chord  
    end do
    do i_l = 1 , size(elems_ll)
      do j = 1 , size(elems_ll)
        df_dGam(i_l,j) = 0.5_wp * elems_ll(i_l)%p%chord * &
                         dcl_v(i_l) * &
                         sum( a_ik(i_l,j,:) * ( &
                         elems_ll(i_l)%p%nor      * cos(a_v(i_l)) &
                       - elems_ll(i_l)%p%tang_cen * sin(a_v(i_l)) ) )
      end do
    end do

    grad_Gam = matmul( M_array , Gam ) + Lam + matmul(transpose(df_dGam), Lam)
    grad_Lam = Gam + f_Gam

!   ! debug ---
!   write(*,*) ' minval , maxval : '
!   write(*,*) ' M_array : ' , minval(M_array ) , maxval(M_array )
!   write(*,*) ' Gam     : ' , minval(Gam     ) , maxval(Gam     )
!   write(*,*) ' Lam     : ' , minval(Lam     ) , maxval(Lam     )
!   write(*,*) ' grad_Gam: ' , minval(grad_Gam) , maxval(grad_Gam)
!   write(*,*) ' grad_Lam: ' , minval(grad_Lam) , maxval(grad_Lam)
!   ! debug ---

!   ! debug ---
!   write(*,*) ' i_l , Gam , Lam , f_gam , grad_Gam , grad_Lam '
!   do i_l = 1 , size(elems_ll)
!     write(*,*) i_l , Gam(i_l) , Lam(i_l) , f_gam(i_l) , grad_Gam(i_l) , grad_Lam(i_l)
!   end do
!   ! debug ---

    diff = 0.0_wp
    !> Update ll intensity
    do i_l = 1 , size(elems_ll)

      diff = max( diff ,               0.001_wp * abs(grad_Gam(i_l) ) )
      !>> Gradient descent
      !> Update ll intensity
      Gam(i_l) = Gam(i_l)            - 0.001_wp * grad_Gam(i_l)
      elems_ll(i_l)%p%mag = Gam(i_l)
      !> Update constraint value
      Lam(i_l) = Lam(i_l)            - 0.001_wp * grad_Lam(i_l)
      
      !>> Other optimization algorithms ???
      ! ...

      !> Update max_mag_ll 
      max_mag_ll = max(max_mag_ll,abs(elems_ll(i_l)%p%mag))

    enddo

!   !> Stopping criterion
!   write(*,*) diff
!   write(*,*) max_mag_ll
!   write(*,*) fp_tol
    if ( diff .lt. fp_tol*max_mag_ll ) exit ! convergence

    Jopt = 0.5_wp*sum( Gam*matmul(M_array,Gam) ) + &
           sum ( Lam * ( Gam + f_Gam ) )
    write(*,*) ' J(',ic,'): ' , Jopt

    write(*,*) ' stop in mod_liftlin.f90, l.624. ' ; stop

  end do

  ! === Load computations === ( same as solve_liftlin() routine )
  ! ...


  ! === Deallocations ==
  deallocate(Gam, Lam, grad_Gam, grad_Lam, df_dGam, f_Gam, Gamma_old)
  deallocate(a_v, c_m, u_v, dcl_v)
  deallocate(a_ik)


end subroutine solve_liftlin_optim

!----------------------------------------------------------------------

!> Solve the lifting line, in an iterative way
!!
!! The lifting line solution is not obtained from the solution of the linear
!! system. It is fully explicit, but by being nonlinear requires an
!! iterative solution.
subroutine solve_liftlin_piszkin( &
                         elems_ll, elems_tot, &
                         elems_impl, elems_ad, &
                         elems_wake, elems_vort, &
                         airfoil_data, it, al_kernel)
 !type(t_expl_elem_p), intent(inout) :: elems_ll(:)
 type(t_liftlin_p), intent(inout) :: elems_ll(:)
 type(t_pot_elem_p),  intent(in)    :: elems_tot(:)
 type(t_impl_elem_p), intent(in)    :: elems_impl(:)
 type(t_expl_elem_p), intent(in)    :: elems_ad(:)
 type(t_pot_elem_p),  intent(in)    :: elems_wake(:)
 type(t_vort_elem_p), intent(in)    :: elems_vort(:)
 type(t_aero_tab),    intent(in)    :: airfoil_data(:)
 real(wp)           , allocatable, intent(inout) :: al_kernel(:,:)
 real(wp)           , allocatable                :: al_kernel_out(:,:)
 real(wp) :: uinf(3)
 integer  :: i_l, j, ic
 real(wp) :: vel(3), v(3), up(3)
 real(wp), allocatable :: vel_w(:,:)
 real(wp) :: unorm, alpha, alpha_2d, alpha_avg
 real(wp) , allocatable :: alpha_avg_v(:) , alpha_avg_new_v(:)
 real(wp) :: cl
 real(wp), allocatable :: aero_coeff(:)
 real(wp), allocatable :: dou_temp(:)

 ! mach and reynolds number for each el
 real(wp) :: mach , reynolds
 ! arrays used for force projection
 real(wp) , allocatable :: a_v(:)   ! size(elems_ll)
 real(wp) , allocatable :: c_m(:,:) ! size(elems_ll) , 3
 real(wp) , allocatable :: u_v(:)   ! size(elems_ll)
 real(wp) , allocatable :: ui_v(:,:) ! size(elems_ll)
 real(wp) , allocatable :: dcl_v(:) ! size(elems_ll)
 real(wp) , allocatable :: f_chord_v(:)

 real(wp) , allocatable :: M(:,:)

 !> sim_param
 ! fixed point algorithm for ll
 real(wp) :: fp_tol , fp_damp , diff
 integer  :: fp_maxIter
 ! stall regularisation: params read as inputs
 logical  :: stall_regularisation
 real(wp):: al_stall 
 integer :: i_do , i , nn_stall , n_stall
 integer :: n_iter_reg
 ! load computation
 logical :: load_avl , adaptive_reg
 real(wp) :: e_l(3) , e_d(3)
 real(wp), allocatable :: Gamma_old(:)

 real(wp) :: max_mag_ll

 type(t_liftlin), pointer :: el

 integer, intent(in) :: it

 real(wp) , allocatable :: diff_v(:)
 character(len=4) :: msg_4
 character(len=max_char_len) :: msg
 character(len=*), parameter :: this_sub_name = 'solve_liftlin_piszkin'

!! debug ---
!do i_l = 1 , size(elems_ll)
!  if ( ( allocated(elems_ll(i_l)%p%neigh  ) ) .and. &
!       (associated(elems_ll(i_l)%p        ) ) ) then 
!    do i = 1 , size(elems_ll(i_l)%p%neigh)
!      if ( associated( elems_ll(i_l)%p%neigh(i)%p ) ) then
!        write(*,*) elems_ll(i_l)%p%neigh(i)%p%id
!      end if
!    end do
!  end if
!end do
!write(*,*) ' stop in mod_liftlin.f90, around l. 370' ; stop
!! debug ---

  allocate( diff_v(sim_param%llMaxIter+1) ) ; diff_v = 0.0_wp

  ! params of the fixed point iterations
  fp_tol     = sim_param%llTol
  fp_damp    = sim_param%llDamp
  fp_maxIter = sim_param%llMaxIter
  ! params for stall regularisation
  stall_regularisation = sim_param%llStallRegularisation
  n_stall    = sim_param%llStallRegularisationNelems
  n_iter_reg = sim_param%llStallRegularisationNiters
  al_stall   = sim_param%llStallRegularisationAlphaStall * pi / 180.0_wp
  ! param for load computation
  load_avl   = sim_param%llLoadsAVL
  ! adaptive ll regularisation
  adaptive_reg = sim_param%llArtificialViscosityAdaptive 

  uinf = sim_param%u_inf

  !> allocate and fill Gamma_old array of the ll intensity at previous dt
  allocate(Gamma_old(size(elems_ll)))
  do i_l = 1 , size(elems_ll)
    Gamma_old(i_l) = elems_ll(i_l)%p%mag
  end do

  !> allocate temporary arrays
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

  ! allocate array containing aoa, aero coeffs and relative velocity
  allocate(  a_v( size(elems_ll)  )) ;   a_v = 0.0_wp
  allocate(  c_m( size(elems_ll),3)) ;   c_m = 0.0_wp
  allocate(  u_v( size(elems_ll)  )) ;   u_v = 0.0_wp
  allocate( ui_v( size(elems_ll),3)) ;  ui_v = 0.0_wp
  allocate(dcl_v( size(elems_ll)  )) ; dcl_v = 0.0_wp
  allocate(alpha_avg_v(    size(elems_ll))) ; alpha_avg_v     = 0.0_wp
  allocate(alpha_avg_new_v(size(elems_ll))) ; alpha_avg_new_v = 0.0_wp
  allocate(al_kernel_out(size(elems_ll),size(elems_ll)))
  al_kernel_out = al_kernel

  ! Remove the "out-of-plane" component of the relative velocity:
  ! 2d-velocity to enter the airfoil look-up-tables
  ! IS THIS LOOP USED? (u_v) seems to be overwritten few lines down)
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

  ! Initial condition on the induced angle ( from previous time steps )
  if ( it .ne. 0 ) then
    do i_l = 1 , size(elems_ll)
      if  ( abs(elems_ll(i_l)%p%alpha) .lt. al_stall * 180.0_wp/pi ) then
        a_v(i_l) = elems_ll(i_l)%p%alpha * pi/180.0_wp
      else
        a_v(i_l) = elems_ll(i_l)%p%alpha_isolated * pi/180.0_wp
      end if
    end do
  else
    do i_l = 1 , size(elems_ll)
      a_v(i_l) = elems_ll(i_l)%p%alpha_isolated* pi/180.0_wp
    end do
  end if

  ! ===== Iterative loop =====
  do ic = 1, fp_maxIter
    diff = 0.0_wp             ! max diff ("norm \infty")

    ! === Average value of AOA (regularisation?) ===
    alpha_avg_v = matmul( al_kernel_out , a_v )

!   ! debug ---
!   if ( ic .eq. 1 ) then
!     write(*,*) '%> First iteration : i_l , al[deg] , avg al[deg]' 
!     do i_l = 1 , size(elems_ll)
!       write(*,*) i_l , a_v(i_l)*180.0_wp/pi , alpha_avg_v(i_l) 
!     end do
!   end if
!   ! debug ---
 
    ! === Update LL intensity ===
!$omp parallel do private(i_l, el, j, v, vel, up, unorm, alpha, alpha_avg, alpha_2d, mach, &
!$omp& reynolds, aero_coeff, cl) schedule(dynamic,4)
    do i_l = 1,size(elems_ll)

      !select type(el => elems_ll(i_l)%p) ; type is(t_liftlin)
      el => elems_ll(i_l)%p

        alpha = a_v(i_l) * 180.0_wp/pi

        ! === Average value of AOA (regularisation?) ===
        alpha_avg = alpha_avg_v(i_l) * 180.0_wp/pi
        
        ! overall relative velocity computed in the centre of the ll elem
        vel = ui_v(i_l,:) + uinf - el%ub + vel_w(:,i_l)
        ! "effective" velocity = proj. of vel in the n-t plane
        !!! up =  el%nor*sum(el%nor*vel) + el%tang_cen*sum(el%tang_cen*vel)
        !!! u_v(i_l) = norm2(up) 
        unorm = u_v(i_l)      ! velocity w/o induced velocity
       
        ! compute local Reynolds and Mach numbers for the section
        ! needed to enter the LUT (.c81) of aerodynamic loads (2d airfoil)
        mach     = unorm / sim_param%a_inf      
        reynolds = sim_param%rho_inf * unorm * & 
                   el%chord / sim_param%mu_inf    

        ! Read the aero coeff from .c81 tables
        call interp_aero_coeff ( airfoil_data,  el%csi_cen, el%i_airfoil , &
                                   (/alpha_avg, mach, reynolds/) , sim_param , &
                                                      aero_coeff , dcl_v(i_l) )
        cl = aero_coeff(1)   ! cl needed for the iterative process

        ! Compute the "equivalent" intensity of the vortex line 
        dou_temp(i_l) = - 0.5_wp * unorm * cl * el%chord

        c_m(i_l,:) = aero_coeff

        el%vel_outplane = sum(el%bnorm_cen*vel)
      !end select

    enddo  ! i_l
!$omp end parallel do

    ! === Update ll intensity ===
    do i_l = 1,size(elems_ll)
      elems_ll(i_l)%p%mag = dou_temp(i_l)
    enddo

    ! === Update AOA and velocity ===
    diff = 0.0_wp
!$omp parallel do private(i_l, el, j, v, vel, up, unorm, alpha, alpha_2d, mach, &
!$omp& reynolds, aero_coeff, cl) schedule(dynamic,4)
    do i_l = 1,size(elems_ll)

      ! compute velocity
      vel = 0.0_wp
      do j = 1,size(elems_ll)
        call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
        vel = vel + v
      enddo
      ui_v(i_l,:) = vel / ( 4.0_wp * pi )

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
        ! =========================================================

        ! === Damped update of the aoa ===
        a_v(i_l) = a_v(i_l) + 1.0_wp/fp_damp * ( alpha * pi/180.0_wp - a_v(i_l) )

! !$omp atomic
!         diff = max( diff, abs(alpha*pi/180.0_wp-a_v(i_l)) )
! !$omp end atomic

    enddo  ! i_l
!$omp end parallel do

    ! === Average value of AOA (regularisation?) ===
    alpha_avg_new_v = matmul( al_kernel_out , a_v )
    diff = maxval( abs( alpha_avg_new_v - alpha_avg_v ) )

    diff_v(ic) = diff

    ! === Update kernel ( for adaptive algorithm ) ===
    ! call build_ll_array_optim( elems_ll , M , al_kernel )
    if ( adaptive_reg ) then
      call update_kernel( elems_ll , al_kernel , dcl_v , al_kernel_out , f_chord_v ) 
    end if

    !> Stopping criterion
    if ( diff .le. fp_tol * pi/180.0_wp ) then
      exit ! convergence
    end if
!   if ( diff .le. 0.1_wp * pi/180.0_wp ) exit ! convergence

  enddo !solver iterations

  ! debug ---
  if ( adaptive_reg ) then
    do i_l = 1 , size(elems_ll)
      write(*,*) f_chord_v(i_l)
    end do
  end if
  ! debug ---

! ! debug ---
! write(*,*) '%> Last iteration : i_l , al[deg] , avg al[deg]' 
! do i_l = 1 , size(elems_ll)
!   write(*,*) i_l , a_v(i_l)*180.0_wp/pi , alpha_avg_new_v(i_l) 
! end do
! ! debug ---

  write(msg_4,'(I4.4)') it
  open(unit=21,file=trim(sim_param%basename_debug)//'_conv_'//msg_4//'.dat')
  do i_l = 1 , ic-1
    write(21,*) diff_v(i_l)
  end do
  close(21)

  ! === Update el % alpha === ( here or updated values below, as in Piszkin? )
  do i_l = 1 , size(elems_ll)
    elems_ll(i_l)%p%alpha = a_v(i_l) * 180.0_wp/pi
  end do

  ! === Update dGamma_dt field ===
  do i_l = 1,size(elems_ll)
    elems_ll(i_l)%p%dGamma_dt = ( elems_ll(i_l)%p%mag - Gamma_old(i_l) ) / &
                                  sim_param%dt 
  end do

  ! === Loads computation ===
  ! compute LL sectional loads from singularity intensities. 
  do i_l = 1,size(elems_ll)

   !select type(el => elems_ll(i_l)%p)
   !type is(t_liftlin)
   el => elems_ll(i_l)%p
    ! avg delta_p = \vec{F}.\vec{n} / A = ( L*cos(al)+D*sin(al) ) / A
    !> steady pressure contribution ( overwritten below )
    el%pres   = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
               ( c_m(i_l,1) * cos(a_v(i_l)) +  c_m(i_l,2) * sin(a_v(i_l)) )

    if ( load_avl ) then

      ! === Kutta-Joukowski theorem ~ AVL ===
      !> Inviscid contribution ~ AVL
      call el % compute_dforce_jukowski ( elems_tot , elems_wake )

      !> Add viscous contribution
      el % dforce = el % dforce + &
           0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * el%area * &
           c_m(i_l,2) * ( sin(el%al_ctr_pt) * el%nor      +  &
                          cos(el%al_ctr_pt) * el%tang_cen )

      ! === Update AOAs and aerodynamic coefficients ===
      !> unit vector in the direction of the relative velocity (_d for drag)
      !  and in the direction of the lift (_l for lift)
      e_l = el%nor * cos( el%al_ctr_pt ) - el%tang_cen * sin( el%al_ctr_pt )
      e_d = el%nor * sin( el%al_ctr_pt ) + el%tang_cen * cos( el%al_ctr_pt )
      e_l = e_l / norm2(e_l)
      e_d = e_d / norm2(e_d)

      !> Unsteady contribution
      el%dforce = el%dforce &
                  - sim_param%rho_inf * el%area * el%dGamma_dt &
                  * e_l ! lift direction

      !> Update aerodynamic coefficients and AOA, and pressure
      c_m(i_l,1) = sum( el % dforce * e_l ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      c_m(i_l,2) = sum( el % dforce * e_d ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      a_v(i_l) = el % al_ctr_pt

      !> overwrite pressure field to take into account unsteady contributions
      el%pres = sum( el%dforce * el%nor ) / el%area

    else

      ! === elementary force = p*n + tangential contribution from L,D ===
      el%dforce = ( el%nor * el%pres + &
                    el%tang_cen * &
                    0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * ( &
                   -c_m(i_l,1) * sin(a_v(i_l)) + c_m(i_l,2) * cos(a_v(i_l)) &
                   ) ) * el%area

      ! === Update AOAs and aerodynamic coefficients ===
      !> unit vector in the direction of the relative velocity (_d for drag)
      !  and in the direction of the lift (_l for lift)
      e_l = el%nor * cos( a_v(i_l) ) - el%tang_cen * sin( a_v(i_l) )
      e_d = el%nor * sin( a_v(i_l) ) + el%tang_cen * cos( a_v(i_l) )
      e_l = e_l / norm2(e_l)
      e_d = e_d / norm2(e_d)

      !> Unsteady contribution
      el%dforce = el%dforce &
                  - sim_param%rho_inf * el%area * el%dGamma_dt &
                  * e_l ! lift direction

      !> Update aerodynamic coefficients and AOA, and pressure
      c_m(i_l,1) = sum( el % dforce * e_l ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      c_m(i_l,2) = sum( el % dforce * e_d ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )

      !> overwrite pressure field to take into account unsteady contributions
      el%pres = sum( el%dforce * el%nor ) / el%area

    end if

    ! elementary moment = 0.5 * rho * v^2 * A * c * cm, 
    ! - around bnorm_cen (always? TODO: check)
    ! - referred to the ref.point of the elem,
    !   ( here, cen of the elem = cen of the liftlin (for liftlin elems) )
    el%dmom = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
                   el%chord * el%area * c_m(i_l,3)

!   el%alpha = a_v(i_l) * 180_wp/pi
    el%vel_2d = u_v(i_l)
    el%aero_coeff = c_m(i_l,:)

   !end select
  end do

  if(sim_param%debug_level .ge. 3) then
    write(msg,*) 'iterations: ',ic; call printout(trim(msg))
    write(msg,*) 'diff',diff*180.0_wp/pi; call printout(trim(msg))
  endif

  ! useful arrays ---
  deallocate(dou_temp, vel_w, Gamma_old)
  deallocate(a_v,c_m,u_v, alpha_avg_v)
  

end subroutine solve_liftlin_piszkin

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

 ! mach and reynolds number for each el
 real(wp) :: mach , reynolds
 ! arrays used for force projection
 real(wp) , allocatable :: a_v(:)   ! size(elems_ll)
 real(wp) , allocatable :: c_m(:,:) ! size(elems_ll) , 3
 real(wp) , allocatable :: u_v(:)   ! size(elems_ll)

 !> sim_param
 ! fixed point algorithm for ll
 real(wp) :: fp_tol , fp_damp , diff
 integer  :: fp_maxIter
 ! stall regularisation: params read as inputs
 logical  :: stall_regularisation
 real(wp):: al_stall 
 integer :: i_do , i , nn_stall , n_stall
 integer :: n_iter_reg
 ! load computation
 logical :: load_avl
 real(wp) :: e_l(3) , e_d(3)
 real(wp), allocatable :: Gamma_old(:)

 real(wp) :: max_mag_ll

 type(t_liftlin), pointer :: el

 integer, intent(in) :: it

 real(wp) , allocatable :: diff_v(:)
 character(len=4) :: msg_4
 character(len=max_char_len) :: msg
 character(len=*), parameter :: this_sub_name = 'solve_liftlin'

!! debug ---
!do i_l = 1 , size(elems_ll)
!  if ( ( allocated(elems_ll(i_l)%p%neigh  ) ) .and. &
!       (associated(elems_ll(i_l)%p        ) ) ) then 
!    do i = 1 , size(elems_ll(i_l)%p%neigh)
!      if ( associated( elems_ll(i_l)%p%neigh(i)%p ) ) then
!        write(*,*) elems_ll(i_l)%p%neigh(i)%p%id
!      end if
!    end do
!  end if
!end do
!write(*,*) ' stop in mod_liftlin.f90, around l. 370' ; stop
!! debug ---

  allocate( diff_v(sim_param%llMaxIter+1) ) ; diff_v = 0.0_wp

  ! params of the fixed point iterations
  fp_tol     = sim_param%llTol
  fp_damp    = sim_param%llDamp
  fp_maxIter = sim_param%llMaxIter
  ! params for stall regularisation
  stall_regularisation = sim_param%llStallRegularisation
  n_stall    = sim_param%llStallRegularisationNelems
  n_iter_reg = sim_param%llStallRegularisationNiters
  al_stall   = sim_param%llStallRegularisationAlphaStall * pi / 180.0_wp
  ! param for load computation
  load_avl   = sim_param%llLoadsAVL

  uinf = sim_param%u_inf

  !> allocate and fill Gamma_old array of the ll intensity at previous dt
  allocate(Gamma_old(size(elems_ll)))
  do i_l = 1 , size(elems_ll)
    Gamma_old(i_l) = elems_ll(i_l)%p%mag
  end do

  !> allocate temporary arrays
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

  ! allocate array containing aoa, aero coeffs and relative velocity
  allocate( a_v( size(elems_ll)  )) ;   a_v = 0.0_wp
  allocate( c_m( size(elems_ll),3)) ;   c_m = 0.0_wp
  allocate( u_v( size(elems_ll)  )) ;   u_v = 0.0_wp

  ! Remove the "out-of-plane" component of the relative velocity:
  ! 2d-velocity to enter the airfoil look-up-tables
  ! IS THIS LOOP USED? (u_v) seems to be overwritten few lines down)
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
        ! =========================================================

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
 
    diff_v(ic) = diff/max_mag_ll

    if ( diff/max_mag_ll .le. fp_tol ) exit ! convergence

  enddo !solver iterations

  write(msg_4,'(I4.4)') it
  open(unit=21,file=trim(sim_param%basename_debug)//'_conv_'//msg_4//'.dat')
  do i_l = 1 , ic-1
    write(21,*) diff_v(i_l)
  end do
  close(21)

  ! === Update dGamma_dt field ===
  do i_l = 1,size(elems_ll)
    elems_ll(i_l)%p%dGamma_dt = ( elems_ll(i_l)%p%mag - Gamma_old(i_l) ) / &
                                  sim_param%dt 
  end do

  ! === Loads computation ===
  ! compute LL sectional loads from singularity intensities. 
  do i_l = 1,size(elems_ll)

   !select type(el => elems_ll(i_l)%p)
   !type is(t_liftlin)
   el => elems_ll(i_l)%p
    ! avg delta_p = \vec{F}.\vec{n} / A = ( L*cos(al)+D*sin(al) ) / A
    !> steady pressure contribution ( overwritten below )
    el%pres   = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
               ( c_m(i_l,1) * cos(a_v(i_l)) +  c_m(i_l,2) * sin(a_v(i_l)) )

    if ( load_avl ) then

      ! === Kutta-Joukowski theorem ~ AVL ===
      !> Inviscid contribution ~ AVL
      call el % compute_dforce_jukowski ( elems_tot , elems_wake )

      !> Add viscous contribution
      el % dforce = el % dforce + &
           0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * el%area * &
           c_m(i_l,2) * ( sin(el%al_ctr_pt) * el%nor      +  &
                          cos(el%al_ctr_pt) * el%tang_cen )

      ! === Update AOAs and aerodynamic coefficients ===
      !> unit vector in the direction of the relative velocity (_d for drag)
      !  and in the direction of the lift (_l for lift)
      e_l = el%nor * cos( el%al_ctr_pt ) - el%tang_cen * sin( el%al_ctr_pt )
      e_d = el%nor * sin( el%al_ctr_pt ) + el%tang_cen * cos( el%al_ctr_pt )
      e_l = e_l / norm2(e_l)
      e_d = e_d / norm2(e_d)

      !> Unsteady contribution
      el%dforce = el%dforce &
                  - sim_param%rho_inf * el%area * el%dGamma_dt &
                  * e_l ! lift direction

      !> Update aerodynamic coefficients and AOA, and pressure
      c_m(i_l,1) = sum( el % dforce * e_l ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      c_m(i_l,2) = sum( el % dforce * e_d ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      a_v(i_l) = el % al_ctr_pt

      !> overwrite pressure field to take into account unsteady contributions
      el%pres = sum( el%dforce * el%nor ) / el%area

    else

      ! === elementary force = p*n + tangential contribution from L,D ===
      el%dforce = ( el%nor * el%pres + &
                    el%tang_cen * &
                    0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * ( &
                   -c_m(i_l,1) * sin(a_v(i_l)) + c_m(i_l,2) * cos(a_v(i_l)) &
                   ) ) * el%area

      ! === Update AOAs and aerodynamic coefficients ===
      !> unit vector in the direction of the relative velocity (_d for drag)
      !  and in the direction of the lift (_l for lift)
      e_l = el%nor * cos( a_v(i_l) ) - el%tang_cen * sin( a_v(i_l) )
      e_d = el%nor * sin( a_v(i_l) ) + el%tang_cen * cos( a_v(i_l) )
      e_l = e_l / norm2(e_l)
      e_d = e_d / norm2(e_d)

      !> Unsteady contribution
      el%dforce = el%dforce &
                  - sim_param%rho_inf * el%area * el%dGamma_dt &
                  * e_l ! lift direction

      !> Update aerodynamic coefficients and AOA, and pressure
      c_m(i_l,1) = sum( el % dforce * e_l ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )
      c_m(i_l,2) = sum( el % dforce * e_d ) / &
        ( 0.5_wp * sim_param%rho_inf * u_v(i_l) ** 2.0_wp * el%area )

      !> overwrite pressure field to take into account unsteady contributions
      el%pres = sum( el%dforce * el%nor ) / el%area

    end if

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
  deallocate(dou_temp, vel_w, Gamma_old)
  deallocate(a_v,c_m,u_v)
  

end subroutine solve_liftlin

!----------------------------------------------------------------------
!> Compute the inviscid contribution of the elementary force on the on
! the actual element using Kutta-Jukowski theorem as implemented in AVL,
! using the velocity computed at the ctr_pt on the LL, in the subroutine
! get_vel_ctr_cp_liftlin. The viscous contribution is performed in
! the subroutine solve_liftlin
! 
subroutine compute_dforce_jukowski_liftlin(this, elems, wake_elems) 
 class(t_liftlin), intent(inout) :: this
 type(t_pot_elem_p),intent(in):: elems(:)
 type(t_pot_elem_p),intent(in):: wake_elems(:)

 real(wp) :: vel_w(3), gam(3)

 call this % get_vel_ctr_pt( elems , wake_elems )

 gam = cross ( this % vel_ctr_pt, this % edge_vec(:,1) )

 this%dforce = sim_param%rho_inf * gam * this%mag

! ! debug ---
! write(*,*) ' this % vel_ctr_pt : ' , this % vel_ctr_pt
! write(*,*) ' sim_param%rho_inf : ' , sim_param%rho_inf 
! write(*,*) ' gam               : ' , gam               
! write(*,*) ' this%mag          : ' , this%mag          
! ! debug ---

end subroutine compute_dforce_jukowski_liftlin

!----------------------------------------------------------------------
!> Compute the velocity and the "induced" incidence angle at ctr_pt
! located on the LL. These quantities are used in the inviscid load 
! computation, using AVL expression (~ VL elements)
!
subroutine get_vel_ctr_pt_liftlin(this, elems, wake_elems)
 class(t_liftlin), intent(inout) :: this
 type(t_pot_elem_p),intent(in):: elems(:)
 type(t_pot_elem_p),intent(in):: wake_elems(:)
 
 real(wp) :: v(3),x0(3)
 integer :: j

 ! Initialisation to zero
 this%vel_ctr_pt = 0.0_wp

 ! Control point at 1/4-fraction of the chord
 x0 = this%cen - this%tang_cen * this%chord / 2.0_wp

 !=== Compute the velocity from all the elements ===
 do j = 1,size(wake_elems)  ! wake panels
 
   call wake_elems(j)%p%compute_vel(x0,sim_param%u_inf,v)
   this%vel_ctr_pt = this%vel_ctr_pt + v
 
 enddo

 do j = 1,size(elems) ! body elements
 
   call elems(j)%p%compute_vel(x0,sim_param%u_inf,v)
   this%vel_ctr_pt = this%vel_ctr_pt + v
 
 enddo

 this%vel_ctr_pt = this%vel_ctr_pt/(4.0_wp*pi) &
               + sim_param%u_inf + this%uvort - this%ub

 this%al_ctr_pt = atan2( sum(this%vel_ctr_pt * this%nor     ) , & 
                         sum(this%vel_ctr_pt * this%tang_cen) )

end subroutine get_vel_ctr_pt_liftlin

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
! this%area = 0.5_wp * norm2(nor) / 0.75_wp 
  this%area = 0.5_wp * norm2(nor)
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
! this%chord = this%chord / 0.75_wp

  ! === Piszkin, Lewinski (1976) LL model for swept wings ===
  ! - cos_lambda = cos(lambda) , where lambda is the local sweep angle
  ! - %cen is overwritten and set = to %cen, because the fmm routines
  !   compute velocity on the %cen of the elems, so far.
  ! - %d_2pi_coslambda = x_CP - x_{1/4*c} * 2*pi * cos(lambda)
  !
  !> sweep angle
  cos_lambda  = norm2( cross( this%tang_cen , -this%edge_uni(:,1) ) )
 
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
