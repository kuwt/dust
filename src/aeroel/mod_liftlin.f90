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
  t_sim_param

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

public :: t_liftlin, update_liftlin, solve_liftlin


!----------------------------------------------------------------------

type, extends(c_expl_elem) :: t_liftlin
  real(wp), allocatable :: tang_cen(:)
  real(wp), allocatable :: bnorm_cen(:)
  real(wp)              :: csi_cen
  integer               :: i_airfoil(2)
  real(wp)              :: chord
  real(wp)              :: ctr_pt(3)
  real(wp)              :: nor_zeroLift(3)
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
subroutine compute_pres_liftlin (this, R_g, sim_param)
 class(t_liftlin) , intent(inout) :: this
 real(wp)         , intent(in)    :: R_g(3,3)
 type(t_sim_param), intent(in)    :: sim_param
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
subroutine compute_dforce_liftlin (this, sim_param)
 class(t_liftlin), intent(inout) :: this
 !type(t_elem_p), intent(in) :: elems(:)
 type(t_sim_param), intent(in) :: sim_param

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
 type(t_expl_elem_p), intent(inout) :: elems_ll(:)
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
                         sim_param, airfoil_data, it)
 type(t_expl_elem_p), intent(inout) :: elems_ll(:)
 type(t_pot_elem_p),  intent(in)    :: elems_tot(:)
 type(t_impl_elem_p), intent(in)    :: elems_impl(:)
 type(t_expl_elem_p), intent(in)    :: elems_ad(:)
 type(t_pot_elem_p),  intent(in)    :: elems_wake(:)
 type(t_vort_elem_p), intent(in)    :: elems_vort(:)
 type(t_sim_param),   intent(in)    :: sim_param
 type(t_aero_tab),    intent(in)    :: airfoil_data(:)
 real(wp) :: uinf(3)
 integer  :: i_l, j, ic
 real(wp) :: vel(3), v(3), up(3)
 real(wp), allocatable :: vel_w(:,:)
 real(wp) :: unorm, alpha
 real(wp) :: cl
 real(wp), allocatable :: aero_coeff(:)
 real(wp), allocatable :: dou_temp(:)

 ! mach and reynolds number for each el
 real(wp) :: mach , reynolds
 ! arrays used for force projection
 real(wp) , allocatable :: a_v(:)   ! size(elems_ll)
 real(wp) , allocatable :: c_m(:,:) ! size(elems_ll) , 3
 real(wp) , allocatable :: u_v(:)   ! size(elems_ll)

 ! fixed point algorithm for ll ----
 real(wp) :: fp_tol , fp_damp , diff
 integer  :: fp_maxIter
 logical  :: stall_regularisation
 integer              :: n_inner_stall
 real(wp) , parameter :: al_stall = 15.0_wp*pi/180.0_wp ! hardcoded

 real(wp) :: max_mag_ll

! debug ----
 integer, intent(in) :: it
 character(len=4) :: it_str
 real(wp) , allocatable :: re_v(:) , ma_v(:)

 real(wp) , allocatable :: vel_ctr(:,:)
 real(wp) , allocatable :: chord_v(:) , cl_v(:)
 real(wp) , allocatable :: el_csi_cen_v(:)
 integer  , allocatable :: el_i_airfoil_m(:,:) 
! debug ----

! newton ---
 real(wp) , allocatable :: ll_mag(:)
 real(wp) , allocatable :: Amat(:,:) , Amat_newton(:,:)
 real(wp) , allocatable :: bvec(:) , dclvec(:) , fvec(:) , clvec(:)
 real(wp) , allocatable :: alvec(:) , d_al(:)
 integer :: nll , ii , ik
 real(wp) , parameter :: newton_tol = 1e-6_wp
 integer  , parameter :: newton_maxit = 100
 integer :: newton_it
 real(wp) :: dclda
 integer, allocatable :: ipiv(:)
 integer :: info
 real(wp) :: res , res_old
 real(wp) :: al0
! newton ---

 
 character(len=max_char_len) :: msg
 character(len=*), parameter :: this_sub_name = 'solve_liftlin'

  ! parameters of the fixed point iterations
  fp_tol     = sim_param%llTol
  fp_damp    = sim_param%llDamp
  fp_maxIter = sim_param%llMaxIter
  stall_regularisation = sim_param%llStallRegularisation
! n_inner_stall        = sim_param%llStallRegularisationNelems
  n_inner_stall = 1

  uinf = sim_param%u_inf

  nll = size(elems_ll)

  allocate(dou_temp(size(elems_ll))) ; dou_temp = 0.0_wp

  !Compute the velocity from all the elements except for liftling elems
  ! and store it outside the loop, since it is constant
  ! === %cen === velocity on the centres
  allocate(vel_w(3,size(elems_ll))) ; vel_w = 0.0_wp 
  do i_l = 1,size(elems_ll)
    do j = 1,size(elems_impl) ! body panels: liftlin, vortlat 
      call elems_impl(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
      vel_w(:,i_l) = vel_w(:,i_l) + v
    enddo
    do j = 1,size(elems_ad) ! actuator disks 
      call elems_ad(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
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
  enddo

  vel_w = vel_w/(4.0_wp*pi)

  ! === %ctr_pt === velocity on the control points
  allocate( vel_ctr(3,size(elems_ll)) ) ; vel_ctr = 0.0_wp
  do i_l = 1,size(elems_ll)
    select type( el => elems_ll(i_l)%p ) ; type is (t_liftlin) 
      do j = 1,size(elems_impl) ! body panels: liftlin, vortlat 
          call elems_impl(j)%p%compute_vel(el%ctr_pt,uinf,v)
          vel_ctr(:,i_l) = vel_ctr(:,i_l) + v
      enddo
      do j = 1,size(elems_ad) ! actuator disks 
          call elems_ad(j)%p%compute_vel(el%ctr_pt,uinf,v)
          vel_ctr(:,i_l) = vel_ctr(:,i_l) + v
      enddo
      do j = 1,size(elems_wake) ! wake panels
          call elems_wake(j)%p%compute_vel(el%ctr_pt,uinf,v)
          vel_ctr(:,i_l) = vel_ctr(:,i_l) + v
      enddo
      do j = 1,size(elems_vort) ! wake vort
          call elems_vort(j)%p%compute_vel(el%ctr_pt,uinf,v)
          vel_ctr(:,i_l) = vel_ctr(:,i_l) + v
      enddo
    end select
  enddo

  vel_ctr = vel_ctr/(4.0_wp*pi)

  ! ------

  allocate( a_v( size(elems_ll)  )) ;   a_v = 0.0_wp
  allocate( c_m( size(elems_ll),3)) ;   c_m = 0.0_wp
  allocate( u_v( size(elems_ll)  )) ;   u_v = 0.0_wp
  allocate(re_v( size(elems_ll)  )) ;  re_v = 0.0_wp
  allocate(ma_v( size(elems_ll)  )) ;  ma_v = 0.0_wp

  ! debug ----
  allocate(chord_v( size(elems_ll)  )) ; chord_v = 0.0_wp
  allocate(   cl_v( size(elems_ll)  )) ;    cl_v = 0.0_wp
  allocate(el_i_airfoil_m( size(elems_ll),2)) ; el_i_airfoil_m = 0.0_wp
  allocate(el_csi_cen_v  ( size(elems_ll)  )) ; el_csi_cen_v   = 0.0_wp
  ! debug ----

  ! Remove the "out-of-plane" component of the relative velocity:
  ! 2d-velocity to enter the airfoil look-up-tables
  do i_l=1,size(elems_ll)
   select type(el => elems_ll(i_l)%p)
   type is(t_liftlin)
     u_v(i_l) = norm2((uinf-el%ub) - &
         el%bnorm_cen*sum(el%bnorm_cen*(uinf-el%ub))) 
   end select
  end do


  ! allocate and initialised to zero other vectors used in newton methods
  allocate(fvec(  nll)) ; fvec   = 0.0_wp
  allocate(dclvec(nll)) ; dclvec = 0.0_wp
  allocate( clvec(nll)) ;  clvec = 0.0_wp
  allocate( alvec(nll)) ;  alvec = 0.0_wp
  ! newton ---

  !Calculate the induced velocity on the airfoil
  do ic = 1, fp_maxIter       !TODO: Refine this iterative process 
    diff = 0.0_wp             ! max diff ("norm \infty")
    max_mag_ll = 0.0_wp
    do i_l = 1,size(elems_ll)

      ! compute velocity
      vel = 0.0_wp
      do j = 1,size(elems_ll)
! ! debug ----
!         write(*,'(A,4(F12.4))') ' elems(j)%mag , elems(l)%cen ' , &
!                      elems_ll(j)%p%mag , elems_ll(i_l)%p%cen
! ! debug ----
        call elems_ll(j)%p%compute_vel(elems_ll(i_l)%p%cen,uinf,v)
        vel = vel + v
! ! debug ----
!         write(*,'(A,3(F12.4))') '        v : ' , v
! ! debug ----
      enddo
! ! debug ----
!       write(*,*) '      vel : ' , vel 
! ! debug ----

      select type(el => elems_ll(i_l)%p)
      type is(t_liftlin)
        vel = vel/(4.0_wp*pi) + uinf - el%ub +vel_w(:,i_l)
        !vel =                 + uinf - el%ub +vel_w(:,i_l)
        !vel = uinf - el%ub
!       up = vel-el%bnorm_cen*sum(el%bnorm_cen*vel)
        up = uinf-el%ub+vel_w(:,i_l) &
            -el%bnorm_cen*sum(el%bnorm_cen*(uinf-el%ub+vel_w(:,i_l)))
        u_v(i_l) = norm2(up) ! overwrite !!!
        up = vel-el%bnorm_cen*sum(el%bnorm_cen*vel)

        ! Compute reference velocity (includes the motion of the body)
        ! to use in LUT (.c81) of aerodynamic loads (2d airfoil)
        !TODO: test these two approximations
        unorm = u_v(i_l)      ! velocity w/o induced velocity
      ! unorm = norm2(up)     ! full velocity
       
        ! Angle of incidence (full velocity)
        alpha = atan2(sum(up*el%nor), sum(up*el%tang_cen))
        alpha = alpha * 180.0_wp/pi  ! .c81 tables defined with angles in [deg]

        ! compute local Reynolds and Mach numbers for the section
        ! needed to enter the LUT (.c81) of aerodynamic loads (2d airfoil)
        mach     = unorm / sim_param%a_inf      
        reynolds = sim_param%rho_inf * unorm * & 
                   el%chord / sim_param%mu_inf    

        ! debug ----
        el_i_airfoil_m(i_l,:) = el%i_airfoil
        el_csi_cen_v(  i_l  ) = el%csi_cen
        ! debug ----

        ! read the aero coeff from .c81 tables
        call interp_aero_coeff ( airfoil_data,  el%csi_cen, el%i_airfoil , &
                                   (/alpha, mach, reynolds/) , sim_param , &
                                                              aero_coeff )
        cl = aero_coeff(1)   ! cl needed for the iterative process
! newton cleaning
!       call interp_aero_coeff ( airfoil_data,  el%csi_cen, el%i_airfoil, &
!                                 (/alpha, mach, reynolds/) , sim_param , &
!                                              aero_coeff , dclda , al0 )
! newton cleaning

! newton cleaning
!       ! >>> update/define the zero lift normal versor for liftlin elems
!       al0 = al0 * pi/180.0_wp
!       el%nor_zeroLift = el%nor * cos(al0) - el%tang_cen * sin(al0) 
! newton cleaning

!       ! debug ----
!       write(*,*) ' cl , cd , cm , dclda ' 
!       write(*,*) aero_coeff , dclda
!       ! debug ----
       
        ! Compute the "equivalent" intensity of the vortex line 
        dou_temp(i_l) = - 0.5_wp * unorm * cl * el%chord
        diff = max(diff,abs(elems_ll(i_l)%p%mag-dou_temp(i_l)))

        ! debug ----
        chord_v(i_l) = el%chord 
        cl_v(i_l) = cl
        ! debug ----

      end select

      re_v(i_l) = reynolds
      ma_v(i_l) = mach
      c_m(i_l,:) = aero_coeff
      a_v(i_l)   = alpha * pi/180.0_wp ! [rad]

      ! --- to be used for the first newton iteration ---
      dclvec(i_l) = dclda   ! [deg]
      ! clvec, alvec set before the Newton loop

    enddo  ! i_l

    ! --- avoid stall on an elem between 2 elems without stall ---
    ! TODO:
    ! - read neighbouring elements and not the prev and next elem in ll numbering <<<<
    ! - decide if using .and. or .or. in the if statement (or user input ?)       <<<<
    ! - use n_inner_stall to check more than one inner elems 
    ! - read al_stall from tables (?) or leave it as an input from the user
    ! al_stall = ... \hardcoded as a parameter
    ! rough ----
    if ( stall_regularisation ) then 
      do i_l = 2 , size(elems_ll) - 1
        if ( a_v(i_l) .gt. al_stall ) then
          if ( ( a_v(i_l-1) .lt. al_stall ) .and. &
               ( a_v(i_l+1) .lt. al_stall ) ) then
             a_v(     i_l) = 0.5_wp * ( a_v(     i_l-1) + a_v(     i_l+1) ) 
             dou_temp(i_l) = 0.5_wp * ( dou_temp(i_l-1) + dou_temp(i_l+1) ) 
           end if
        elseif ( a_v(i_l) .lt. -al_stall ) then
          if ( ( a_v(i_l-1) .gt. -al_stall ) .and. &
               ( a_v(i_l+1) .gt. -al_stall ) ) then
             a_v(     i_l) = 0.5_wp * ( a_v(     i_l-1) + a_v(     i_l+1) ) 
             dou_temp(i_l) = 0.5_wp * ( dou_temp(i_l-1) + dou_temp(i_l+1) ) 
           end if
        end if
      end do
    end if
    ! rough ---- 

    ! Update ll intensity
    do i_l = 1,size(elems_ll)
!     ! debug ----
!     write(*,*) elems_ll(i_l)%p%mag
!     ! debug ----
      elems_ll(i_l)%p%mag = ( dou_temp(i_l)+ fp_damp*elems_ll(i_l)%p%mag )&
                             /(1.0_wp+fp_damp)
      max_mag_ll = max(max_mag_ll,abs(elems_ll(i_l)%p%mag))
    enddo

!   ! debug ----
!   write(*,*) ' Iter n.' , ic 
!   write(*,'(A)') ' i_l  ,   u_v  ,   a_v  ,  chord ,   cl  , dou_temp : ' 
!   do i_l = 1,size(elems_ll)
!      write(*,'(I4.4,5(F10.4),2(I4.2),F10.4)') &
!          i_l , u_v(i_l) , a_v(i_l) , chord_v(i_l) , cl_v(i_l) , dou_temp(i_l) , &
!          el_i_airfoil_m(i_l,:) , el_csi_cen_v(i_l)
!   end do
!   ! debug ----

    if ( diff/max_mag_ll .le. fp_tol ) exit

  enddo

  ! Loads computation ------------
  do i_l = 1,size(elems_ll)
   select type(el => elems_ll(i_l)%p)
   type is(t_liftlin)
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


   end select
  end do

  if(sim_param%debug_level .ge. 3) then
    write(msg,*) 'iterations: ',ic; call printout(trim(msg))
    write(msg,*) 'diff',diff; call printout(trim(msg))
    write(msg,*) 'diff/max_mag_ll:',diff/max_mag_ll; call printout(trim(msg))
  endif


! newton cleaning
!   ! === Newton's algorithm to assign u_rel.n = 0 ===
!   
!   ! === Amat, bvec === 
!   ! allocate and fill constant vectors and array for newton method
!   allocate(Amat(nll,nll), bvec(nll)) ; Amat = 0.0_wp ; bvec = 0.0_wp
! 
!   ! set all the intensity of the liftling lines to 1.0_wp in order to
!   ! compute the influence of the unitary vortices on the control points 
!   ! needed to build matrix Amat
!   allocate(ll_mag(nll)) 
!   do i_l = 1 , nll
!     ll_mag(i_l) = elems_ll(i_l)%p%mag 
!     elems_ll(i_l)%p%mag = 1.0_wp
!   end do
!  
!   ! fill constant matrix Amat and vector bvec 
!   do ii = 1 , nll
!     select type( el=>elems_ll(ii)%p ) ; type is (t_liftlin)
!       bvec(ii) = sum( ( uinf - el%ub + vel_ctr(:,ii) ) * el%nor_zeroLift )
!       do ik = 1 , nll
!         call elems_ll(ik)%p%compute_vel(el%ctr_pt,uinf,v)
!         Amat(ii,ik) = sum( v / (4.0_wp*pi) * el%nor_zeroLift )
!       end do 
!     end select
!   end do 
! 
! ! ! debug ----
! ! write(it_str,'(I4.4)') it 
! ! open(unit=21, file='./Debug/ll_Amat_'//trim(it_str)//'.dat')
! ! do i_l = 1 , nll
! !   write(21,*) Amat(i_l,:)
! ! end do
! ! close(21)
! ! ! debug ----
! 
!   do ik = 1 , nll
!     select type( el=>elems_ll(ik)%p ) ; type is (t_liftlin)
!       Amat(:,ik) = Amat(:,ik) * (-0.5*u_v(ik)*el%chord)
!     end select
!   end do
!   
!   do i_l = 1,nll  !Initialise ll magnitude to the value at previous timestep
!     elems_ll(i_l)%p%mag = ll_mag(i_l) 
!   end do
!   ! === Amat, bvec === 
! 
! ! ! debug ----
! ! write(*,*) ' Amat : '
! ! do ii = 1 , nll
! !   write(*,*) Amat(ii,:)
! ! end do
! ! ! debug ----
! 
!   clvec = c_m(:,1)                  ! from previous algorithm / initial guess
!   alvec = a_v(:)*180.0_wp/pi        ! from previous algorithm / initial guess
!   fvec = matmul(Amat,clvec) + bvec  ! residue 
! 
!   ! debug ----
!   write(*,*) '    a_v    ,    alvec    ,    clvec    ,    dclvec    ,    fvec    ,    mag ' 
!   do i_l = 1 , nll
!    select type(el => elems_ll(i_l)%p) ; type is (t_liftlin)
!     write(*,'(6E14.4)') a_v(i_l) , alvec(i_l) , clvec(i_l) , dclvec(i_l) , fvec(i_l) , el%mag
!    end select 
!   end do
!   ! debug ----
! 
!   ! === Newton's iterations ===
!   if  ( it .ge. 2000 ) then
!   ! dclvec set in the previous part of the algorithm
! 
!   write(*,*) ' === Before Newton''s iterations === '
!   write(*,*) '      norm2(res)  : ' ,      norm2(fvec)
!   write(*,*) ' maxval(abs(res)) : ' , maxval(abs(fvec)) 
! 
!   allocate(Amat_newton(nll,nll)) ! ; Amat_newton = Amat ! since dgsev destroys the input mat
!   allocate(d_al(nll))            ! ; d_al = 0.0_wp
!   allocate(ipiv(nll))            ! needed by dgesv
! 
!   res = norm2(fvec) ; res_old = 2.0_wp * res
!   newton_it = 1
!   do while( ( res .gt. newton_tol ) .and. &
!             ( newton_it .lt. newton_maxit ) .and. & 
!             ( res .lt. res_old ) )
! 
!     ! fill Amat_newton ----
!     do ik = 1 , nll
!       Amat_newton(:,ik) = Amat(:,ik) * dclvec(ik)
!     end do
! 
! !   ! debug ----
! !   write(*,*) ' Amat_newton : '
! !   do ii = 1 , nll
! !     write(*,*) Amat_newton(ii,:)
! !   end do
! !   ! debug ----
! 
!     ! Compute the increment
!     d_al = -fvec
!     call dgesv(nll,1,Amat_newton,nll,ipiv,d_al,nll,info) 
! 
!     if ( info .ne. 0 ) then
!       write(msg,*) 'error while solving linear system, Lapack DGESV error code ', info
!       call error(this_sub_name, this_mod_name, trim(msg))
!     end if
! 
!     ! Update the solution 
!     alvec = alvec + 0.15_wp * d_al
! 
!     do i_l = 1 , nll
! !     write(*,*) ' i_l , al , ma , re : ' , &
! !                  i_l , alvec(i_l) , ma_v(i_l) , re_v(i_l)
!       select type( el=>elems_ll(i_l)%p ) ; type is (t_liftlin)
!         call interp_aero_coeff ( airfoil_data,  el%csi_cen, el%i_airfoil, &
!                       (/alvec(i_l), ma_v(i_l), re_v(i_l)/) , aero_coeff , dclda , al0 )
!       end select
!       c_m(i_l,:)  = aero_coeff     ! overwrite ! 
!       clvec( i_l) = aero_coeff(1)
!       dclvec(i_l) = dclda
!     end do
! 
!     fvec = matmul(Amat,clvec) + bvec  ! residue 
! 
!     res_old = res
!     res = norm2(fvec)
!     
!     write(*,*) ' n.it , res : ' , newton_it , norm2(fvec)
! 
!     newton_it = newton_it + 1
! 
!   end do
! 
! ! write(*,*) ' alvec: '
! ! write(*,*)   alvec   
! 
!   a_v = alvec * pi/180.0_wp ! overwrite
! 
!   ! Loads computation ------------
!   ! Is it necessary to divide the loads by (cos(alpha))^2 ?
!   do i_l = 1,size(elems_ll)
!    select type(el => elems_ll(i_l)%p) ; type is(t_liftlin)
!     ! singualrity intensity
!     el%mag    = -0.5_wp * u_v(i_l) * clvec(i_l) * el%chord
! 
!     ! avg delta_p = \vec{F}.\vec{n} / A = ( L*cos(al)+D*sin(al) ) / A
!     el%pres   = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
!                ( c_m(i_l,1) * cos(a_v(i_l)) +  c_m(i_l,2) * sin(a_v(i_l)) )
!     ! elementary force = p*n + tangential contribution from L,D
!     el%dforce = ( el%nor * el%pres + &
!                   el%tang_cen * &
!                   0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * ( &
!                  -c_m(i_l,1) * sin(a_v(i_l)) + c_m(i_l,2) * cos(a_v(i_l)) &
!                  ) ) * el%area
!     ! elementary moment = 0.5 * rho * v^2 * A * c * cm, 
!     ! - around bnorm_cen (always? TODO: check)
!     ! - referred to the ref.point of the elem,
!     !   ( here, cen of the elem = cen of the liftlin (for liftlin elems) )
!     el%dmom = 0.5_wp * sim_param%rho_inf * u_v(i_l)**2.0_wp * &
!                    el%chord * el%area * c_m(i_l,3)
!     end select
!   end do
! 
!   write(*,*) ' alvec , ll%mag , u , cl ,  cd ,  cm ,  pres , dforce' 
!   do i_l = 1 , nll
!     write(*,*) alvec(i_l) , elems_ll(i_l)%p%mag , u_v(i_l) , &
!                c_m(i_l,:) , elems_ll(i_l)%p%pres , elems_ll(i_l)%p%dforce
!   end do
! 
! ! stop
!   end if
! 
! 
! 
!   ! === Newton's algorithm to assign u_rel.n = 0 ===
! 
! 
! ! Debug ----
! 
! 
! 
! ! do i_l = 1 , size(elems_ll)
! !   vel_ctr_rel(:,i_l) = vel_ctr(:,i_l) + uinf - elems_ll(i_l)%p%ub
! ! end do
! 
! ! Debug ----
! 
! ! debug ----
!   write(it_str,'(I4.4)') it 
!   open(unit=21,file='./Debug/liftlin_bc_it_'//trim(it_str)//'.dat')
!   do i_l = 1 , size(elems_ll)
!     select type( el => elems_ll(i_l)%p ) 
!     type is(t_liftlin)
!       write(21,*) el%cen , a_v(i_l) , c_m(i_l,:) , re_v(i_l) , ma_v(i_l) ! , &
! !         sum( vel_ctr_rel(:,i_l) * el%nor )
!     end select
!   end do
!   close(21)
! 
!   deallocate( vel_ctr )
! ! debug ----
! 
! ! newton ---
!   deallocate( Amat , bvec , fvec , clvec , dclvec )
! ! missing deallocate ...
! newton ---
! newton cleaning

  deallocate(dou_temp, vel_w)
  deallocate(a_v,c_m,u_v)

end subroutine solve_liftlin

!----------------------------------------------------------------------

subroutine calc_geo_data_liftlin(this, vert)
 class(t_liftlin), intent(inout) :: this
 real(wp), intent(in) :: vert(:,:)

 integer :: is, nsides
 real(wp):: nor(3), tanl(3)

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

  this%area = 0.5_wp * norm2(nor) * 4.0_wp/3.0_wp
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

! this%bnorm_cen = cross(this%tang_cen, this%nor)  ! old
  this%bnorm_cen = this%ver(:,2) - this%ver(:,1)
  this%bnorm_cen = this%bnorm_cen / norm2(this%bnorm_cen)

  this%chord = sum(this%edge_len((/2,4/)))*0.5_wp

  this%ctr_pt = this%cen + this%tang_cen * this%chord / 2.0_wp

  ! overwrite nor
  this%nor = cross( this%bnorm_cen , this%tang_cen )
  this%nor = this%nor / norm2(this%nor)

  ! *** 2019-02-27 ***
  ! ...before 2019-02-27 at the beginning of the routine
  ! local tangent unit vector as in PANAIR
  tanl = 0.5_wp * ( this%ver(:,nsides) + this%ver(:,1) ) - this%cen

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
