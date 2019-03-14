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


!> Module to treat several parts of the pressure integral equation

module mod_pressure_equation


use mod_param, only: &
  wp, nl, max_char_len, pi

use mod_sim_param, only: &
  t_sim_param, sim_param

use mod_math, only: &
  cross

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_aeroel, only: &
  t_impl_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_linsys_vars, only: t_linsys

use mod_geometry, only: &
  t_geo

use mod_wake, only: &
  t_wake

!----------------------------------------------------------------------

implicit none

public :: dump_linsys_pres, press_normvel_der, initialize_pressure_sys, &
          assemble_pressure_sys, solve_pressure_sys

private

character(len=*), parameter :: this_mod_name = 'mod_pressure_equation'

!----------------------------------------------------------------------

contains

!----------------------------------------------------------------------

!> Initialize the linear system of the pressure
!!
subroutine  initialize_pressure_sys(linsys, geo, elems)
 type(t_linsys), intent(inout), target :: linsys
 type(t_geo), intent(in) :: geo
 type(t_impl_elem_p), intent(inout) :: elems(:)
 integer :: ie

  !Allocate matrix and rhs for pressure integral equation
  allocate( linsys%idSurfPan(geo%nSurfPan) )
  linsys%idSurfPan    = geo%idSurfPan
  allocate( linsys%idSurfPanG2L(geo%nSurfPan) )
  linsys%idSurfPanG2L = geo%idSurfPanG2L
  allocate( linsys%A_pres(geo%nSurfPan,geo%nSurfPan) )
  linsys%A_pres = 0.0_wp
  allocate( linsys%b_pres(geo%nSurfPan) )
  linsys%b_pres = 0.0_wp
  allocate( linsys%res_pres(geo%nSurfPan) )
  linsys%res_pres = 0.0_wp
  allocate( linsys%b_matrix_pres_debug(geo%nSurfPan,geo%nSurfPan) )
  linsys%b_matrix_pres_debug = 0.0_wp
  allocate( linsys%b_static_pres(geo%nstatic_SurfPan,geo%nstatic_SurfPan) ) 
  linsys%b_static_pres = linsys%b_static( & 
         geo%idSurfPan(1:geo%nstatic_SurfPan), &
         geo%idSurfPan(1:geo%nstatic_SurfPan) )

  do ie = 1 , geo%nSurfpan
    select type( el => elems(geo%idSurfPan(ie))%p ) ; class is(t_surfpan)
      el%pres_sol => linsys%res_pres(ie) 
    end select
  end do

end subroutine initialize_pressure_sys


!----------------------------------------------------------------------

!> Assemble the pressure equation system
subroutine assemble_pressure_sys(linsys, geo, elems, wake)
 type(t_linsys), intent(inout) :: linsys
 type(t_geo), intent(in) :: geo
 type(t_impl_elem_p), intent(in) :: elems(:)
 type(t_wake), intent(in) ::  wake

 integer :: ie, ip, iw, p1, p2, inext, is
 integer :: ipp(4) , iww(4), ntot
 real(wp) :: elcen(3), dist(3), dist2(3)
 real(wp) :: Pinf, rhoinf, uinf(3)

  ! Free-stream conditions
  uinf   = sim_param%u_inf
  Pinf   = sim_param%P_inf
  rhoinf = sim_param%rho_inf
  ntot = linsys%rank


  ! Pressure integral equation +++++++++++++++++++++++++++++++++++++++++

  ! Slicing --------------------
  linsys%A_pres = linsys%A( geo%idSurfPan , geo%idSurfPan )
  ! Remove Kutta condtion ------
  do ie = 1 , geo%nSurfpan
    select type( el => elems(geo%idSurfPan(ie))%p ) ; class is(t_surfpan)
      call el%correct_pressure_kutta( &
            (/wake%pan_p, wake%rin_p/), wake%pan_gen_elems_id, linsys,uinf,ie,1,ntot)
    end select
  end do
  
  ! rhs: ....

  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
  ! Assemble the RHS of the linear system for the Bernoulli polynomial     !
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
  ! RHS = B_infty +                                 (a) far field
  !       + \oint_{Sb} { G ( du/dt - ... ) } +      (b) time-dependent
  !       +  \int_{V } { DG . \omega x U } +        (c) rotational
  !       +  \int_{Sb} { viscous terms }            (d) viscous
   
  ! RHS assembled as the sum of:
  ! (b) time-dependent source contribution: computed and assembled in the
  !     %build_row routines above
  ! do nothing here
  ! (c) rotational effects (up to now, ignoring ring wakes)
!$omp parallel do private(ie, iw, dist, elcen) schedule(dynamic,32) 
  do ie = 1 , geo%nSurfpan
    elcen = elems( geo%idSurfpan(ie) )%p%cen
    ! (c.1) particles ( part_p )
    if ( allocated( wake%part_p ) ) then
      do iw = 1 , size( wake%part_p )
        if  ( .not. ( wake%part_p(iw)%p%free ) ) then

          dist = wake%part_p(iw)%p%cen - elcen
          linsys%b_pres(ie) = linsys%b_pres(ie) + & 
            sum( dist * cross( wake%part_p(iw)%p%dir , wake%part_p(iw)%p%vel ) ) * &
                 wake%part_p(iw)%p%mag / ( (sqrt(sum(dist**2)+sim_param%VortexRad**2))**3.0_wp ) 
! singular  >>>  wake%part_p(iw)%p%mag / ( norm2(dist)**3.0_wp ) 
! Rosenhead >>>  wake%part_p(iw)%p%mag / ( (sqrt(sum(dist**2)+sim_param%VortexRad**2))**3.0_wp ) 
                 !EXPERIMENTAL: adding vortex rosenhead regularization
        end if
      end do 
    end if
  enddo
!$omp end parallel do

!$omp parallel do private(ie, iw, dist, dist2, p1, p2, ipp, iww,ip, is, inext, elcen) schedule(dynamic, 64) 
  do ie = 1 , geo%nSurfpan
    elcen = elems( geo%idSurfpan(ie) )%p%cen
    ! (c.2) line elements ( end_vorts )
    do iw = 1 , wake%n_pan_stripes
      if( associated( wake%end_vorts(iw)%mag ) ) then
        dist  = wake%end_vorts(iw)%ver(:,1) - elcen 
        dist2 = wake%end_vorts(iw)%ver(:,2) - elcen
        linsys%b_pres(ie) = linsys%b_pres(ie) - &
          0.5_wp * wake%end_vorts(iw)%mag * sum( wake%end_vorts(iw)%edge_vec * &
               ( cross(dist , wake%end_vorts(iw)%ver_vel(:,1) ) /(norm2(dist )**3.0_wp) + &
                 cross(dist2, wake%end_vorts(iw)%ver_vel(:,2) ) /(norm2(dist2)**3.0_wp) ) )  
      end if
    end do
    ! (c.3) constant surface doublets = vortex rings ( wake_panels )
    do ip = 1 , wake%pan_wake_len
      do iw = 1 , wake%n_pan_stripes
 
        p1 = wake%i_start_points(1,iw) 
        p2 = wake%i_start_points(2,iw)
 
        ipp = (/ ip , ip , ip+1, ip+1 /)
        iww = (/ p1 , p2 , p2  , p1   /)
 
        do is = 1 , wake%wake_panels(iw,ip)%n_ver ! do is = 1 , 4 
          
          inext = mod(is,wake%wake_panels(iw,ip)%n_ver) + 1 
          dist  = wake%wake_panels(iw,ip)%ver(:,is   ) - elcen 
          dist2 = wake%wake_panels(iw,ip)%ver(:,inext) - elcen
          
          linsys%b_pres(ie) = linsys%b_pres(ie) - &
            0.5_wp * wake%wake_panels(iw,ip)%mag * sum( wake%wake_panels(iw,ip)%edge_vec(:,is) * &
                 ( cross(dist , wake%pan_w_vel(:,iww(is   ),ipp(is   )) ) /(norm2(dist )**3.0_wp) + &
                   cross(dist2, wake%pan_w_vel(:,iww(inext),ipp(inext)) ) /(norm2(dist2)**3.0_wp) ) )   
!           ! old ----
!           0.5_wp * wake%end_vorts(iw)%mag * sum( wake%wake_panels(iw2,iw)%edge_vec(:,is) * &
!                ( cross(dist , wake%wake_panels(iw2,iw)%ver_vel(:,is   ) ) /(norm2(dist )**3.0_wp) + &
!                  cross(dist2, wake%wake_panels(iw2,iw)%ver_vel(:,inext) ) /(norm2(dist2)**3.0_wp) ) )  
!                ( cross(dist , wake%pan_w_vel(:,iw2  , ) ) /(norm2(dist )**3.0_wp) + &
!                  cross(dist2, wake%pan_w_vel(:,iw2+1, ) ) /(norm2(dist2)**3.0_wp) ) )  
!           ! old ----
          
        end do 
      end do
    end do
    ! (c.4) rings from actuator disks
    ! TODO: ...

    ! (d) viscous effects
    ! TODO: to be implemented

  end do
!$omp end parallel do

  ! (a) far field contribution
  linsys%b_pres = linsys%b_pres + &
                    4*pi * ( Pinf + 0.5_wp * rhoinf * norm2(uinf) ** 2.0_wp ) ! H_inf
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !
  ! END Assemble the RHS of the linear system for the Bernoulli polynomial !
  ! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ !

  ! Pressure integral equation +++++++++++++++++++++++++++++++++++++++++
end subroutine assemble_pressure_sys

!----------------------------------------------------------------------

!> Solve the linear system for the pressure
!!
!! At the moment the linear system is solved directly by means of a call to 
!! a standard lapack solver.
subroutine solve_pressure_sys(linsys)
 type(t_linsys), intent(inout) :: linsys

 real(wp) , allocatable :: A_tmp(:,:)
 integer, allocatable :: IPIV(:)
 integer              :: INFO
 integer              :: ARNK
 character(len=max_char_len) :: msg
 character(len=*), parameter :: this_sub_name = 'solve_linsys_pressure'
 !TODO: cleanup output vars and calls
 !real(KIND=8) :: elapsed_time ! time_ini , time_fin
 !integer :: count1 , count_rate1 , count_max1 
 !integer :: count2 , count_rate2 , count_max2 

  ARNK = size(linsys%A_pres,1)
  if ( size(linsys%A_pres,2) .ne. ARNK ) then
    write(msg,*) ' matrix of the pressure integral equation A_pres is not square'
    call error(this_sub_name, this_mod_name, trim(msg))
  end if

  allocate(A_tmp(ARNK,ARNK) ) ; A_tmp = linsys%A_pres
  allocate(IPIV(ARNK))

  linsys%res_pres = linsys%b_pres 

  !call system_clock(count1,count_rate1,count_max1)
  call dgesv(ARNK,1,A_tmp,ARNK,IPIV, &
             linsys%res_pres ,ARNK,INFO) 
  !call system_clock(count2,count_rate2,count_max2)

  if ( INFO .ne. 0 ) then
    write(msg,*) 'error while solving linear system, Lapack DGESV error code ',&
                  INFO
    call error(this_sub_name, this_mod_name, trim(msg))
  end if

  deallocate(A_tmp,IPIV) 

end subroutine solve_pressure_sys

!----------------------------------------------------------------------

!> Compute the normal velocity derivative
!!
!! compute the time derivative of the normal component of the velocity on 
!!  surfpan to be used in the source rhs of the Bernoulli integral equation.
!! surf_vel_SurfPan_old should be saved at the end of the time step 
subroutine press_normvel_der(geo, elems, surf_vel_SurfPan_old)
 type(t_geo), intent(in) :: geo
 type(t_impl_elem_p), intent(inout) :: elems(:)
 real(wp), intent(in) :: surf_vel_SurfPan_old(:,:)

 integer :: i_el, i_e
 real(wp) :: GradS_Un(3), DivS_U
 character(len=*), parameter :: this_sub_name = 'press_normvel_der'

  do i_el = 1 , geo%nSurfPan

    select type ( el => elems(geo%idSurfPan(i_el))%p ) ; class is ( t_surfpan )

      el%dUn_dt = sum( el%nor * ( el%ub - & 
             surf_vel_SurfPan_old( i_el , : ) ) ) / sim_param%dt
!            surf_vel_SurfPan_old( geo%idSurfPanG2L(i_el) , : ) ) ) / sim_param%dt ! <<< mod-2018-12-21

      ! Compute GradS_Un
      GradS_Un = 0.0_wp
      do i_e = 1 , el%n_ver
        if ( associated(el%neigh(i_e)%p) ) then !  .and. &
          select type(el_neigh=>el%neigh(i_e)%p) ; class is (t_surfpan)
            GradS_Un = GradS_Un + &
              matmul( geo%refs( geo%components(elems(i_el)%p%comp_id)%ref_id )%R_g , &
                el%pot_vel_stencil(:,i_e) * ( &
                  sum(el%nor* (el_neigh%surf_vel - el%surf_vel) ) ) )
          end select
        else
!         select type(el_neigh=>el%neigh(i_e)%p) ; class is (t_surfpan)
            GradS_Un = GradS_Un + &
              matmul( geo%refs( geo%components(elems(i_el)%p%comp_id)%ref_id )%R_g , &
               el%pot_vel_stencil(:,i_e) * ( &
                  sum(el%nor* ( - 2.0_wp * el%surf_vel) ) ) )
!         end select
        end if
      end do
!     GradS_Un = GradS_Un - el%nor * sum(el%nor*GradS_Un) ! tangential projection

      ! Compute DivS_U
      DivS_U = 0.0_wp
      do i_e = 1 , el%n_ver
        if ( associated(el%neigh(i_e)%p) ) then !  .and. &
          select type(el_neigh=>el%neigh(i_e)%p) ; class is (t_surfpan)
            DivS_U = DivS_U + &
               sum( &
                 matmul( geo%refs( geo%components(elems(i_el)%p%comp_id)%ref_id )%R_g ,   &
                                                            el%pot_vel_stencil(:,i_e) ) * &
                    ( el_neigh%surf_vel - el%surf_vel )   )
          end select
        else  
!         select type(el_neigh=>el%neigh(i_e)%p) ; class is (t_surfpan)
            DivS_U = DivS_U + &
              sum( &
                matmul( geo%refs( geo%components(elems(i_el)%p%comp_id)%ref_id )%R_g ,   & 
                                                           el%pot_vel_stencil(:,i_e) ) * &
                   ( - 2.0_wp * el%surf_vel ) )
!         end select
        end if
      end do  

      ! Compute "source intensity" of Bernoulli equations
      el%bernoulli_source = + el%dUn_dt & !    n . DU/Dt 
         - sum( GradS_Un * ( el%ub ))   & !  - GradS_Un . el%ub 
         + DivS_U * sum(el%ub*el%nor)     !  + Un * Div_S U

    end select
  end do

end subroutine press_normvel_der


!> Dump the matrix and rhs of the pressure linear system to file
!!
!! TODO: consider moving these functionalities to the i/o modules
subroutine dump_linsys_pres(linsys , filen_A , filen_b , filen_b_debug )
 type(t_linsys), intent(in) :: linsys
 character(len=*) , intent(in) :: filen_A , filen_b
 character(len=*) , optional , intent(in) :: filen_b_debug

 integer :: fid
 integer :: i1 


 fid = 23
 open(unit=fid, file=trim(adjustl(filen_A)) )
 do i1 = 1 , size(linsys%A_pres,1)
  write(fid,*) linsys%A_pres(i1,:)
 end do
 close(fid)
 
 fid = 24
 open(unit=fid, file=trim(adjustl(filen_b)) )
 do i1 = 1 , size(linsys%b_pres,1)
  write(fid,*) linsys%b_pres(i1)
 end do
 close(fid)

 if ( present( filen_b_debug ) ) then
   fid = 24
   open(unit=fid, file=trim(adjustl(filen_b_debug)) )
   do i1 = 1 , size(linsys%b_matrix_pres_debug,1)
    write(fid,*) linsys%b_matrix_pres_debug(i1,:)
   end do
   close(fid)
 end if

end subroutine dump_linsys_pres

!----------------------------------------------------------------------

end module mod_pressure_equation
