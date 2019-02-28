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


!> Module to treat several parts of the pressure integral equation

module mod_pressure_equation


use mod_param, only: &
  wp, nl, max_char_len, pi

use mod_sim_param, only: &
  t_sim_param, sim_param

use mod_handling, only: &
  error, warning, info, printout, dust_time, t_realtime

use mod_aeroel, only: &
  t_impl_elem_p

use mod_surfpan, only: &
  t_surfpan

use mod_linsys_vars, only: t_linsys

use mod_geometry, only: &
  t_geo

!----------------------------------------------------------------------

implicit none

public :: dump_linsys_pres, press_normvel_der

private

character(len=*), parameter :: this_mod_name = 'mod_pressure_equation'

!----------------------------------------------------------------------

contains

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
