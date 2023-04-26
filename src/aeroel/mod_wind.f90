
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
!! Copyright (C) 2018-2022 Politecnico di Milano,
!!                           with support from A^3 from Airbus
!!                    and  Davide   Montagnani,
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
!!          Andrea Colli
!!=========================================================================

module mod_wind

use mod_param, only: &
  wp, pi

use mod_sim_param, only: &
  sim_param
!----------------------------------------------------------------------

implicit none

public :: variable_wind

private

!----------------------------------------------------------------------
contains
!----------------------------------------------------------------------

function variable_wind(pos, time) result(wind)

 real(wp), intent(in) :: pos(3)
 real(wp), intent(in) :: time

 real(wp) :: wind(3)

 real(wp) :: gust_origin(3), gust_front_direction(3), gust_front_speed, &
             gust_u_des, gust_perturb_direction(3), gust_gradient, gust_time

 real(wp) :: s

  wind = sim_param%u_inf

  if (sim_param%use_gust) then
    select case(trim(sim_param%GustType))
      case('AMC')
        gust_origin = sim_param%gust_origin
        gust_front_direction = sim_param%gust_front_direction/ &
                                            norm2(sim_param%gust_front_direction)
        gust_front_speed = sim_param%gust_front_speed
        gust_u_des = sim_param%gust_u_des
        gust_perturb_direction = sim_param%gust_perturb_direction/&
                                          norm2(sim_param%gust_perturb_direction)
        gust_gradient = sim_param%gust_gradient
        gust_time = sim_param%gust_time

        ! penetration distance
        ! distance from the gust front, negative for the gust approaching
        s = -sum((pos-(gust_origin+gust_front_speed*gust_front_direction*&
                                         (time-gust_time)))*gust_front_direction)

        if (s .ge. 0.0_wp .and. s .lt. 2.0_wp*gust_gradient) then
          wind = wind + gust_u_des/2*(1.0_wp-cos(pi*s/gust_gradient))
        end if

      case('linear') ! for testing
        gust_origin = sim_param%gust_origin
        gust_front_direction = sim_param%gust_front_direction/&
                                            norm2(sim_param%gust_front_direction)
        gust_front_speed = sim_param%gust_front_speed
        gust_u_des = sim_param%gust_u_des
        gust_perturb_direction = sim_param%gust_perturb_direction/&
                                          norm2(sim_param%gust_perturb_direction)
        gust_gradient = sim_param%gust_gradient
        gust_time = sim_param%gust_time

        s = sum((pos-(gust_origin+sim_param%u_inf*time))*sim_param%u_inf)/&
                                                           norm2(sim_param%u_inf)

        wind = wind + real(pos(1),wp)*(/0.0_wp, 0.0_wp, 0.1_wp/)!s*gust_u_ds/gust_gradient*(/0.0, 0.0, 0.1/)
      case default
    end select
  end if

end function variable_wind

!----------------------------------------------------------------------


!----------------------------------------------------------------------

end module mod_wind
